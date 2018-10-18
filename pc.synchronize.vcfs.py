
from datetime import datetime
from glob import glob
from argparse import ArgumentParser
import requests
import os
import logging
import subprocess

logging.basicConfig(level=logging.INFO, format='%(message)s', handlers=[logging.FileHandler('synchronize.log'), logging.StreamHandler()])

'''
	pc.synchronize.vcfs.py
	October 4th, 2018
	Dennis Kao, Center for Computational Medicine

	Description:
		Script to match and upload vcf files on our local file system to users on PhenomeCentral.

	Implementation details:
		requests: 
			gets a list of user data from PhenomeCentral in JSON format
		matching algorithm:
			exact_match: user's external_id is identical to a vcf file's name minus the extension (123_456R matches 123_456R.vcf)
			prefix_match: user's external_id is the prefix to a vcf (123_456 matches 123_456_replacement_2017.06.12.vcf)
		uploading vcfs:
			run a script created by PhenomeCentral team for each match of user id to vcf file

	Output:
		Stats on vcfs on filesystem, vcf-user matching result, vcfs on PhenomeCentral

	System requirements:
		Python 3.5
		UNIX file system
'''

def get_user_data(pc_credentials):

	def parse_auth(pc_credentials):
		with open(pc_credentials) as f:
			line = f.readline().strip('\n').split(':')
			username = line[0]
			password = line[1]
			return (username, password)

	user_data_url = "https://phenomecentral.org/get/PhenomeCentral/patient_list_script?action=list"

	with requests.Session() as session:
		session.auth = parse_auth(pc_credentials)
		user_data = session.get(user_data_url).json()["data"]
		
	return user_data

def scan_filesystem_vcf(results_dir, exhaustive_search=False):

	'''
	TODO: Use file creation date when deciding between 2 identically named files instead of modification date -- os.path.getmtime()

	Note on exhaustive_search:
	Our file structure is extremely consistent. 99.9999% of the time, vcfs that you want to find are detected by the regex defined below.
	VCF's which fall outside of this condition are mostly vcfs from an older pipeline version, or vcfs uploaded from collaborates and not created from our pipelines.
	'''

	vcf_paths = []
	vcf_dict = {}
	duplicates = {}

	if not exhaustive_search:		
		vcf_paths = glob('{}/*/*/*.vcf'.format(results_dir)) # equivalent to: ls results_dir/*/*/*.vcf
	else:
		for dirpath, dirname, files in os.walk(results_dir):
			for filename in files:
				if filename.endswith('.vcf'):
					vcf_paths.append(os.path.join(dirpath, filename))

	# exclude some directories known to contain junk
	vcf_paths = [ p for p in vcf_paths if "/misc/" not in p and "/miscx/" not in p and "/calx/" not in p and "/database/" not in p ]

	# detect and filter out identically named vcf files
	for path in vcf_paths:
		filename = os.path.basename(path)
		if filename in duplicates.keys():
			duplicates[filename].append(path)
		elif filename in vcf_dict.keys():
			duplicates[filename] = []
			duplicates[filename].append(path)
			duplicates[filename].append(vcf_dict[filename])
			vcf_dict.pop(filename)
		else:
			vcf_dict[filename] = path

	if duplicates:
		logging.warning('[DUPLICATES] Identically named vcfs detected. Attempting to resolve based on earliest modification date.')
		for filename, paths in duplicates.items():
			paths = sorted(paths, key=os.path.getmtime)
			if os.path.getmtime(paths[-1]) == os.path.getmtime(paths[-2]):
				logging.warning('Could not resolve duplicates of {} based on earliest modification date. Script won\'t upload any of these files: {}'.format(filename, str(paths)))
			else:
				vcf_dict[filename] = paths[-1]
				logging.warning('Using {} out of all files: {}'.format(vcf_dict[filename], str(paths)))

	return vcf_dict

def match_vcf2user(vcfs, no_vcf_users):

	def print_matched(message, users, matched):
		header = '\t'.join(['ID', 'EXTERNAL_ID', 'VCF'])
		external_id = { user["id"]:user["external_id"] for user in users if user["id"] in matched.keys() }
		logging.info(message)
		logging.info(header)
		for id, paths in matched.items():
			logging.info('\t'.join([id, external_id[id], str(paths)]))

	def remove_matched(vcf_dict, users, matched_dict):
		return {vcf:path for vcf, path in vcf_dict.items() if path not in matched_dict.values()}, [user for user in users if user["id"] not in matched_dict.keys()]

	def match_exact(vcf_dict, users):
		def vcf_extension(id):
			return u'{}.vcf'.format(id)
		return { user["id"]:vcf_dict[vcf_extension(user["external_id"])] for user in users if vcf_extension(user["external_id"]) in vcf_dict.keys() }

	def match_prefix(vcf_dict, users):
		prefix_match = {}
		multiple_prefix_match_hits = {}

		for user in users:
			external_id = user["external_id"]
			pc_id = user["id"]
			for vcf, filepath in vcf_dict.items():
				if vcf.startswith(external_id) or vcf.startswith(external_id.strip()) or vcf.startswith(external_id.strip('-')) or vcf.startswith(external_id.strip('_')):
					if pc_id in multiple_prefix_match_hits.keys():
						multiple_prefix_match_hits[pc_id].append(filepath)
					elif pc_id not in prefix_match.keys():
						prefix_match[pc_id] = filepath
					else:
						multiple_prefix_match_hits[pc_id] = []
						multiple_prefix_match_hits[pc_id].append(filepath)
						multiple_prefix_match_hits[pc_id].append(prefix_match[pc_id])
						prefix_match.pop(pc_id)

		if multiple_prefix_match_hits:
			print_matched("{}: These files will not be uploaded. External_id is the prefix to multiple VCFs:".format(match_prefix.__name__), users, multiple_prefix_match_hits)

		return prefix_match

	matched = {}
	has_external_id = [ user for user in no_vcf_users if user["external_id"] is not None and not user["external_id"].isspace() and len(user["external_id"]) != 0 ]

	for matching_algo in [ match_exact, match_prefix ]:
		newly_matched = matching_algo(vcfs, has_external_id)

		if newly_matched:
			message = '{}: managed to match {} users'.format(matching_algo.__name__, str(len(newly_matched)))
			print_matched(message, has_external_id, newly_matched)
		else:
			logging.info('{}: did not manage to match any users to vcfs'.format(matching_algo.__name__))

		matched.update(newly_matched)
		vcfs, has_external_id = remove_matched(vcfs, has_external_id, matched)

	unmatched = [ u for u in no_vcf_users if u["id"] not in matched.keys() ] # consented, no vcf, and no matches with vcf on filesystem
	return matched, unmatched

def upload_to_pc(pc_credentials, matched):
	with open("failed_uploads.txt", "w") as failed:
		for id, path in matched.items():
			logging.info('Uploading {} to PhenomeCentral user {}'.format(path, id))
			cmd = 'python ~/cre/pc.upload_vcf.py {} {} {} phenomecentral.org'.format(id, path, pc_credentials)
			try:
				p = subprocess.check_output(cmd, shell=True, timeout=10)
				logging.info(p)
			except Exception as e:
				failed.write('{}\t{}\t{}\n'.format(id, path, str(e)))

def synchronize(pc_credentials, results_dir, upload_flag):

	logging.info("Scanning filesystem for vcf's at: {}".format(results_dir))
	vcfs = scan_filesystem_vcf(results_dir)
	total_vcfs = len(vcfs)

	logging.info("Getting user data from PhenomeCentral")
	user_data = get_user_data(pc_credentials)

	granted = [ u for u in user_data if u["genetic_consent"] == "granted" ]
	not_granted = [ u for u in user_data if u["genetic_consent"] == "not granted" ]
	no_vcf = [ u for u in granted if u["vcf"] == "none" ] # consented, but no vcf
	has_vcf = [ u for u in granted if u["vcf"] != "none" ] # consented and has vcf
	already_uploaded = [ u for u in has_vcf if u["vcf"] in vcfs.keys() ] # consented, has vcf, and vcf matches with one on filesystem

	logging.info("Matching users to vcfs ...")
	matched, unmatched = match_vcf2user(vcfs, no_vcf)

	logging.info('')
	logging.info("Filesystem statistics:")
	logging.info("Total vcfs: %d" % total_vcfs)
	logging.info("VCFs already uploaded to PhenomeCentral: %d" % len(already_uploaded))
	logging.info("Remaining VCFs: %d" % (total_vcfs-len(already_uploaded)))
	logging.info('')
	logging.info("PhenomeCentral user statistics:")
	logging.info("Total users:\t%d" % len(user_data))
	logging.info("Granted consent:\t%d" % len(granted))
	logging.info("Not granted:\t%d" % len(not_granted))
	logging.info("Granted consent and has a vcf uploaded:\t%d" % len(has_vcf))
	logging.info("Granted consent and no vcf uploaded:\t%d" % len(no_vcf))
	logging.info("Granted consent, no vcf uploaded, and has a matched vcf:\t%d" % len(matched))
	logging.info("Granted consent, no vcf uploaded, and no matched vcf:\t%d" % len(unmatched))
	logging.info('')

	if upload_flag:
		upload_to_pc(pc_credentials, matched)

if __name__ == '__main__':

	parser = ArgumentParser()
	parser.add_argument('-pc_credentials', help='Text file containing PhenomeCentral login credentials in the format: username:password', required=True)
	parser.add_argument('-results_dir', help='Path to the results directory containing vcf files. Our typical directory structure: results/10x/1031R/1031R_1234.vcf', required=True)
	parser.add_argument('--upload', help='Upload VCF file\'s identified through matching to PhenomeCentral', action='store_true')
	args = parser.parse_args()

	logging.info('started vcf synchronization on {}'.format(datetime.now().strftime("%Y-%m-%d_%H:%M:%S")))
	synchronize(args.pc_credentials, args.results_dir, args.upload)
	logging.info('synchronization finished on {}'.format(datetime.now().strftime("%Y-%m-%d_%H:%M:%S")))
