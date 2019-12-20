import os
import sys
import glob
import argparse
import pprint
import subprocess
from collections import OrderedDict

ALLOWED_INPUT_EXTENSIONS = OrderedDict({"*.fq": "fastq","*.fastq": "fastq","*.fq.gz":"fastq","*.fastq.gz":"fastq","*.cram":"cram","*.bam":"bam"})

parser = argparse.ArgumentParser()

parser.add_argument("-u", "--upload-dir", help="Directory containing uploaded files for samples", required=True)
parser.add_argument("-d", "--processing-dir", help="Directory to create pipeline folders in", required=True)
parser.add_argument("-s", "--sample-list", help="File with SampleIDs. One per line. SampleIDs must be in the format famID_sampleName", required=True)


args = parser.parse_args()

if args.upload_dir is not None:
    if os.path.exists(args.upload_dir) and os.path.isdir(args.upload_dir):
        UPLOAD_DIRECTORY = args.upload_dir
    else:
        print(args.upload_dir+" is not a valid directory")
        sys.exit(0)


if args.processing_dir is not None:
    if os.path.exists(args.processing_dir) and os.path.isdir(args.processing_dir):
        PROCESSING_DIRECTORY = args.processing_dir
        os.chdir(PROCESSING_DIRECTORY)
        SUBMIT_SCRIPTS_DIRECTORY = os.path.join(PROCESSING_DIRECTORY,"submit_scripts")
        SUBMIT_LOG_DIRECTORY = os.path.join(PROCESSING_DIRECTORY,"submit_logs") 
        for directory in (SUBMIT_SCRIPTS_DIRECTORY, SUBMIT_LOG_DIRECTORY):
            if not os.path.exists(directory):
                subprocess.run(['mkdir',directory])
    else:
        print(args.processing_dir+" is not a valid directory")
        sys.exit(0)

if args.sample_list is not None:
    if os.path.exists(args.sample_list) and os.path.isfile(args.sample_list):
        SAMPLE_FILE = args.sample_list
    else:
        print(args.sample_list + " is not a valid file")
        sys.exit(0)


def check_and_return_input_files(family_id, sample_name):
    
        input_files = []
        input_file_type = None
        sample_id = '_'.join([family_id, sample_name]) 
        for extension,filetype in ALLOWED_INPUT_EXTENSIONS.items(): 
            for dirpath,dirnames,filename in os.walk(UPLOAD_DIRECTORY):
                for subdirectory in dirnames:
                    if subdirectory == sample_id: ## upload folders must be named in sample_id format 
                        input_files = glob.glob(os.path.join(dirpath,subdirectory,extension))
                        if filetype == "fastq" and len(input_files) >0:
                            input_files = list(filter(lambda elem: '_R1_' in elem or '_R2_' in elem, input_files))
                            if len(input_files) == 0:
                                print("Fastq files for this sample do not have _R1_ or _R2_ in them")
                        if len(input_files) >0:
                            input_file_type = filetype
                            break
                if len(input_files) >0:
                    break
            if len(input_files) >0:
                break
        return (input_files,input_file_type)

def prepare_input_files(inputs):

    family_id, sample_name = inputs['sample_id'].split("_")
    job_name = inputs['sample_id'] + "_" + inputs['command']
    job_file = os.path.join(SUBMIT_SCRIPTS_DIRECTORY,job_name+".sh")
    input_files_directory  = os.path.join(PROCESSING_DIRECTORY,family_id,"input") 
    cram_files_directory = os.path.join(PROCESSING_DIRECTORY,family_id,"cram")

    with open(job_file,"w") as submit_fh:
        
        #print_headers here
        #submit_fh.write("#PBS -l walltime=20:00:00,nodes=1:ppn=1\n#PBS -joe \n#PBS -o {1}\n#PBS -d {0}\n#PBS -l vmem=30g,mem=30g\n\n".format(PROCESSING_DIRECTORY,SUBMIT_LOG_DIRECTORY))
        submit_fh.write("#PBS -l walltime=20:00:00,nodes=1:ppn=1\n#PBS -o {1}\n#PBS -d {0}\n#PBS -l vmem=30g,mem=30g\n\n".format(PROCESSING_DIRECTORY,SUBMIT_LOG_DIRECTORY))
        submit_fh.write("mkdir -p {0}\n".format(input_files_directory))
        if inputs['input_files'][0].endswith('.cram'):
            submit_fh.write("mkdir -p {0}\n".format(cram_files_directory))
 
        if inputs['input_file_type'] == 'fastq':
            inp_link = None
            for input_file in inputs['input_files']:
                if '_R1_' in input_file:
                    inp_link = family_id + "_" + sample_name + "_1.fq"
                if '_R2_' in input_file:
                    inp_link = family_id + "_" + sample_name + "_2.fq"
                if input_file.endswith(".gz"):
                    inp_link +=".gz"
                submit_fh.write("ln -s {} {}\n".format(input_file,os.path.join(input_files_directory,inp_link)))
        
        if inputs['input_file_type'] == 'bam':
            inp_link = inputs['sample_id'] + ".bam"
            submit_fh.write("ln -s {} {}\n".format(inputs['input_files'][0],os.path.join(input_files_directory,inp_link)))

        if inputs['input_file_type'] == 'cram':
            inp_link = inputs['sample_id'] + ".cram"
            submit_fh.write("ln -s {} {}\n\n".format(inputs['input_files'][0],os.path.join(cram_files_directory,inp_link)))
            submit_fh.write("module load java\n\n")
            submit_fh.write("cramtools fastq -Xmx10g -F {} --skip-md5-check -z -I {} -R /hpf/largeprojects/ccmbio/naumenko/tools/bcbio/genomes/Hsapiens/hg19/seq/hg19.fa\n\n".format(inputs['sample_id'],os.path.join(cram_files_directory,inp_link)))        
            submit_fh.write("mv {} {}\n".format(os.path.join(PROCESSING_DIRECTORY,inputs['sample_id']+"_1.fastq.gz"),os.path.join(input_files_directory,inputs['sample_id']+"_1.fq.gz")))
            submit_fh.write("mv {} {}\n".format(os.path.join(PROCESSING_DIRECTORY,inputs['sample_id']+"_2.fastq.gz"),os.path.join(input_files_directory,inputs['sample_id']+"_2.fq.gz")))
            
        submit_fh.close()
        return job_file

def prepare_config_files(inputs):
    
    family_id, sample_name = inputs['sample_id'].split("_")
    job_name = family_id + "_" + inputs['command']
    job_file = os.path.join(SUBMIT_SCRIPTS_DIRECTORY,job_name+".sh")
    
    with open(job_file,"w") as submit_fh:
        
        #print_headers here
        submit_fh.write("#PBS -l walltime=1:00:00,nodes=1:ppn=1\n#PBS -joe\n#PBS -o {1}\n#PBS -d {0}\n#PBS -l vmem=20g,mem=20g\n\n".format(PROCESSING_DIRECTORY,SUBMIT_LOG_DIRECTORY))
        submit_fh.write("~/cre/cre.prepare_bcbio_run.sh "+ str(family_id))
    submit_fh.close()
    return job_file

def submit_job(**params):
   
    family_id, sample_name = params['sample_id'].split("_") 
    job_to_submit = None
    if params['command'] == 'prepare':
        job_to_submit = prepare_input_files(params)
    elif params['command'] == 'config':
        job_to_submit = prepare_config_files(params)
    elif params['command'] == 'bcbio':
        job_to_submit = ' ~/cre/bcbio.pbs -v project='+ str(family_id)
    elif params['command'] == 'cre':
        job_to_submit = ' ~/cre/cre.sh -v family='+ str(family_id)+',cleanup=1'

    if job_to_submit is not None:

        if len(params['previous_job_ids']) == 0:
            output = subprocess.run(["qsub",job_to_submit], stdout=subprocess.PIPE)
        else:
            previous_job_id_str = ":".join(params['previous_job_ids'])
            output = subprocess.run(["qsub -W depend=afterok:"+previous_job_id_str+" "+job_to_submit],shell=True, stdout=subprocess.PIPE)
        job_id = output.stdout.decode('utf-8').rstrip()
    
        return str(job_id)


steps = ["prepare","config","bcbio","cre"]

families = {}
with open(SAMPLE_FILE,"r") as inp_fh:

    for line in inp_fh:
      
        line = line.rstrip()
        if '_' not in line:
            print(line + " not in correct format. sample_id must be in the format famid_samplename")
            sys.exit(0)

        family_id, sample_name = line.split("_")
        
        if family_id not in families:
            families[family_id] =[]
        families[family_id].append(sample_name)

inp_fh.close()

for family_id in families:

    print("Processing family "+ str(family_id))
    job_ids = []
    run_bcbio = 1
    for sample_name in families[family_id]:
        sample_id = family_id+"_"+sample_name
        input_files, input_file_type = check_and_return_input_files(family_id, sample_name)

        if len(input_files) == 0:
            print("No recognized input files for sample "+sample_id+ ". Not setting up this family")
            run_bcbio = 0
            continue
        if (input_file_type == 'bam' or input_file_type == 'cram') and len(input_files)>1:
            print("Multiple bam/cram files exist for sample "+sample_id+". Not setting up this family")
            run_bcbio = 0 
            continue
        job_ids.append(submit_job(sample_id=sample_id,previous_job_ids=[],command="prepare",input_files=input_files, input_file_type=input_file_type))
    
    if run_bcbio:
        for step in steps[1:]:
            job_id = submit_job(sample_id=sample_id, previous_job_ids=job_ids, command=step, input_files=input_files, input_file_type=input_file_type)
            job_ids = [job_id] 
