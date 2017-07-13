#!/usr/bin/env python

"""
Usage: $0 PATIENT_ID VCF CREDENTIALS_FILE

Upload the VCF file as an attachment for the given patient. For example:

  $0 P0000102 path/to/file.vcf path/to/credentials.txt

    Phenomecentral API documentation
    https://phenotips.org/DevGuide/PatientsRESTfulAPI1.3#read-eid201
    
    propotype by Orion Buske 2017-07-13
"""

from __future__ import with_statement

import sys
import os
import logging

from httplib import HTTPConnection
from urllib import urlencode
from base64 import encodestring
import xml.etree.ElementTree as ET

# Set globally after parsing arguments
CREDENTIALS = None

def load_credentials(filename):
    global CREDENTIALS
    logging.info('Using credentials in file: {0}'.format(filename))
    with open(filename) as ifp:
        CREDENTIALS = ifp.read().strip()

class Request(object):
    def __init__(self, method, resource, content_type='text/plain', **kwargs):
        auth = encodestring(CREDENTIALS).replace('\n', '')
        self._headers = {
            'Authorization': 'Basic {0}'.format(auth),
            'Content-Type': content_type,
            'Accept': 'application/xml',
            }
        self._method = method
        self._href = resource
        self._kwargs = kwargs
        self._conn = None

    def __enter__(self):
        self._conn = HTTPConnection('localhost', 8080)
        #self._conn.set_debuglevel(1)
        self._conn.request(self._method, self._href, headers=self._headers, **self._kwargs)
        return self._conn.getresponse()

    def __exit__(self, *args, **kwargs):
        if self._conn:
            self._conn.close()


def get_vcf_object(record_id):
    href = '/rest/wikis/xwiki/spaces/data/pages/{0}/objects/PhenoTips.VCF'.format(record_id)
    result = {}
    with Request('GET', href) as response:
        if response.status == 500:
            # No VCF object for patient
            logging.info('No PhenoTips.VCF object for record: {0}'.format(record_id))
            return result

        assert response.status == 200, \
            'ERROR getting PhenoTips.VCF object for record: {0}'.format(record_id)

        logging.info('Found PhenoTips.VCF object for record: {0}'.format(record_id))
        root = ET.fromstring(response.read())

    # Server is Python 2.6.6, so does not support namespaces argument
    # namespaces = {'xwiki': 'http://www.xwiki.org'}
    summary = root.find('{http://www.xwiki.org}objectSummary')

    for link in summary.findall('{http://www.xwiki.org}link'):
        if link.get('rel') == 'http://www.xwiki.org/rel/object':
            href = link.get('href')
            logging.info('Found VCF object for {0}: {1}'.format(record_id, href))
            result['href'] = href
            break

    result['headline'] = summary.findtext('{http://www.xwiki.org}headline', '').strip()
    return result

def set_vcf_properties(record_id, vcf_object, filename, reference_genome='GRCh37'):
    method = 'PUT'
    href = vcf_object.get('href')
    if not href:
        # Create VCF missing properties object
        method = 'POST'
        href = '/rest/wikis/xwiki/spaces/data/pages/{0}/objects'.format(record_id)

    filebase = os.path.basename(filename)
    params = {
        'className': 'PhenoTips.VCF',
        'property#filename': filebase,
        'property#reference_genome': reference_genome,
        }

    content_type = 'application/x-www-form-urlencoded'

    with Request(method, href, content_type=content_type, body=urlencode(params)) as response:
        assert response.status in [201, 202], \
            'Unexpected response ({0}) from setting VCF properties'.format(response.status)

        if response.status == 201:
            href = response.getheader('Location')
            logging.info('Created VCF object: {0}'.format(href))
            return href
        elif response.status == 202:
            logging.info('Updated VCF object: {0}'.format(href))
            return vcf_object

def delete_old_vcf_file(record_id):
    attachment_name = '{0}.vcf'.format(record_id)
    href = '/rest/wikis/xwiki/spaces/data/pages/{0}/attachments/{1}'.format(record_id, attachment_name)
    with Request('HEAD', href) as response:
        if response.status in [200, 500]:
            logging.info('Found linked old-style VCF: {0}'.format(attachment_name))
        else:
            return

    with Request('DELETE', href) as response:
        assert response.status in [204, 500], \
            'Unexpected response ({0}) from deleting VCF file'.format(response.status)
        logging.info('Deleted attachment: {0}'.format(attachment_name))

    with Request('DELETE', href) as response:
        assert response.status == 404, \
            'Failed to delete existing old-style VCF file'.format(response.status)

def upload_vcf_file(record_id, filename):
    attachment_name = os.path.basename(filename)
    href = '/rest/wikis/xwiki/spaces/data/pages/{0}/attachments/{1}'.format(record_id, attachment_name)
    with open(filename) as file:
        with Request('PUT', href, content_type='text/plain', body=file) as response:
            assert response.status in [201, 202], \
                'Unexpected response ({0}) from uploading VCF file'.format(response.status)
            if response.status == 201:
                logging.info('Created attachment: {0}'.format(attachment_name))
            elif response.status == 202:
                logging.info('Updated attachment: {0}'.format(attachment_name))

def script(record_id, vcf_filename, credentials_file):
    load_credentials(credentials_file)

    vcf_object = get_vcf_object(record_id)

    vcf_basename = os.path.basename(vcf_filename)
    if vcf_object and vcf_object.get('headline', '') == vcf_basename:
        # Filename already appears to be properly linked to record, skip
        logging.info('Record {0} already has linked VCF file: {1}'.format(record_id, vcf_basename))
    else:
        href = set_vcf_properties(record_id, vcf_object, vcf_filename)

        if href:
            upload_vcf_file(record_id, vcf_filename)

    ## TODO: remove once PC is clean
    ## delete_old_vcf_file(record_id)


def parse_args(args):
    from optparse import OptionParser
    usage = "usage: %prog [options] PATIENT_ID VCF CREDENTIALS_FILE"
    parser = OptionParser(usage=usage)
    (options, args) = parser.parse_args()

    if len(args) != 3:
        parser.error('Invalid number of arguments')

    return (options, args)

def main(args=sys.argv[1:]):
    options, args = parse_args(args)
    kwargs = dict(options.__dict__)
    logging.basicConfig(level=logging.DEBUG)
    script(*args, **kwargs)

if __name__ == '__main__':
    sys.exit(main())
