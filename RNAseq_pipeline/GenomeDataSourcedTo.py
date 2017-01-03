#!/bin/bash/python
#encoding:utf-8

import sys

class GenomeDataSourcedTo():
    """
    this class would init with one input and two outputs names specifying output list of urls and output list of hash code
    """
    def __init__(self, o_data_url, o_md5sum):
        self._o_data_url = o_data_url
        self._o_md5sum = o_md5sum

    def from_file(self, filename):
        fastq_url_list = open(self._o_data_url, 'w')
        md5sum_code_list = open(self._o_md5sum, 'w')
        opened_inputfile = open(filename, 'r')

        for line in opened_inputfile:
            if 'md5sum' in line:
                md5sum_code_line = line.split(' ')[1]
		separator_idx = md5sum_code_line.find('&')
		md5sum_code = md5sum_code_line[:separator_idx] + '\n'
                md5sum_code_list.write(md5sum_code)
            if 'href' in line and 'gz' in line:
                fastq_url = line.split('"')[1] + "\n"
                fastq_url_list.write(fastq_url)

        fastq_url_list.close()
        md5sum_code_list.close()
        opened_inputfile.close()

    def from_url(self, url):
        pass


if __name__ == '__main__':
    genomeDataSource = GenomeDataSourcedTo(sys.argv[2], sys.argv[3])
    genomeDataSource.from_file(sys.argv[1])

