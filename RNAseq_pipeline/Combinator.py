#/bin/bash/python
# encoding:utf-8


import sys


class Combinator(object):
    def __init__(self, output, input_list):
        self.filenames = []
        with open(input_list, 'r') as filenames:
            for filename in filenames:
                self.filenames.append(filename.rstrip())

        self.output = output
        self.content = []

    def _file2list(self, filename, deliminator='\t'):
        """
        @ params: str, filename
        @ params: str, deliminator, default set to '\t'
        @ returns: list, content
        """
        ctn_list = []
        # content list
        opened_file = open(filename, 'r')
        for line in opened_file:
            linelist = line.rstrip().split(deliminator)
            # split() method takes a str and outputs a new list object
            # that's why there is no need to make a copy
            ctn_list.append(linelist)
        
        opened_file.close()

        return ctn_list

    def _find_special_char(self, string, char, first_any_last):
        """
        @ params: str char special characters
        @ params: str first_any_last indicator of which one you are looking for
        @ returns: int idx of the specific occurrence of that special char in given string
        """
        pass

        
    def combine_featureCounts(self, skip_lines=1, last_column_only=True, extract_sample_name=True):
        """
        @ params: int, skip_lines, how many lines do you want to skip; default is 1 (header);
        @ params: bool specify the range to extract content; default is True;
        @ params: bool tell if you want to take the base of the sample name
        @ params: *args files to be combined
        @ returns: a list of combined content
        """
        if len(self.filenames) < 2:
            raise Exception('Error: not enough featureCounts output found.')
        
        ctn_combined = []

        for i in xrange(len(self.filenames)):
            ctn_list = self._file2list(self.filenames[i])
            # remove skipped lines
            for _ in xrange(skip_lines):
                ctn_list.pop(0)
            # take the base of sample name field
            samplename = ctn_list[0][-1]  # first line last element
            for j in xrange(len(samplename)-1, -1, -1):
                if samplename[j] == '/':
                    basename = samplename[j+1:]
                    break
            ctn_list[0][-1] = basename
            # do we need other columns at all ?
            if i > 0 and last_column_only:
                for line_idx, linelist in enumerate(ctn_list):
		    # make the element a list even there's only one element
		    # so their dimensions are in a uniform
                    ctn_list[line_idx] = [linelist[-1]]
                ctn_combined.append(ctn_list)
            else:
                ctn_combined.append(ctn_list)
            # re-organize the content
        total_samples = len(ctn_combined)
        total_lines = len(ctn_combined[0])
        for line_idx in xrange(total_lines):
            newline = []
            for sample_idx in xrange(total_samples):
                newline += ctn_combined[sample_idx][line_idx]
            self.content.append(newline)

        return None
    
    def write_out(self, deliminator='\t'):
        opened_output = open(self.output, 'w')
        if len(self.content) == 0:
            raise Exception("No content stored.")
        for line in self.content:
            newline = ''
            for element in line:
                newline += element + deliminator
            newline += '\n'
            opened_output.write(newline)
        opened_output.close()



if __name__ == '__main__':
    c = Combinator(sys.argv[2], sys.argv[1])
    c.combine_featureCounts()
    c.write_out()
