#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 04 15:30:35 2016



Please feel free to contact me for any question.
--
Zewei Song
University of Minnesota
Dept. Plant Pathology
songzewei@outlook.com
www.songzewei.org
"""

def main():
    function_list = []
    with open('function_list.txt', 'rU') as f:
        for line in f:
            line = line.strip('\n')
            function_list.append(line)
    
    
    with open('temp_main.txt', 'wb') as f:
        for function in function_list:
            f.write('if args.{0}:\n'.format(function))
            f.write('\timport {0} as function\n'.format(function))
            f.write('\tfunction.main(sub_args)\n')
            f.write('\n')

if __name__ == '__main__':
    main()
