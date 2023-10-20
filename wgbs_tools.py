#https://github.com/nloyfer/wgbs_tools
#Genome configuration
#wgbstools init_genome hg19

import string,sys,os
import subprocess

def step1():
    region, co = [], 0
    sample = 'GSM6810004_CNVS-NORM-110000264-cfDNA-WGBS-Rep1.pat.gz'
    fp = open('region2.txt','r') #chr13:27928900-27929266
    for line in fp:
        if line.strip():
            region.append(line.strip())
    fp.close()

    with open('GSM6810004.txt','a') as output_file:
        for i in region:
            co += 1
            print("Region",co," ",i)

            command = "/BIO5/wgbs_tools/wgbstools vis "+sample+" -r "+i
            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            stdout, stderr = process.communicate()

            if 'Methylation' in stdout.decode():
                output_lines = [line for line in stdout.decode().split('\n') if 'Methylation average:' in line]
                output_file.write('\n'.join(output_lines)+'\n')
            else:
                output_file.write('No results\n')

