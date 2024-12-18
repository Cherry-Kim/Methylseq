import string,sys,os
  
#https://github.com/brentp/bwa-meth
#https://github.com/huishenlab/dupsifter/blob/main/example/README.md
#-------------------------
## Installation
#wget https://github.com/brentp/bwa-meth/archive/master.zip
#unzip master.zip
#cd bwa-meth-master/
#sudo python setup.py install
#-------------------------

def STEP0_RefIndex(REF):
    #Creating the Reference Index
    os.system('/BIO1/bwa-meth-master/bwameth.py index '+REF)

def STEP1(REF,sample):
    #https://github.com/huishenlab/dupsifter/blob/main/example/README.md
    #Aligning Reads to the Reference
    #os.system('/BIO1/bwa-meth-master/bwameth.py -t 96 --reference '+REF+' '+sample+'_1.fastq.gz '+sample+'_2.fastq.gz | samtools view -b - | samtools sort -@ 96 -o '+sample+'.bam && samtools index '+sample+'.bam')
    os.system('samtools flagstat '+sample+'.bam > '+sample+'_flagstat.txt')

def STEP2(sample):
    if not os.path.isdir('TEMP_PICARD'):
        os.mkdir('TEMP_PICARD')
    os.system('java -Xmx16g -jar /BIO1/picard.jar MarkDuplicates TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT I='+sample+'.bam O='+sample+'.dedup.bam M='+sample+'.duplicate_metrics REMOVE_DUPLICATES=true AS=true CREATE_INDEX=true')


def STEP3(REF, sample):
    # Create a pileup VCF of DNA methylation and genetic information
    os.system('/BIO1/biscuit/biscuit pileup -@ 96 -o '+sample+'_pileup.vcf '+REF+' '+sample+'.dedup.bam && bgzip -@ 96 '+sample+'.vcf && tabix -p vcf '+sample+'.vcf.gz')
'''
    os.system('/BIO1/biscuit/biscuit pileup -@ 96 -o bwameth_pileup.vcf '+REF+' bwameth.sorted.markdup.bam')
    os.system('bgzip -@ 96 bwameth_pileup.vcf')
    os.system('tabix -p vcf bwameth_pileup.vcf.gz')
'''

    #Methylation percentage
    #os.system('/BIO1/biscuit/biscuit vcf2bed -c '+sample+'.vcf.gz > '+sample+'_beta_m_u.bed')

def main():
    REF = '/BIO1/REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    #STEP0_RefIndex(REF)
    file_list = os.listdir('./')
    file_ls = [file for file in file_list if file.endswith('_1.fastq.gz')]
    file_ls.sort()
    for i in file_ls:
        sample = i.split('_1.fastq.gz')[0]
        print(sample)
        #STEP1(REF,sample)
        STEP2(sample)
        #STEP3(REF, sample)
main()
