import string,sys,os
  
#https://huishenlab.github.io/biscuit/biscuitsifter/
#-------------------------
## Download Source Code and Compile
#git clone --recursive https://github.com/huishenlab/biscuit.git
#cd biscuit
#make
#-------------------------

def STEP1_RefIndex(REF,sample):
    #Creating the Reference Index
    #os.system('/BIO1/biscuit/biscuit index '+REF)

    #Aligning Reads to the Reference
    #BISCUIT can be used to map, duplicate mark, sort, and index the provided reads.
    #dupsifter is a command line tool for marking PCR duplicates in both WGS and WGBS datasets.
    #os.system('/BIO1/biscuit/biscuit align -@ 96 '+REF+' '+sample+'_1.fastq.gz '+sample+'_2.fastq.gz | /BIO1/dupsifter/dupsifter_v1.0.0_linux_amd64 /BIO1/REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa | samtools sort -@ 96 -o '+sample+'.bam -O BAM - && samtools index '+sample+'.bam')

    os.system('samtools flagstat '+sample+'.bam > '+sample+'_flagstat.txt')

def STEP2(sample):
    if not os.path.isdir('TEMP_PICARD'):
        os.mkdir('TEMP_PICARD')
    os.system('java -Xmx16g -jar /BIO1/picard.jar MarkDuplicates TMP_DIR=TEMP_PICARD VALIDATION_STRINGENCY=LENIENT I='+sample+'.bam O='+sample+'.dedup.bam M='+sample+'.duplicate_metrics REMOVE_DUPLICATES=true AS=true CREATE_INDEX=true')

def STEP3(REF, sample):
    #Generating Pileup for a Single Sample
    os.system('/BIO1/biscuit/biscuit pileup -@ 96 -o '+sample+'_pileup.vcf '+REF+' '+sample+'.dedup.bam && bgzip -@ 96 '+sample+'_pileup.vcf && tabix -p vcf '+sample+'_pileup.vcf.gz')

    # Extract DNA methylation into BED format
    os.system('/BIO1/biscuit/biscuit vcf2bed -c '+sample+'_pileup.vcf.gz > '+sample+'_beta_m_u.bed')


def main():
    REF = '/BIO1/REF/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
    file_list = os.listdir('./')
    file_ls = [file for file in file_list if file.endswith('_1.fastq.gz')]
    file_ls.sort()
    for i in file_ls:
        sample = i.split('_1.fastq.gz')[0]
        print(sample)
        #STEP1_RefIndex(REF,sample)
        STEP2(sample)
        #STEP3(REF, sample)
main()
