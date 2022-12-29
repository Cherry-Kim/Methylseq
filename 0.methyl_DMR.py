import string,glob,sys,os

def STEP1_METHYL(sample,REF):
    os.system('trim_galore --paired  --gzip -o trim_galore/ '+sample+'_R1.fastq.gz '+sample+'_R2.fastq.gz  -q 20 --length 20')
    os.system('/BIO1/Bismark-master/bismark --bowtie2 --bam '+REF+' -1 trim_galore/'+sample+'_R1_val_1.fq.gz -2 trim_galore/'+sample+'_R2_val_2.fq.gz -p 48')
    os.system('/BIO1/Bismark-master/deduplicate_bismark --bam '+sample+'_R1_val_1_bismark_bt2_pe.bam')
    os.system('/BIO1/Bismark-master/bismark_methylation_extractor '+sample+'_R1_val_1_bismark_bt2_pe.deduplicated.bam --bedGraph')

def STEP2_DMR(PATH,Group1,Group2,Group3,in2_T,in1_N,h2_T,h1_N):
    #2-1. Convert .bedgraph
    #1.Remove first line 2. Print not the "KI" lines 3. Add "chr" 
    G = [Group1,Group2,Group3]
    for k in G:
        print(k)
        file_list = os.listdir(PATH+k)
        print(len(file_list))
        bed_list = [file for file in file_list if file.endswith('_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph')]
        co = 0
        for i in bed_list:
            co +=1
            sample = i.split('_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph')[0]
            print (co,sample)
            os.system("""sed '1d' """+PATH+k+"""/"""+sample+"""_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph | grep -v -P "KI" | grep -v -P "GL" | awk '{print "chr"$1 "\t" $2 "\t" $3 "\t" $4}' > """+sample+"""_bismark.dedup.chr.bedGraph""")
            os.system('/BIO1/bedtools2/bin/sortBed -i '+sample+'_bismark.dedup.chr.bedGraph > '+sample+'_bismark.dedup.chr.sort.bedGraph')
            os.system('mv *_bismark.dedup.chr.bedGraph '+PATH+k)
            os.system('mv *chr.sort.bedGraph '+PATH+k)
    
    #2-2. Generate metilene input files
    os.system('cp '+PATH+Group1+'/*_bismark.dedup.chr.sort.bedGraph /BIO6/REAL_ANALYSIS/2.metilene/')
    os.system('cp '+PATH+Group2+'/*_bismark.dedup.chr.sort.bedGraph /BIO6/REAL_ANALYSIS/2.metilene/')
    os.chdir('/BIO6/REAL_ANALYSIS/2.metilene/')
    os.system('perl /BIO1/metilene_v0.2-8/metilene_input.pl --in1 '+in2_T+' --in2 '+in1_N+' --h1 '+Group1+' --h2 '+Group2+' --out metilene_'+Group1+'_'+Group2+'.input')

    #fpout = open('head.txt','w')
    #fpout.write('\t'.join(['chrom','pos',h2_T,h1_N]) +'\n')
    #os.system('cat head.txt metilene_'+Group1+'_'+Group2+'.input > metilene_'+Group1+'_'+Group2+'.head.input')

    #2-3. Identufy DMRs
    os.system('/BIO1/metilene_v0.2-8/metilene -t 48  -a '+Group1+' -b '+Group2+' metilene_'+Group1+'_'+Group2+'.input | sort -V -k1,1 -k2,2n >  metilene_'+Group1+'_'+Group2+'.output')
    
    #2-4. Filter DMRs
    os.system('/BIO1/metilene_v0.2-8/metilene_output.pl -q  metilene_'+Group1+'_'+Group2+'.output -o metilene_'+Group1+'_'+Group2+'.filter -a '+Group1+' -b '+Group2)
    fpout = open('head2.txt','w')
    fpout.write('\t'.join(['Chr','start','stop','q-value','Mean_meth.diff.','#CpG','Mean_meth.Tumor','Mean_meth.Normal']) +'\n')
    os.system('cat head2.txt metilene_'+Group1+'_'+Group2+'.filter_qval.0.05.out > metilene_'+Group1+'_'+Group2+'.filter_qval.0.05.head.out')

    #2-5. rGREAT input file
    os.system("""sed '1 i\chr\tstart\tend' metilene_Tumor_Normal.filter_qval.0.05.bedgraph | awk '{print $1 "\t" $2 "\t" $3}' > metilene_Tumor_Normal.filter_qval.0.05.rGREAT.input.bed""")
    os.system('/usr/bin/Rscript ./rGREAT.r')

def main():
    REF = '/BIO5/'  #Homo_sapiens.GRCh38.dna.primary_assembly.fa
    co = 0
    fp = glob.glob('*_R1.fastq.gz')
    for fname in fp:
        co += 1
        sample = fname.split('_R1.fastq.gz')[0]
        print( "### SAMPLE ",co, sample)
        #STEP1_METHYL(sample,REF)
main()

def main2():
    PATH,Group1,Group2,Group3 = '/BIO6/REAL_ANALYSIS/bedGraph/','Tumor','Normal','WBC'
    co = 0
    T,N,H =[],[],[]
    fp1 = open('mapfile.txt','r')
    fp1.readline()
    for line in fp1:
        line_temp = line[:-1].split('\t')   #['1', '1. 49624-104 (Normal)', '49624-104_Normal_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz', '32', 'Female', 'White', 'Normal Colorectal Mucosa']
        if line_temp[6] == 'Tumor':
            s = line_temp[2].split('_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz')[0]
            T.append(s)
        elif line_temp[6] == 'Normal Colorectal Mucosa':
            s1 = line_temp[2].split('_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz')[0]
            N.append(s1)
        elif line_temp[6] == 'Healthy Control WBC':
            s2 = line_temp[2].split('_R1_val_1_bismark_bt2_pe.deduplicated.bedGraph.gz')[0]
            H.append(s2)
#    print (len(T))
    h2_T =  "\t".join(map(lambda z: z+"_T", T)) #header #49624-104_Tumor_T
    in2_T =  ",".join(map(lambda z: z+"_bismark.dedup.chr.sort.bedGraph", T))   #input #49624-104_Tumor_bismark.dedup.chr.sort.bedGraph
#    print (len(N))
    h1_N =  "\t".join(map(lambda z: z+"_N", N))
    in1_N =  ",".join(map(lambda z: z+"_bismark.dedup.chr.sort.bedGraph", N))
#    print (len(H))
    h3_H =  "\t".join(map(lambda z: z+"_H", H))
    in3_H =  ",".join(map(lambda z: z+"_bismark.dedup.chr.sort.bedGraph", H))
    STEP2_DMR(PATH,Group1,Group2,Group3,in2_T,in1_N,h2_T,h1_N)
main2()
