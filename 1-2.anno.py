import string,sys,os

def state4():
    os.chdir('/BIO6/REAL_ANALYSIS/2.metilene')
    fp = open('metilene_Tumor_Normal.filter_qval.0.05.anno.txt','r')
    fpout = open('metilene_Tumor_Normal.filter_qval.0.05.anno2.txt','w')
    hd = fp.readline()
    head = hd.rstrip().split('\t')
    temp1,temp2 = '\t'.join(head[0:3]),'\t'.join(head[6:])
    fpout.write(temp1+'\t'+"annotation"+'\t'+""+'\t'+temp2+'\n')
    for line in fp:
        line_temp = line.rstrip().split('\t')
        temp = '\t'.join(line_temp[6:])
        if 'Intron' in line_temp[5]:
            b = line_temp[5].split('(')[1][:-1]
            if 'intron 1 of' in line_temp[5]:
                fpout.write('\t'.join([line_temp[0],str(line_temp[1]),str(line_temp[2])])+'\t'+"1st Intron"+'\t'+b+'\t'+temp+'\n')
            else:
                fpout.write('\t'.join([line_temp[0],str(line_temp[1]),str(line_temp[2])])+'\t'+"Other Intron"+'\t'+b+'\t'+temp+'\n')
        elif 'Exon' in line_temp[5]:
            d = line_temp[5].split('(')[1][:-1]
            if 'exon 1 of' in line_temp[5]:
                fpout.write('\t'.join([line_temp[0],str(line_temp[1]),str(line_temp[2])])+'\t'+"1st Exon"+'\t'+d+'\t'+temp+'\n')
            else:
                fpout.write('\t'.join([line_temp[0],str(line_temp[1]),str(line_temp[2])])+'\t'+"Other Exon"+'\t'+d+'\t'+temp+'\n')

        else:
            fpout.write('\t'.join([line_temp[0],str(line_temp[1]),str(line_temp[2])])+'\t'+line_temp[5]+'\t'+""+'\t'+temp+'\n')

def main():
    state4()
main()
