#!usr/bin/env/ python3


file_out = open("longueur.txt",'w')

with open('GCF_006496715.1_Aalbo_primary.1_genomic.gff','r') as files:
    for lines in files.readlines():
        if str(lines)[0] != "#" :
            line=lines.split('\t')

            if line[2]=='gene' :
                start=line[3]
                stop=line[4]
                info=line[8].split(';')
                name=info[0].replace('ID=gene-',"")
                name=name.replace("'","")
                gene = name + '\t' + str(int(stop)-int(start)) + '\n'
                file_out.write(gene)