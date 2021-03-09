#! /usr/bin/env nextflow

params.help = null

log.info ""
log.info "--------------------------------------------------------"
log.info "  DE 1.0 : Pipeline for the DE "
log.info "--------------------------------------------------------"

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE : sudo nextflow run Lipinski-B/DE-nf -profile docker --input /data/ --STAR_Index /STARIndex/ --output /output/"
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run DE.nf [-r vX.X -profile singularity] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info ""
    log.info "--input                      FOLDER                      Folder where you can find your data (fasta/fastq files)."
    log.info "--output                      FOLDER                      Folder where you want to find your result."
    log.info ""
    log.info "Optional arguments:"
    log.info "--<OPTION>                      <TYPE>                      <DESCRIPTION>"
    log.info "--STAR_Index                      FOLDER                      Folder where you can find the STAR index."
    log.info "--annotation                      FOLDER                      Folder where you can find the annotation to use."
    log.info "--R                      STRING                      'on'/'off' : Chose to use or not the standard R analyses from the pipeline."
    log.info "--thread                      INT                      Number of thread to use."
    

    log.info ""
    log.info "Flags:"
    log.info "--<FLAG>                                                    <DESCRIPTION>"
    log.info ""
    exit 0
} else {
    /* Software information */
    log.info "help:                               ${params.help}"
}


// -- Path :
params.input = null
params.output = null
params.GTF = "/home/boris/Bureau/projet/projetS2/data/GCF_006496715.1_Aalbo_primary.1_genomic.gtf"

// -- Option :
params.R = "off"
params.thread = 1
params.STAR_Index = "off"
params.FNA = "/home/boris/Bureau/projet/projetS2/data/GCF_006496715.1_Aalbo_primary.1_genomic.fna"


// -- Pipeline :
process Mapping{ 
  publishDir params.output+'/mapping/', mode: 'copy'
  cpus params.thread
  
  input:
  file data from Channel.fromPath(params.input+'*').collect()
  //file data from Channel.fromPath(params.STAR_Index).collect()

  output:
  file "*Aligned.out.sam" into Mapping_sam
  file "*Log.out" into Mapping_Log
  
  shell:
  if(params.STAR_Index=="off") {
    '''
    mkdir STARIndex_last/
    STAR --runThreadN !{params.thread} \
      --runMode genomeGenerate \
      --genomeDir STARIndex_last/ \
      --genomeFastaFiles !{params.FNA} \
      --sjdbGTFfile !{params.GTF} \
      --sjdbOverhang 74 \
      --genomeSAsparseD 3

    mkdir data/
    mv *gz data/
    cd data/ 

    #Mapping analyse :
    ulimit -v 27598325157
    for file in *; do
      STAR \
      --runThreadN !{params.thread} \
      --genomeDir ../STARIndex_last/ \
      --readFilesCommand gunzip -c \
      --readFilesIn $file \
      --outFileNamePrefix $file \
      --outSAMunmapped Within
    done

    mv * ../.
    cd ..
    rm -r data/
    '''
  } else {
    '''
    #Mapping analyse :
    ulimit -v 27598325157
    for file in *; do
      STAR \
      --runThreadN !{params.thread} \
      --genomeDir !{params.STAR_Index} \
      --readFilesCommand gunzip -c \
      --readFilesIn $file \
      --outFileNamePrefix $file \
      --outSAMunmapped Within
    done
    '''
    }
}


process Intersection{ 
  publishDir params.output+'/intersect/', mode: 'copy'
  
  input:
  file data from Mapping_sam
  
  output:
  file "*.txt" into Intersect
  
  shell:
  '''
  #Intersection analyse :
  for file in *; do
    htseq-count --stranded=yes --nprocesses=!{params.thread} --mode=union $file !{params.GTF} > ${file}_intersect.txt
  done
  '''
}

process Merge_result{ 
  publishDir params.output+'/merge/', mode: 'copy'
  
  input:
  file data from Intersect
  
  output:
  file "finale.txt" into Result
  
  shell:
  '''
  #Differancial expression analyses :
  files=(*)
  awk '{print $1}' ${files[0]} > AAAA.txt

  for file in *_intersect.txt; do
    awk '{print $2}' $file > ${file}.tmp
    rm $file
    mv ${file}.tmp $file 
    tail -n +2 "$file" > "$file.tmp" && mv "$file.tmp" "$file"
    echo "${file%%.*}" > ${file}.name
    cat ${file}.name ${file} > ${file}.tmp && mv ${file}.tmp ${file}
    rm ${file}.name
  done

  paste -d "\t" * > finale.txt
  rm AAAA.txt
  '''
}


if(params.R=="on"){
  process DE{ 
    publishDir params.output+'/R/', mode: 'copy'
    
    input:
    file data from Result.collect()
    //file data from Channel.fromPath(params.input+'data/Metadata.xls').collect()
    
    output:
    file "*.pdf" into Result_DE
    
    shell:
    '''
    Rscript !{baseDir}/bin/DE.r finale.txt !{baseDir}/data/Metadata.xls
    '''
    }
}
