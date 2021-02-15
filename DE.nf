#! /usr/bin/env nextflow

params.help = null

log.info ""
log.info "--------------------------------------------------------"
log.info "  DE 1.0 : Pipeline for the DE "
log.info "--------------------------------------------------------"

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE : sudo nextflow run Lipinski-B/DE-nf -profile docker "
    log.info "  USAGE : sudo nextflow run ../DE-nf/DE.nf -profile docker --input /media/bobo/Seagate/all/ --STAR_Index ~/media/bobo/Seagate/STARIndex"
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run DE.nf [-r vX.X -profile singularity] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info ""
    log.info "--input                      FOLDER                      Folder where you can find your data (fasta/fastq files)"
    log.info ""
    log.info "Optional arguments:"
    log.info "--<OPTION>                      <TYPE>                      <DESCRIPTION>"
    log.info "--STAR_Index                      FOLDER                      Folder where you can find the STAR index"
    log.info "--annotation                      FOLDER                      Folder where you can find your annotation"
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

// -- Option :
params.STAR_Index = null
params.thread = 1
params.annotation = 1

// -- Path :
params.input = null
params.result = null

// -- Pipeline :
process Mapping{ 
  //publishDir params.result+'!{params.result}/mapping/', mode: 'move'
  
  input:
  file data from Channel.fromPath(params.input+'*').collect()

  output:
  file "*.Aligned.out.sam" into Mapping
  file "*.Log.out" into Mapping_Log
  
  shell:
  '''
  #Mapping analyse :
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


process Intersection{ 
  publishDir params.result+'!{params.result}/intersection/', mode: 'move'
  
  input:
  //file data from Mapping
  
  output:
  //file "STAR/" into Intersect
  
  shell:
  '''
  #Intersection analyse :
  #htseq-count -m union --stranded=no -r pos $4 $Annotation > $5
  #htseq-count --version
  '''
}

process DE{ 
  publishDir params.result+'!{params.result}/', mode: 'copy'
  
  input:
  //file data from Intersect
  
  output:
  //file "STAR/" into Result
  
  shell:
  '''
  #Differancial expression analyses :
  '''
}
