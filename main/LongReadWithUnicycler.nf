#!/user/bin/env nextflow

nextflow.enable.dsl=2

//Initial task of pipeline
//Run trimming for adapters, put through filtlong, and run fastqc

//input params
params.out
params.sample
//true or false
params.isolate


println """\
         H Y B R I D  A S S E M B L Y  W I T H  U N I C Y C L E R
         ========================================================
         Output Directory               : ${params.out    }
         Sample Name or ID              : ${params.sample }
         Isolate non chromosomal contigs: ${params.isolate}
         """
         .stripIndent()
//Channels
lr_ch = Channel.fromPath("*_long.fastq.gz", checkIfExists: true)
sr_ch = Channel.fromFilePairs("*{_1,_2,_R1,_R2}_001.{fastq,fq}.gz", checkIfExists: true)

workflow{
  //precprocess and concatenate fastqs
  concatenate(lr_ch)

  //trim short reads
  trimSR(sr_ch)

  //trims adapters
  trim(concatenate.out.concatTrimmed)

  //obtains read metrics
  nanoPlot(trim.out.poreOut)

  //filters concatenated fastq
  filtLong(trim.out.poreOut)

  //assembly
  uniCyclerAssembly(trimSR.out.srTrimmed.combine(filtLong.out.filtered))

  //IsolateContigs
  pullCircularizedContigs(uniCyclerAssembly.out.assembly)

  //AssemblyQC
  AssemblyQC(uniCyclerAssembly.out.assembly)
}

process trimSR{
  tag{"Trimmomatic on Short Reads"}

  publishDir("${params.out}/SrTrimmed", mode: 'copy')

  input:
  tuple val(sample),path(shReads)

  output:
  path("*${sample}{_1,_2}.fastq.gz"), emit: srTrimmed
  path("*.txt")

  script:
  """
  java -jar /Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 4 ${shReads} -baseout ${sample}_short.fastq.gz SLIDINGWINDOW:10:20 MINLEN:100 > ${sample}.trim.stats.txt
  mv ${sample}_short_1P.fastq.gz ${sample}_1.fastq.gz
  mv ${sample}_short_2P.fastq.gz ${sample}_2.fastq.gz
  """

}

process concatenate{
  tag{"Concatenate Fastqs"}
  label 'process_low'

  input:
  path(barcode)

  output:
  path("${params.sample}_long.fastq.gz"), emit: concatTrimmed

  script:
  """
  mkdir ${params.out}
  cat ${barcode} > ${params.sample}_long.fastq.gz
  """
}

process trim{
  //run poreChop on multiple files from input
  tag{"PoreChop"}
  label 'process_low'

  publishDir("${params.out}/trimmed", mode: 'copy')

  input:
  path(reads)

  output:
  path("${params.sample}_long_trimmed.fastq.gz"), emit: poreOut

  script:
  """
  ~/Porechop/porechop-runner.py  -t 8 -i ${reads} -o ${params.sample}_long_trimmed.fastq.gz
  """
}

process nanoPlot{
  tag{"NanoPlot"}
  label 'process_low'

  publishDir("${params.out}", mode: 'copy')

  input:
  path(reads)

  output:
  path("*")

  script:
  """
  NanoPlot --fastq ${reads} --N50 --tsv_stats -o nanoPlotStats
  """
}

process filtLong{
  tag{"FiltLong"}
  label 'process_low'

  publishDir("${params.out}/filtered", mode: 'copy')

  input:
  path(reads)

  output:
  path("${params.sample}_filtered.fastq.gz"), emit: filtered

  script:
  """
  filtlong --keep_percent 95 ${reads} | gzip > ${params.sample}_filtered.fastq.gz
  """
}

process uniCyclerAssembly{
  tag{"Unicycler Hybrid Asssembler"}
  label 'process_low'

  publishDir("${params.out}", mode: 'copy')

  input:
  path(reads)

  output:
  path("*")
  path("${params.sample}_Assembly/assembly.fasta"), emit: assembly

  script:
  """
  unicycler -1 ${reads[0]} -2 ${reads[1]} -l ${reads[2]} -o ${params.sample}_Assembly
  """
}

process AssemblyQC{
  tag{"QUAST"}
  label 'process_low'

  publishDir("${params.out}", mode: 'copy')

  input:
  path(assembly)

  output:
  path("quastOut/*")

  script:
  """
  quast.py -o quastOut ${assembly}
  """
}

process pullCircularizedContigs{
  tag{"Isolate Circularized Contigs"}

  publishDir("${params.out}", mode: 'copy')

  input:
  path(assembly)

  output:
  path("*.fasta")

  when:
  params.isolate = "true"

  script:
  """
  #!/usr/bin/python

  import sys
  import os.path
  #Assembly Filename
  assemFile = "${assembly}"
  #Name of output file
  outfileName = "IsolatedContigs.fasta"
  headerFound = False
  secondHeaderFound = False
  try:
      with open(assemFile, 'r') as infile:
          lines = infile.readlines()
          #If the user wants all contigs other than the first, presumably
          #chromosomal contig, this flag will grab all of them and 
          #add them to a new file.
          for line in lines:
              if secondHeaderFound:
                  outfile.write(line)
              #Looks for the second contig as that is the second
              #smallest, excluding the largest, the chromosomal contig
              if ">2" in line:
                  secondHeaderFound = True
                  outfile = open(outfileName, 'a')
                  outfile.write(line)

  except:
      print("Assembly file does not exist.")
  """
}
