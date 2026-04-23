#!/usr/bin/env nextflow

// Copyright (C) 2026 IARC/WHO
// This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
// This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
// See the GNU General Public License for more details <http://www.gnu.org/licenses/>.

nextflow.enable.dsl=2

// ---------------------------
// PARAMETERS
// ---------------------------

params.help = null
params.bams = null
params.bam_csv = null
params.fasta= null
params.fai = null

//Header for the IARC tools - logo generated using the following page : http://patorjk.com/software/taag  (ANSI logo generator)
def IARC_Header (){
     return  """
#################################################################################
# ██╗ █████╗ ██████╗  ██████╗██████╗ ██╗ ██████╗ ██╗███╗   ██╗███████╗ ██████╗  #
# ██║██╔══██╗██╔══██╗██╔════╝██╔══██╗██║██╔═══██╗██║████╗  ██║██╔════╝██╔═══██╗ #
# ██║███████║██████╔╝██║     ██████╔╝██║██║   ██║██║██╔██╗ ██║█████╗  ██║   ██║ #
# ██║██╔══██║██╔══██╗██║     ██╔══██╗██║██║   ██║██║██║╚██╗██║██╔══╝  ██║   ██║ #
# ██║██║  ██║██║  ██║╚██████╗██████╔╝██║╚██████╔╝██║██║ ╚████║██║     ╚██████╔╝ #
# ╚═╝╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝╚═════╝ ╚═╝ ╚═════╝ ╚═╝╚═╝  ╚═══╝╚═╝      ╚═════╝  #
# Nextflow pipelines for cancer genomics.########################################
"""
}

//header for bam2cram tool (uses ANSI colors to make a short tool description)
//useful url: http://www.lihaoyi.com/post/BuildyourownCommandLinewithANSIescapecodes.html
def tool_header (){
        return """
        BAM\u001b[32;1m 2\u001b[33;1m CRAM\u001b[31;1m (${workflow.manifest.version})
        """
}

// ---------------------------
// PARAMETER CHECKS
// ---------------------------

if (!params.fasta || !params.fai)
    exit 1, "The reference fasta file and/or its index are missing"

if( !params.bams && !params.bam_csv ) {
    exit 1, "No bams are provided!"
}

// ---------------------------
// PROCESSES
// ---------------------------

process BAM2CRAM {

    tag "${sample}-2cram"
    publishDir "${params.output_folder}/CRAM", mode: 'copy'

    input:
        tuple val(sample), file(bam), file(bam_index)
        file(ref)
        file(ref_index)

    output:
        tuple val(sample), file("${sample}.cram"), file("${sample}.cram.crai"), emit: cramfiles

    script:
    """
    samtools view -C -T ${ref} ${bam} -o ${sample}.cram
    samtools index ${sample}.cram
    """
}

process STATS_BAMS{

    tag "${sample}-bam_stats"
    publishDir "${params.output_folder}/qc/bam", mode: 'copy'

    input:
        tuple val(sample), file(bam), file(bam_index)

    output:
        tuple val(sample), file("${sample}.bam.flagstat"), file("${sample}.bam.stats")

    script:
    """
    samtools flagstat ${bam} > ${sample}.bam.flagstat
    samtools stats ${bam} > ${sample}.bam.stats
    """
}

process STATS_CRAMS {

    tag "${sample}-cram_stats"
    publishDir "${params.output_folder}/qc/cram", mode: 'copy'

    input:
        tuple val(sample), file(cram), file(cram_index)

    output:
        tuple val(sample), file("${sample}.cram.flagstat"), file("${sample}.cram.stats")

    script:
    """
    samtools flagstat ${cram} > ${sample}.cram.flagstat
    samtools stats ${cram} > ${sample}.cram.stats
    """
}

// Process to check if bam_qc and cram_qc are identical

process CHECK_CONVERSION {

    tag "${sample}-check"
    publishDir "${params.output_folder}/qc/check", mode: 'copy'

    input:
        tuple val(sample), file(cram), file(cram_index), file(bam), file(bam_index)
        tuple val(sample), file(c_fs), file(c_stat), file(b_fs), file(b_stat)

    output:
        file("${sample}_check.report.txt")

    script:
    """
    if diff ${c_fs} ${b_fs} > /dev/null
    then
        fs_test="OK"
    else
        fs_test="fail"
    fi
    #Remove the line of stat that is diferent because of file extension
    grep -v "# The command line was:" ${c_stat} > ${c_stat}.no_cmd_line
    grep -v "# The command line was:" ${b_stat} > ${b_stat}.no_cmd_line

    if diff ${c_stat}.no_cmd_line ${b_stat}.no_cmd_line > /dev/null
    then
        s_test="OK"
    else
        s_test="fail"
    fi

    c_size=\$(du -Hh ${cram} | cut -f1)
    b_size=\$(du -Hh ${bam} | cut -f1)

    echo "ID\tflagstat\tstats\tBAM_size\tCRAM_size" > ${sample}_check.report.txt
    echo "${sample}\t\$fs_test\t\$s_test\t\$b_size\t\$c_size" >> ${sample}_check.report.txt
    """
}

// ---------------------------
// WORKFLOW
// ---------------------------

workflow {

  	log.info IARC_Header()
    log.info tool_header()
// --------------------------------------------------
// INFO / HELP
// --------------------------------------------------

log.info ""
log.info "----------------------------------------------------------------------------------------------------------------"
log.info "  bam2cram-nf : Nextflow pipeline to convert bam to cram  "
log.info "----------------------------------------------------------------------------------------------------------------"
log.info "Copyright (C) IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it under certain conditions; see LICENSE for details."
log.info "----------------------------------------------------------------------------------------------------------------"
log.info ""

if (params.help) {
  log.info ''
	log.info '-------------------------------------------------------------'
  log.info 'USAGE: '
  log.info 'nextflow run iarcbioinfo/bam2cram-nf --bams input_folder --fasta ref.fa --fai ref.fa.fai -profile singularity'
	log.info '-------------------------------------------------------------'
	log.info ''
  log.info 'Mandatory arguments:'
  log.info '    --bams   FOLDER                  Input folder containing BAM files and indexes.'
  log.info '    --fasta          FILE                  Path to fasta reference to encode the CRAM file'
  log.info '    --fai          FILE                  Path to fasta reference index.'
  log.info 'Input alternatives:'
  log.info '    --bam_csv          FILE                  Path to file with tabular data for each sample to process [label bam index ]'
  log.info 'Output:'
  log.info '    --output_folder   FOLDER                  Ouput folder [default: ./results]'
  log.info ''
  exit 0
}
else {
      /* Software information */
   log.info "bams = ${params.bams}"
   log.info "bam_csv = ${params.bam_csv}"
   log.info "fasta         = ${params.fasta}"
   log.info "output_folder= ${params.output_folder}"
   log.info "help=${params.help}"
 }

// SET INPUTS

    ch_fasta = Channel.value(file(params.fasta))
        .ifEmpty{ exit 1, "Fasta file not found: ${params.fasta}" }

    ch_fai = Channel.value(file(params.fai))
        .ifEmpty{ exit 1, "Fasta index file not found: ${params.fai}" }

if(params.bam_csv) {
    inputbams = Channel
        .fromPath(params.bam_csv)
        .ifEmpty { exit 1, "CSV file not found: ${params.bam_csv}" }
        .splitCsv(header: true, sep: '\t', strip: true)
        .map { row -> tuple(row.label, file(row.bam), file(row.index)) }
        .ifEmpty { exit 1, "CSV was empty - no input files supplied" }

    } else {
        if(file(params.bams).listFiles().findAll { it.name ==~ /.*bam/ }.size() > 0){

            bams = Channel
                .fromPath("${params.bams}/*.bam")
                .map { path -> tuple(path.name.replace(".bam",""), path) }

            bams_index = Channel
                .fromPath("${params.bams}/*.bam.bai")
                .map { path -> tuple(path.name.replace(".bam.bai",""), path) }

            inputbams = bams.join(bams_index)

        } else {
            println "ERROR: input folder ${params.bams} contains no BAM files"
            System.exit(1)
        }
    }

    // Duplicate channels explicitly
    bam4stats  = inputbams
    bam4check  = inputbams


    // ---------------------------
    // RUN PROCESSES
    // ---------------------------

    cram_out = BAM2CRAM(inputbams, ch_fasta, ch_fai)
	cram4check = cram_out

    bam_qc   = STATS_BAMS(bam4stats)
    cram_qc  = STATS_CRAMS(cram_out)

    cram_bam_list = cram4check.join(bam4check)
    cram_bam_qc   = cram_qc.join(bam_qc)

    report_qc = CHECK_CONVERSION(cram_bam_list, cram_bam_qc)

    report_qc
        .collectFile(
            name: 'bam2cram_summary.txt',
            storeDir: params.output_folder,
            seed: 'ID\tflagstat\tstats\tBAM_size\tCRAM_size\n',
            newLine: false,
            skip: 1
        )

}

