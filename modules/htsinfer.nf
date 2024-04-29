#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process HTSINFER_SE {
    label 'htsinfer'
    tag "${sample}"
    publishDir "${params.htsinfer_out_dir}/${sample}", mode: 'copy', pattern: "*"

    input:
    tuple val(sample), path(fq1)

    output:
    tuple val(sample), path('*')

    script:
        """
        htsinfer \
            --output-directory ${params.htsinfer_out_dir}/${sample} \
            --temporary-directory $TMP \
            --cleanup-regime KEEP_ALL \
            --records ${params.records} \
            --threads ${params.threads} \
            --read-layout-min-match-percentage ${params.read_min_percent} \
            --read-layout-min-frequency-ratio ${params.read_min_freq} \
            --library-source-min-match-percentage  ${params.lib_min_percent} \
            --library-source-min-frequency-ratio  ${params.lib_min_freq} \
            --library-type-max-distance  ${params.lib_max_dist} \
            --library-type-mates-cutoff  ${params.lib_mates_cutoff} \
            --read-orientation-min-mapped-reads  ${params.read_orient_min_mapped} \
            --read-orientation-min-fraction  ${params.read_orient_min_freq} \
            --verbosity DEBUG \
            ${fq1} \
            > ${sample}_result.json \
            2> ${sample}_error.txt
        """
}

process HTSINFER_PE {
    label 'htsinfer'
    tag "${sample}"
    publishDir "${params.htsinfer_out_dir}/${sample}", mode: 'copy', pattern: "*"

    input:
    tuple val(sample), path(fq1), path(fq2)

    output:
    tuple val(sample), path('*')

    script:
        """
        htsinfer \
            --output-directory ${params.htsinfer_out_dir}/${sample} \
            --temporary-directory $TMP \
            --cleanup-regime KEEP_ALL \
            --records ${params.records} \
            --threads ${params.threads} \
            --read-layout-min-match-percentage ${params.read_min_percent} \
            --read-layout-min-frequency-ratio ${params.read_min_freq} \
            --library-source-min-match-percentage  ${params.lib_min_percent} \
            --library-source-min-frequency-ratio  ${params.lib_min_freq} \
            --library-type-max-distance  ${params.lib_max_dist} \
            --library-type-mates-cutoff  ${params.lib_mates_cutoff} \
            --read-orientation-min-mapped-reads  ${params.read_orient_min_mapped} \
            --read-orientation-min-fraction  ${params.read_orient_min_freq} \
            --verbosity DEBUG \
            ${fq1} ${fq2} \
            > ${sample}_result.json \
            2> ${sample}_error.txt
        """
}
