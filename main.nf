#!/usr/bin/env nextflow
params.vcfs = "s3://bdtest2003/tgen/sample1000a.vcf.gz"
params.bedTSV = "s3://bdtest2003/tgen/GRCh38/v3.1-GRCh38-all-stratifications.tsv"
params.outputFolder = "s3://bdtest2003/tgen/results"

nextflow.enable.dsl=2

process tabix {
    input:
        tuple val(meta), path(vcf)

    output:
        tuple val(meta), path("*vcf.gz"), emit: vcf
        path "*vcf.gz*", emit: publishFiles

   // container = 'tabix'
    cpu = 2
    mem = 2
    publishDir path: "${params.outputFolder}", mode: 'move', overwrite: 'true'
    
    script:
    """
        bgzip ${vcf}
        tabix ${vcf}.gz
    """
}

workflow METADATA {
    
    take:
        f
        
    main:
        fileTuple = f.map{ [ [ id: it.simpleName, baseName: it.baseName, size: it.size(), folder: it.parent, fileType: it.extension, fullFile: it ], it ] }

    emit:
        fileTuple

}

process sam_index {
    input:
        tuple val(sample), path(bam)

    output:
        path "${bam}*" 

    //container = 'bcftools'
    cpu = 6
    mem = 2

    publishDir path: "${params.outputFolder}", mode: 'move', overwrite: 'true'

    script:
        """
            samtools index "${bam}" -@ ${task.cpus} "${bam}.bai" 
        """

}

process intersectBam {

    input:
        tuple val(sample), path(bam), val(bed), path(bedFile)

    output:
        path "${sample.id}----${bed}.bam"

    tag "${bam}----${bed}"

   // container = 'bedtools'
    cpu = 2
    mem = 2

    script:
    """
        intersectBed -a ${bam} -b ${bedFile} > "${sample.id}----${bed}.bam"
    """
}

process intersectVcf {

    input:
        tuple val(sample), path(vcf), val(bed), path(bedFile)

    output:
        path "${sample.id}----${bed}.vcf"

 //   container = 'bedtools'
    cpu = 8
    mem = 32
    tag "${vcf}----${bed}"
    
    script:
    """
        intersectBed -a ${vcf} -b ${bedFile} -header > "${sample.id}----${bed}.vcf"
    """
}

workflow {
    beds  = Channel.fromPath(params.bedTSV).splitCsv(header: false, sep: "\t")
    beds.map { [ bed: it[0], bedFile: it[1] ] }
//    bams = Channel.fromPath(params.bams) | METADATA
    vcfs = Channel.fromPath(params.vcfs) | METADATA
//    bams = bams.combine(beds).map { [ sample: it[0], bam: it[1], bed: it[2], bedFile: 's3://bdtest2003/tgen/GRCh38/' + it[3] ] }
    vcfs = vcfs.combine(beds).map { [ sample: it[0], vcf: it[1], bed: it[2], bedFile: 's3://bdtest2003/tgen/GRCh38/' + it[3] ] }

//    sam_index(intersectBam(bams) | METADATA)
    tabix(intersectVcf(vcfs) | METADATA)
}
