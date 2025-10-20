nextflow.enable.dsl=2

params.samplesheet = params.samplesheet ?: 'samples.csv'     // genotype,condition,sample,R1,R2
params.genome_fa   = params.genome_fa   ?: 'GRCm39.primary_assembly.genome.fa'
params.gtf         = params.gtf         ?: 'gencode.vM32.primary_assembly.annotation.gtf'
params.adapters    = params.adapters    ?: 'TruSeq3-PE.fa'
params.outdir      = params.outdir      ?: 'results'
params.star_index  = params.star_index  ?: "${params.outdir}/star_index"
params.threads     = params.threads     ?: 8

// Optional: set containers (Biocontainers) or conda envs
process.container = ''
// process.conda = 'envs/bulk_rnaseq.yml'     // if you prefer conda

// Publish dirs
process.publishDir = { "${params.outdir}/${task.process}".toString() }

// ------------------------------------------------------------
// Helpers
// ------------------------------------------------------------
Channel.fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row ->
        def id = [row.genotype, row.condition, row.sample].join('_')
        tuple(id, file(row.R1), file(row.R2))
    }
    .set { READS }

Channel.value(file(params.genome_fa)).set { GENOME_FA }
Channel.value(file(params.gtf)).set { GTF }

// Build STAR index once
process STAR_INDEX {
    tag "STAR index"
    cpus Math.max(params.threads as int, 4)
    publishDir "${params.outdir}/STAR_INDEX", mode: 'copy'

    input:
    path GENOME_FA
    path GTF

    output:
    path "genomeDir", emit: idx

    script:
    """
    mkdir -p genomeDir
    STAR --runThreadN ${task.cpus} \
         --runMode genomeGenerate \
         --genomeDir genomeDir \
         --genomeFastaFiles ${GENOME_FA} \
         --sjdbGTFfile ${GTF} \
         --sjdbOverhang 100
    """
}

// FASTQC raw
process FASTQC_RAW {
    tag "$id"
    cpus 4
    publishDir "${params.outdir}/fastqc_raw", mode: 'copy'

    input:
    tuple val(id), path(R1), path(R2)

    output:
    tuple val(id), path(R1), path(R2), path("*.zip"), path("*.html")

    script:
    """
    fastqc -t ${task.cpus} ${R1} ${R2}
    """
}

// Trimmomatic PE
process TRIMMOMATIC {
    tag "$id"
    cpus 8
    publishDir "${params.outdir}/trimmed", mode: 'copy'

    input:
    tuple val(id), path(R1), path(R2), path(qz), path(qh)
    path params.adapters

    output:
    tuple val(id),
          path("${id}_R1.trim.fq.gz"),
          path("${id}_R1.unpaired.fq.gz"),
          path("${id}_R2.trim.fq.gz"),
          path("${id}_R2.unpaired.fq.gz")

    script:
    """
    trimmomatic PE -threads ${task.cpus} -phred33 \
      ${R1} ${R2} \
      ${id}_R1.trim.fq.gz ${id}_R1.unpaired.fq.gz \
      ${id}_R2.trim.fq.gz ${id}_R2.unpaired.fq.gz \
      ILLUMINACLIP:${params.adapters}:2:30:10 \
      LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}

// STAR align (sorted BAM)
process STAR_ALIGN {
    tag "$id"
    cpus Math.max(params.threads as int, 8)
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(id),
          path(R1trim),
          path(R1un),
          path(R2trim),
          path(R2un)
    path star_idx from STAR_INDEX.out.idx

    output:
    tuple val(id), path("${id}.Aligned.sortedByCoord.out.bam")

    script:
    """
    STAR --runThreadN ${task.cpus} \
         --genomeDir ${star_idx} \
         --readFilesIn ${R1trim} ${R2trim} \
         --readFilesCommand zcat \
         --outSAMtype BAM SortedByCoordinate \
         --outFileNamePrefix ${id}.
    mv ${id}.Aligned.sortedByCoord.out.bam ${id}.Aligned.sortedByCoord.out.bam
    """
}

// samtools index
process SAMTOOLS_INDEX {
    tag "$id"
    publishDir "${params.outdir}/bam", mode: 'copy'

    input:
    tuple val(id), path(bam)

    output:
    tuple val(id), path(bam), path("${bam}.bai")

    script:
    """
    samtools index ${bam}
    """
}

// samtools flagstat
process SAMTOOLS_FLAGSTAT {
    tag "$id"
    publishDir "${params.outdir}/flagstat", mode: 'copy'

    input:
    tuple val(id), path(bam), path(bai)

    output:
    path("${id}.flagstat.txt")

    script:
    """
    samtools flagstat ${bam} > ${id}.flagstat.txt
    """
}

// featureCounts (all BAMs together)
process FEATURECOUNTS {
    tag "featureCounts"
    cpus Math.max(params.threads as int, 8)
    publishDir "${params.outdir}/featurecounts", mode: 'copy'

    input:
    path GTF
    tuple val(id), path(bam), path(bai) from SAMTOOLS_INDEX.out.collect()

    output:
    path "featureCounts.txt"
    path "featureCounts.txt.summary"

    script:
    def bamList = bam.collect{ it }.join(' ')
    """
    featureCounts -p -T ${task.cpus} -a ${GTF} -o featureCounts.txt ${bamList}
    """
}

// MultiQC at the end
process MULTIQC {
    tag "multiqc"
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path("*") from Channel.of(
        file("${params.outdir}/fastqc_raw"),
        file("${params.outdir}/trimmed"),
        file("${params.outdir}/bam"),
        file("${params.outdir}/flagstat"),
        file("${params.outdir}/featurecounts")
    )

    output:
    path "multiqc_report.html"
    path "multiqc_data"

    script:
    """
    multiqc ${params.outdir} --force
    """
}

// ------------------------------------------------------------
// Workflow
// ------------------------------------------------------------
workflow {
    take: READS, GENOME_FA, GTF

    main:
    raw = FASTQC_RAW(READS)
    trimmed = TRIMMOMATIC(raw.out, file(params.adapters))
    aligned = STAR_ALIGN(trimmed.out, STAR_INDEX(GENOME_FA, GTF))
    indexed = SAMTOOLS_INDEX(aligned.out)
    SAMTOOLS_FLAGSTAT(indexed.out)
    FEATURECOUNTS(GTF, indexed.out)
    MULTIQC()
}
