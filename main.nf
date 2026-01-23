nextflow.enable.dsl=2

// Parameters
params.star_index  = null 
params.docker      = true
params.tool        = "STAR"  

// STAR PARAMS
params.max_forks   = 5


// Files and Directories
params.reads       = "data/*_{1,2}.fastq.gz"
params.genome_url  = "https://ftp.ensembl.org/pub/release-111/fasta/sus_scrofa/dna/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa.gz"
params.gtf_url     = "https://ftp.ensembl.org/pub/release-111/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.111.gtf.gz"
params.transcriptome_url = "https://ftp.ensembl.org/pub/release-115/fasta/sus_scrofa/cdna/Sus_scrofa.Sscrofa11.1.cdna.all.fa.gz"
params.outdir      = "results"


// FASTQC Process
process FASTQC {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*.{html,zip}")

    script:
    """
    fastqc ${reads} 
    """
}

// Building of the Index (saved in the results/index directory)
process BUILD_STAR_INDEX {
    container 'quay.io/biocontainers/star:2.6.1d--0'
    publishDir "${params.outdir}/index", mode: 'copy'

    input:
    path fasta_gz
    path gtf_gz

    output:
    path "star_index"

    script:
    """
    zcat ${fasta_gz} > genome.fa
    zcat ${gtf_gz} > genome.gtf


    mkdir star_index
    STAR --runMode genomeGenerate \
         --genomeDir star_index \
         --genomeFastaFiles genome.fa \
         --sjdbGTFfile genome.gtf \
         --runThreadN 8
    
    rm genome.fa genome.gtf
    """
}

process STAR_ALIGN {
    container 'quay.io/biocontainers/star:2.6.1d--0'
    publishDir "${params.outdir}/alignment", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index_dir 

    output:
    tuple val(sample_id), path("${sample_id}Aligned.out.bam")

    script:
    """
    STAR --genomeDir $index_dir \
         --readFilesIn $reads \
         --readFilesCommand zcat \
         --outFileNamePrefix ${sample_id} \
         --outSAMtype BAM Unsorted \
         --runThreadN 8
    """
}

process KALLISTO_INDEX {
    container 'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1'
    publishDir "${params.outdir}/index/kallisto_index", mode: 'copy'

    input:
    path transcriptome

    output:
    path "transcriptome.idx"

    script:
    """
    zcat ${transcriptome} > transcriptome.fa
    kallisto index -i transcriptome.idx transcriptome.fa
    rm transcriptome.fa
    """
}

process KALLISTO_QUANT {
    container 'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1'
    publishDir "${params.outdir}/counts/kallisto", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    tuple val(sample_id), path("${sample_id}")

    script:
    """
    kallisto quant -i ${index} -o ${sample_id} -b 100 ${reads}
    """
}


process FEATURECOUNTS {
    container 'quay.io/biocontainers/subread:2.0.1--h7d7f7ad_1'
    publishDir "${params.outdir}/counts/STAR", mode: 'copy'
    
    input:
    tuple val(sample_id), path(bam)
    path gtf_gz

    output:
    tuple val(sample_id), path("${sample_id}/${sample_id}.featureCounts.txt")
    tuple val(sample_id), path("${sample_id}/${sample_id}.featureCounts.txt.summary")

    script:
    """
    zcat ${gtf_gz} > annotation.gtf

    featureCounts -p -s 2 -T 8 \\
                  -a annotation.gtf \\
                  -o ${sample_id}/${sample_id}.featureCounts.txt \\
                  ${bam}

    rm annotation.gtf
    """
}

workflow {
    read_pairs_ch = channel.fromFilePairs(params.reads)

    FASTQC(read_pairs_ch)

    if (params.tool == "STAR") 
        {
         run_star(read_pairs_ch)
        }
    else if (params.tool == "KALLISTO")
        {
         run_kallisto(read_pairs_ch)
        }
    else
        {
         log.error "Unsupported tool specified: ${params.tool}. Please choose either 'STAR' or 'KALLISTO'."
        }
}

workflow run_star {
    take: reads_ch

    main:
    gtf_ch = channel.fromPath(params.gtf_url).collect()
    // Check if STAR index is provided
    if (params.star_index) {
        // Existing Index used
        index_ch = channel.fromPath(params.star_index, checkIfExists: true).collect()
    } else {
        // Build Index from Remote Sources
        fasta_ch = channel.fromPath(params.genome_url).collect()
        index_ch = BUILD_STAR_INDEX(fasta_ch, gtf_ch).collect()
    }

    bam_ch = STAR_ALIGN(reads_ch, index_ch)
    FEATURECOUNTS(bam_ch, gtf_ch)
}

workflow run_kallisto {
    take: reads_ch
    
    main:
    index_ch = KALLISTO_INDEX(file(params.transcriptome_url)).collect()
    KALLISTO_QUANT(reads_ch, index_ch)
}