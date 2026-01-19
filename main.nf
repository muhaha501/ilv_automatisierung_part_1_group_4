nextflow.enable.dsl=2

/* * PARAMETER
 * Wir nutzen jetzt das Transkriptom statt Genom+GTF
 */
params.reads = "$baseDir/data/*_{1,2}.fastq.gz"
params.transcriptome = "$baseDir/genome_Kallisto/Sus_scrofa.Sscrofa11.1.cdna.all.fa"
params.outdir = "$baseDir/results"

/*
 * 1. INDEXIEREN (Kallisto)
 * Braucht wenig RAM.
 */
process KALLISTO_INDEX {
    tag "Indexing"
    publishDir "${params.outdir}/kallisto_index", mode: 'copy'

    input:
    path transcriptome

    output:
    path "transcriptome.idx"

    script:
    """
    kallisto index -i transcriptome.idx ${transcriptome}
    """
}

/*
 * 2. QUALITÄTSKONTROLLE (FastQC)
 * 
 */
process FASTQC {
    tag "FASTQC on $sample_id"
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    path "*.{html,zip}"

    script:
    """
    fastqc -q ${reads}
    """
}

/*
 * 3. QUANTIFIZIERUNG (Kallisto)
 * Erledigt Mapping und Zählen in einem Schritt.
 */
process KALLISTO_QUANT {
    tag "Quantifying $sample_id"
    publishDir "${params.outdir}/kallisto_counts", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)
    path index

    output:
    path "${sample_id}"

    script:
    // -b 100: Erstellt 100 Bootstraps (wichtig für statistische Genauigkeit)
    // ${reads} wird automatisch zu "read1 read2" expandiert
    """
    kallisto quant -i ${index} -o ${sample_id} -b 100 ${reads}
    """
}

workflow {
    // Channels erstellen
    read_pairs_ch = Channel.fromFilePairs(params.reads, checkIfExists: true)
    
    // 1. Index bauen
    index_ch = KALLISTO_INDEX(file(params.transcriptome))
    
    // 2. QC (kann parallel laufen)
    FASTQC(read_pairs_ch)
    
    // 3. Quantifizieren (wartet auf Index)
    KALLISTO_QUANT(read_pairs_ch, index_ch)
}
