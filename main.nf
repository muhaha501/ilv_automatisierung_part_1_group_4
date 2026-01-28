nextflow.enable.dsl=2


// Process to download SRA data using prefetch + fasterq-dump (only if not already present)
process DOWNLOAD_SRA {
    container 'quay.io/biocontainers/sra-tools:3.0.3--h87f3376_0'
    publishDir "${params.fastq_dir}", mode: 'copy'
    
    maxForks 3  // Limit concurrent downloads to avoid overwhelming NCBI servers
    
    input:
    val(srr_id)
    
    output:
    // Emit R1 and R2 explicitly; we convert to a pair later in the workflow
    tuple val(srr_id), path("${srr_id}_1.fastq.gz"), path("${srr_id}_2.fastq.gz")
    
    when:
    // Only run if files don't already exist
    !file("${params.fastq_dir}/${srr_id}_1.fastq.gz").exists() || !file("${params.fastq_dir}/${srr_id}_2.fastq.gz").exists()
    
    script:
    """
    echo "Downloading ${srr_id} from SRA..."
    
    # Download SRA file
    prefetch ${srr_id} --max-size 50G
    
    # Convert to FASTQ (paired-end)
    fasterq-dump ${srr_id} --split-files --threads 4
    
    # Compress output files
    gzip ${srr_id}_1.fastq
    gzip ${srr_id}_2.fastq
    
    # Clean up SRA cache
    rm -rf ${srr_id}
    
    echo "Download complete for ${srr_id}"
    """
}

// FASTQC Process
process FASTQC {
    container 'quay.io/biocontainers/fastqc:0.11.9--0'
    publishDir "${params.outdir}/fastqc", mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    // Emit both HTML and ZIP explicitly
    tuple val(sample_id), path("*.html"), path("*.zip")

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
    each index_dir 

    output:
    tuple val(sample_id), path("${sample_id}Aligned.out.bam")

    script:
    """
    STAR --genomeDir ${index_dir} \
         --readFilesIn ${reads} \
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
    each index

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
    each gtf_gz

    output:
    tuple val(sample_id), path("${sample_id}.featureCounts.txt"), emit: counts
    tuple val(sample_id), path("${sample_id}.featureCounts.txt.summary"), emit: summary

    script:
    """
    zcat ${gtf_gz} > annotation.gtf

    featureCounts -p -s 2 -T 8 \\
                  -a annotation.gtf \\
                  -o ${sample_id}.featureCounts.txt \\
                  ${bam}

    rm annotation.gtf
    """
}

// ============================================
// POST-PROCESSING: QC, Visualization & DESeq2
// ============================================

// Bulk RNA-seq QC and Exploratory Analysis
process RNASEQ_QC_ANALYSIS {
    container 'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61571f297571a-0'
    publishDir "${params.outdir}/analysis/qc", mode: 'copy'

    input:
    path count_files

    output:
    path "count_matrix.csv", emit: counts
    path "sample_metadata.csv", emit: metadata
    path "figures/*"
    path "qc_report.txt"
    path "qc_data.pkl"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from sklearn.decomposition import PCA
    import pickle
    import os

    plt.switch_backend('Agg')
    os.makedirs("figures", exist_ok=True)

    # Sample metadata
    metadata = {
        'SRR12835320': 'CTR', 'SRR12835321': 'CTR', 'SRR12835322': 'CTR',
        'SRR12835323': 'CTR', 'SRR12835324': 'CTR', 'SRR12835325': 'CTR',
        'SRR12835326': 'PFA', 'SRR12835327': 'PFA', 'SRR12835328': 'PFA',
        'SRR12835329': 'PFA', 'SRR12835330': 'PFA', 'SRR12835331': 'PFA',
        'SRR12835332': 'OO', 'SRR12835333': 'OO', 'SRR12835334': 'OO',
        'SRR12835335': 'OO', 'SRR12835336': 'OO', 'SRR12835337': 'OO',
        'SRR12835338': 'PFA_OO', 'SRR12835339': 'PFA_OO', 'SRR12835341': 'PFA_OO',
        'SRR12835342': 'PFA_OO', 'SRR12835343': 'PFA_OO', 'SRR12835344': 'PFA_OO',
    }
    colors = {'CTR': '#1f77b4', 'PFA': '#ff7f0e', 'OO': '#2ca02c', 'PFA_OO': '#d62728'}

    # Load count files
    count_files = [f for f in os.listdir('.') if f.endswith('.featureCounts.txt')]
    
    # Build count matrix
    counts = {}
    for f in count_files:
        sample = f.replace('.featureCounts.txt', '')
        df = pd.read_csv(f, sep='\\t', comment='#')
        counts[sample] = df.set_index('Geneid').iloc[:, -1]
    
    count_matrix = pd.DataFrame(counts)
    samples = [s for s in count_matrix.columns if s in metadata]
    count_matrix = count_matrix[samples]
    
    # Create metadata dataframe
    meta_df = pd.DataFrame({'group': [metadata[s] for s in samples]}, index=samples)
    
    # Save CSV
    count_matrix.to_csv("count_matrix.csv")
    meta_df.to_csv("sample_metadata.csv")

    # Normalize
    log_cpm = np.log2(count_matrix / count_matrix.sum() * 1e6 + 1)

    # QC Report
    lib_sizes = count_matrix.sum()
    genes_detected = (count_matrix > 0).sum()
    
    with open("qc_report.txt", "w") as f:
        f.write(f"Samples: {len(samples)}\\n")
        f.write(f"Genes: {len(count_matrix)}\\n")
        f.write(f"Library sizes: {lib_sizes.min():.0f} - {lib_sizes.max():.0f}\\n")
        f.write(f"Genes detected: {genes_detected.min()} - {genes_detected.max()}\\n")

    # Plot 1: Library sizes
    fig, ax = plt.subplots(figsize=(12, 4))
    c = [colors[metadata[s]] for s in samples]
    ax.bar(samples, lib_sizes / 1e6, color=c)
    ax.set_ylabel("Counts (millions)")
    ax.set_title("Library Sizes")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig("figures/01_library_sizes.png", dpi=150)
    plt.close()

    # PCA analysis
    pca = PCA(n_components=5)
    pca_data = pca.fit_transform(log_cpm.T)
    pca_df = pd.DataFrame(pca_data, index=samples, columns=[f'PC{i+1}' for i in range(5)])
    pca_df['group'] = [metadata[s] for s in samples]

    # Plot 2: PCA PC1 vs PC2
    fig, ax = plt.subplots(figsize=(8, 6))
    for group in colors:
        mask = pca_df['group'] == group
        ax.scatter(pca_df.loc[mask, 'PC1'], pca_df.loc[mask, 'PC2'], c=colors[group], label=group, s=80)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.legend()
    ax.set_title("PCA: PC1 vs PC2")
    plt.tight_layout()
    plt.savefig("figures/02_pca_pc1_pc2.png", dpi=150)
    plt.close()

    # Plot 3: PCA PC1 vs PC3
    fig, ax = plt.subplots(figsize=(8, 6))
    for group in colors:
        mask = pca_df['group'] == group
        ax.scatter(pca_df.loc[mask, 'PC1'], pca_df.loc[mask, 'PC3'], c=colors[group], label=group, s=80)
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC3 ({pca.explained_variance_ratio_[2]*100:.1f}%)")
    ax.legend()
    ax.set_title("PCA: PC1 vs PC3")
    plt.tight_layout()
    plt.savefig("figures/03_pca_pc1_pc3.png", dpi=150)
    plt.close()

    # Plot 4: PCA PC2 vs PC3
    fig, ax = plt.subplots(figsize=(8, 6))
    for group in colors:
        mask = pca_df['group'] == group
        ax.scatter(pca_df.loc[mask, 'PC2'], pca_df.loc[mask, 'PC3'], c=colors[group], label=group, s=80)
    ax.set_xlabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_ylabel(f"PC3 ({pca.explained_variance_ratio_[2]*100:.1f}%)")
    ax.legend()
    ax.set_title("PCA: PC2 vs PC3")
    plt.tight_layout()
    plt.savefig("figures/04_pca_pc2_pc3.png", dpi=150)
    plt.close()

    # Plot 5: Variance explained
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(range(1, 6), pca.explained_variance_ratio_ * 100, color='steelblue')
    ax.set_xlabel("PC")
    ax.set_ylabel("Variance Explained (%)")
    ax.set_title("PCA Variance Explained")
    plt.tight_layout()
    plt.savefig("figures/05_pca_variance.png", dpi=150)
    plt.close()

    # Plot 6: Sample correlation
    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(log_cpm.corr(), cmap='RdYlBu_r', vmin=0.8, vmax=1, ax=ax)
    ax.set_title("Sample Correlation")
    plt.tight_layout()
    plt.savefig("figures/06_correlation.png", dpi=150)
    plt.close()

    # Save objects as pickle
    qc_data = {
        'count_matrix': count_matrix,
        'log_cpm': log_cpm,
        'metadata': meta_df,
        'pca_model': pca,
        'pca_data': pca_df,
        'lib_sizes': lib_sizes,
        'genes_detected': genes_detected
    }
    with open("qc_data.pkl", "wb") as f:
        pickle.dump(qc_data, f)

    print("QC complete!")
    """
}

// Differential Expression Analysis
process DIFFERENTIAL_EXPRESSION {
    container 'quay.io/biocontainers/mulled-v2-8186960447c5cb2faa697666dc1e6d919ad23f3e:3127fcae6b6bdaf8181e21a26ae61571f297571a-0'
    publishDir "${params.outdir}/analysis/deseq2", mode: 'copy'

    input:
    path count_matrix
    path metadata

    output:
    path "de_results/"
    path "figures/"
    path "de_data.pkl"

    script:
    """
    #!/usr/bin/env python3
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sns
    from scipy import stats
    import pickle
    import os

    plt.switch_backend('Agg')
    os.makedirs("de_results", exist_ok=True)
    os.makedirs("figures", exist_ok=True)

    # Load data
    counts = pd.read_csv("count_matrix.csv", index_col=0)
    meta = pd.read_csv("sample_metadata.csv", index_col=0)
    
    colors = {'CTR': '#1f77b4', 'PFA': '#ff7f0e', 'OO': '#2ca02c', 'PFA_OO': '#d62728'}

    # Filter genes (at least 10 counts in 3+ samples)
    keep = (counts >= 10).sum(axis=1) >= 3
    counts = counts.loc[keep]

    # Normalize (simple CPM)
    cpm = counts / counts.sum() * 1e6
    log_cpm = np.log2(cpm + 1)

    # DE function
    def compare_groups(g1, g2):
        s1 = meta[meta['group'] == g1].index
        s2 = meta[meta['group'] == g2].index
        
        results = []
        for gene in counts.index:
            v1 = log_cpm.loc[gene, s1].values
            v2 = log_cpm.loc[gene, s2].values
            
            log2fc = v2.mean() - v1.mean()
            _, pval = stats.mannwhitneyu(v1, v2, alternative='two-sided') if (v1.std() > 0 or v2.std() > 0) else (0, 1)
            
            results.append({'gene': gene, 'log2FC': log2fc, 'pvalue': pval})
        
        df = pd.DataFrame(results)
        df['padj'] = np.minimum(1, df['pvalue'] * len(df) / (df['pvalue'].rank()))
        return df.sort_values('pvalue')

    # Run comparisons
    comparisons = [('CTR', 'PFA'), ('CTR', 'OO'), ('CTR', 'PFA_OO'), ('OO', 'PFA_OO')]
    all_results = {}

    for g1, g2 in comparisons:
        name = f"{g2}_vs_{g1}"
        df = compare_groups(g1, g2)
        df.to_csv(f"de_results/{name}.csv", index=False)
        all_results[name] = df
        
        sig = df[(df['padj'] < 0.05) & (df['log2FC'].abs() > 1)]
        print(f"{name}: {len(sig)} DE genes")

    # Volcano plots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    for ax, (name, df) in zip(axes.flat, all_results.items()):
        sig = (df['padj'] < 0.05) & (df['log2FC'].abs() > 1)
        ax.scatter(df['log2FC'], -np.log10(df['pvalue']), c='gray', alpha=0.3, s=5)
        ax.scatter(df.loc[sig, 'log2FC'], -np.log10(df.loc[sig, 'pvalue']), c='red', s=5)
        ax.axvline(1, ls='--', c='gray', lw=0.5)
        ax.axvline(-1, ls='--', c='gray', lw=0.5)
        ax.set_xlabel('log2 Fold Change')
        ax.set_ylabel('-log10(p-value)')
        ax.set_title(name)
    plt.tight_layout()
    plt.savefig("figures/01_volcano.png", dpi=150)
    plt.close()

    # Top genes heatmap
    top_genes = all_results['PFA_OO_vs_CTR'].head(30)['gene'].tolist()
    if len(top_genes) >= 5:
        data = log_cpm.loc[top_genes]
        data_z = (data.T - data.mean(axis=1)) / data.std(axis=1)
        
        fig, ax = plt.subplots(figsize=(12, 8))
        sns.heatmap(data_z.T, cmap='RdBu_r', center=0, ax=ax)
        ax.set_title("Top 30 DE Genes (PFA_OO vs CTR)")
        plt.tight_layout()
        plt.savefig("figures/02_heatmap.png", dpi=150)
        plt.close()

    # Summary bar plot
    fig, ax = plt.subplots(figsize=(8, 5))
    names = list(all_results.keys())
    up = [len(all_results[n][(all_results[n]['padj'] < 0.05) & (all_results[n]['log2FC'] > 1)]) for n in names]
    down = [len(all_results[n][(all_results[n]['padj'] < 0.05) & (all_results[n]['log2FC'] < -1)]) for n in names]
    
    x = np.arange(len(names))
    ax.bar(x - 0.2, up, 0.4, label='Up', color='red', alpha=0.7)
    ax.bar(x + 0.2, down, 0.4, label='Down', color='blue', alpha=0.7)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=30, ha='right')
    ax.set_ylabel('DE Genes')
    ax.legend()
    plt.tight_layout()
    plt.savefig("figures/03_summary.png", dpi=150)
    plt.close()

    # Save all results as pickle
    de_data = {
        'counts_filtered': counts,
        'cpm': cpm,
        'log_cpm': log_cpm,
        'metadata': meta,
        'de_results': all_results,
        'comparisons': comparisons
    }
    with open("de_data.pkl", "wb") as f:
        pickle.dump(de_data, f)

    print("DE analysis complete!")
    """
}

process KALLISTO_TO_COUNTS {
    container 'quay.io/biocontainers/kallisto:0.46.2--h4f7b962_1'

    input:
    tuple val(sample_id), path(kallisto_dir)

    output:
    path "${sample_id}.featureCounts.txt"

    script:
    """
    echo "# Program:kallisto; Command:quant" > ${sample_id}.featureCounts.txt
    
    echo -e "Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t${sample_id}" >> ${sample_id}.featureCounts.txt
    
    tail -n +2 ${kallisto_dir}/abundance.tsv | awk -F'\\t' '{OFS="\\t"; print \$1, "NA", "NA", "NA", "NA", \$2, \$4}' >> ${sample_id}.featureCounts.txt
    """
}

workflow {
    // Determine input source: SRA download or local files
    if (params.download_sra) {
        // Generate SRA IDs from SRR12835320 to SRR12835344 (excluding SRR12835340)
        sra_ids_ch = Channel
            .from(12835320..12835344)
            .filter { it != 12835340 }  // SRR12835340 doesn't exist in the dataset
            .map { "SRR${it}" }
        
        // Download only files that don't exist yet
        downloaded_ch = DOWNLOAD_SRA(sra_ids_ch)
            .map { sid, r1, r2 -> tuple(sid, [r1, r2]) }  // make it look like fromFilePairs output
        
        // Combine downloaded files with existing local files
        existing_ch = Channel.fromFilePairs("${params.fastq_dir}/*_{1,2}.fastq.gz")
        read_pairs_ch = existing_ch.mix(downloaded_ch)
        
        log.info "Downloading SRA data: SRR12835320 to SRR12835344"
        log.info "Existing files in ${params.fastq_dir} will be reused"
    } else {
        // Use local files only
        read_pairs_ch = Channel.fromFilePairs(params.reads)
        log.info "Using local FASTQ files: ${params.reads}"
    }

    FASTQC(read_pairs_ch)

    if (params.tool == "STAR") {
        run_star(read_pairs_ch)
    }
    else if (params.tool == "KALLISTO") {
        run_kallisto(read_pairs_ch)
    }
    else {
        log.error "Unsupported tool specified: ${params.tool}. Please choose either 'STAR' or 'KALLISTO'."
    }
}

workflow run_star {
    take: reads_ch

    main:
    // Prepare GTF channel (in DSL2, channels can be reused directly)
    gtf_ch = Channel.fromPath(params.gtf_url)

    // Check if STAR index is provided
    if (params.star_index) {
        // Use existing index
        index_ch = Channel.fromPath(params.star_index, checkIfExists: true)
    } else {
        // Build index from remote sources
        fasta_ch = Channel.fromPath(params.genome_url)
        index_ch = BUILD_STAR_INDEX(fasta_ch, gtf_ch)
    }

    bam_out = STAR_ALIGN(reads_ch, index_ch)
    fc_out  = FEATURECOUNTS(bam_out, gtf_ch)

    // Collect all count files (first output channel: counts)
    all_counts = fc_out.counts
        .map { sample_id, fc_file -> fc_file }
        .collect()

    run_post_processing(all_counts)
}

workflow run_kallisto {
    take: reads_ch
    
    main:
    transcriptome_ch = Channel.fromPath(params.transcriptome_url)
    index_ch         = KALLISTO_INDEX(transcriptome_ch)
    quant_out = KALLISTO_QUANT(reads_ch, index_ch)

    kallisto_counts = KALLISTO_TO_COUNTS(quant_out)

    all_counts = kallisto_counts.collect()

    run_post_processing(all_counts)
}

workflow run_post_processing {
    take: all_counts

    main:
    // Run QC and exploratory analysis
    qc_results = RNASEQ_QC_ANALYSIS(all_counts)
    
    // Run differential expression analysis
    DIFFERENTIAL_EXPRESSION(
        qc_results.counts,
        qc_results.metadata
    )
}