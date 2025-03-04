import os
import pandas as pd

#Paths
study = "/results/rr/study/hg38s/study267-1000Genomes_ONT/1000Genomes_dorado/"
studytemp = "/results/rr/study/hg38s/study271-temp_EC/1000Genomes_dorado/"
ref = "/home/ecoopman/study279_methylmap_files/hg38.no_alt.fa" #samtools faidx hg38.no_alt.fa

bamfiles = pd.read_table("/home/ecoopman/study279_methylmap_files/first_100_plus_bam_download_limited.txt").set_index("SAMPLE_ID", drop=False)

def get_BAM(wildcards):
    return bamfiles.loc[wildcards.sample_id, "URL"]

# Define the maximum number of concurrent downloads
max_downloads = 10
chrom=["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY","chrM"]
parts = ["1_25000000", "25000001_50000000", "50000001_75000000", "75000001_100000000", "100000001_125000000", "125000001_150000000", "150000001_175000000", "175000001_200000000", "200000001_225000000", "225000001_250000000","250000001_275000000","275000001_300000000"]

rule all:
    input:
        expand(os.path.join(study, "modkit_files/{sample_id}/H_1.bed"), sample_id=bamfiles.index),
        expand(os.path.join(study, "modkit_files_reorg/{sample_id}/{sample_id}_H{hp}.bed.gz"), sample_id=bamfiles.index, hp=[1,2]),
        expand(os.path.join(study, "modkit_files_reorg_splitchrom/{sample_id}/{sample_id}_H{hp}_chrM.tsv"), sample_id=bamfiles.index, hp=[1,2]),
        expand(os.path.join(study, "modkit_files_reorg_splitchrom_splitchunks/{sample_id}/{sample_id}_H{hp}_{chrom}_275000001_300000000.tsv"), sample_id=bamfiles.index, hp=[1,2], chrom=chrom),
        expand(os.path.join(study, "modkit_files_reorg_splitchrom_splitchunks_reorg/{sample_id}/{sample_id}_H{hp}_{chrom}_{parts}.tsv.gz"), sample_id=bamfiles.index, hp=[1,2], chrom=chrom, parts=parts),
        expand(os.path.join(study, "modkit_files_joinchunck/{chrom}/joined_{chrom}_{parts}.tsv.gz"), chrom=chrom, parts=parts),
        expand(os.path.join(study, "modkit_files_joinchunck_joinchrom/{chrom}/joined_{chrom}.tsv.gz"), chrom=chrom),
        os.path.join(study, "1000Genomes.tsv.gz"),
        os.path.join(study, "1000Genomes_sorted.tsv.gz"),
        os.path.join(study, "1000Genomes_header.tsv.gz"),
        os.path.join(study, "modified_header.tsv.gz"),
        os.path.join(study, "1000Genomes_all.tsv.gz"),
        os.path.join(study, "1000Genomes_all.tsv.gz.tbi")

rule download_and_modkit:
    output:
        bam_file = temp(os.path.join(studytemp, "bam_files/{sample_id}.bam")),
        bai_file = temp(os.path.join(studytemp, "bam_files/{sample_id}.bam.bai")),
        bed_file = os.path.join(study, "modkit_files/{sample_id}/H_1.bed"),
        log_file = os.path.join(study, "logfiles/{sample_id}.log")
    params:
        ref = ref,
        get_BAM = get_BAM
    resources:
        download_slots = max_downloads
    shell:
        """
        # Download BAM file
        wget -O {output.bam_file} {params.get_BAM}
        wget -O {output.bai_file} {params.get_BAM}.bai
        
        # Run modkit
        outdir=$(echo {output.bed_file} | sed 's|/[^/]*$||')
        modkit pileup {output.bam_file} $outdir \
        --ref {params.ref} \
        --cpg \
        --ignore "h" \
        --only-tabs \
        --partition-tag HP \
        --prefix H \
        --log-filepath {output.log_file}
        """

rule reorganize_output:
    input:
        bedfile = os.path.join(study, "modkit_files/{sample_id}/H_{hp}.bed"),
    output:
        bedfile_reorg = os.path.join(study, "modkit_files_reorg/{sample_id}/{sample_id}_H{hp}.bed.gz")
    log:
        os.path.join(study, "logfiles/modkit_file_reorg/{sample_id}_{hp}.log")
    shell:
        """
        python 1000Genomes_reorgoutput.py --input {input.bedfile}  2> {log} | bgzip > {output.bedfile_reorg}
        """

rule split_into_chromosomes:
    input:
        os.path.join(study, "modkit_files_reorg/{sample_id}/{sample_id}_H{hp}.bed.gz")
    output:
        os.path.join(study, "modkit_files_reorg_splitchrom/{sample_id}/{sample_id}_H{hp}_chrM.tsv")
    log:
        os.path.join(study, "logfiles/split_chrom_{sample_id}_{hp}.log")
    shell:
        """
        outdir=$(dirname {output})
        python 1000Genomes_split_per_chrom.py --input {input} --outputdir $outdir 2> {log}
        """

rule split_chrom_in_chunks:
    input:
        os.path.join(study, "modkit_files_reorg_splitchrom/{sample_id}/{sample_id}_H{hp}_{chrom}.tsv")
    output: 
        os.path.join(study, "modkit_files_reorg_splitchrom_splitchunks/{sample_id}/{sample_id}_H{hp}_{chrom}_275000001_300000000.tsv")
    log:
        os.path.join(study, "logfiles/split_chunk_{sample_id}_{hp}_{chrom}.log")
    shell:
        """
        outdir=$(dirname {output})
        python 1000Genomes_split_chrom_in_chunks.py --input {input} --outputdir $outdir 2> {log}
        """

rule reorg_files:
    input:
        os.path.join(study, "modkit_files_reorg_splitchrom_splitchunks/{sample_id}/{sample_id}_H{hp}_{chrom}_{parts}.tsv")
    output:
        os.path.join(study, "modkit_files_reorg_splitchrom_splitchunks_reorg/{sample_id}/{sample_id}_H{hp}_{chrom}_{parts}.tsv.gz")
    log:
        os.path.join(study, "logfiles/reorg_files_{sample_id}_{hp}_{chrom}_{parts}.log")
    shell:
        """
        python 1000Genomes_positioncolumn.py --input {input} 2> {log} | bgzip > {output}
        """

rule join_chunks:
    input:
        expand(os.path.join(study, "modkit_files_reorg_splitchrom_splitchunks_reorg/{sample_id}/{sample_id}_H{hp}_{{chrom}}_{{parts}}.tsv.gz"), sample_id=bamfiles.index, hp=[1,2]),
    output:
        os.path.join(study, "modkit_files_joinchunck/{chrom}/joined_{chrom}_{parts}.tsv.gz") #niet meer sorted -> changed code
    log:
        os.path.join(study, "logfiles/joinchunks/join_chunks_{chrom}_{parts}.log")
    shell:
        """
        python 1000Genomes_joinchunks.py --input {input} 2> {log} | bgzip > {output}
        """


rule join_perchrom:
    input:
        expand(os.path.join(study, "modkit_files_joinchunck/{{chrom}}/joined_{{chrom}}_{parts}.tsv.gz"), parts=parts),
    output:
        os.path.join(study, "modkit_files_joinchunck_joinchrom/{chrom}/joined_{chrom}.tsv.gz")
    log:
        os.path.join(study, "logfiles/join_perchrom_{chrom}.log")
    shell:
        """
        python 1000Genomes_joinchrom.py --input {input} 2> {log} | bgzip > {output}
        """


rule join_all:
    input:
        expand(os.path.join(study, "modkit_files_joinchunck_joinchrom/{chrom}/joined_{chrom}.tsv.gz"), chrom=chrom),
    output:
        os.path.join(study, "1000Genomes.tsv.gz")
    log:
        os.path.join(study, "logfiles/join_all.log")
    shell: 
        """
        zcat {input} | awk '!/^chrom/ || FNR == 1' | bgzip  > {output} 2> {log}
        """


rule sort_1000Genomes:
    input:
        os.path.join(study, "1000Genomes.tsv.gz")
    output:
        os.path.join(study, "1000Genomes_sorted.tsv.gz")
    params:
        tem_dir = "/home/ecoopman/temp"
    log:
        os.path.join(study, "logfiles/process_1000Genomes.log")
    shell: 
        """
        (zcat {input} | tail -n +2 | sort -k1,1V -k2,2n --temporary-directory={params.tem_dir} | bgzip > {output})  2> {log} 
        """
        


rule header_1000Genomes:
    input:
        os.path.join(study, "1000Genomes.tsv.gz")
    output:
        os.path.join(study, "1000Genomes_header.tsv.gz")
    log:
        os.path.join(study, "logfiles/header_1000Genomes.log")
    shell: 
        """
        zcat {input} | sed -n '1p' | bgzip > {output} 2> {log}
        """


rule paste_header:
    input:
        header = os.path.join(study, "1000Genomes_header.tsv.gz"),
        data = os.path.join(study, "1000Genomes_sorted.tsv.gz")
    output:
        modified_header = os.path.join(study, "modified_header.tsv.gz"),
        table_all = os.path.join(study, "1000Genomes_all.tsv.gz")
    log:
        os.path.join(study, "logfiles/paste_header.log")
    shell: 
        """
        zcat {input.header} | awk -v OFS='\t' 'NR==1{{$2="position"; print $0; next}} {{print $0}}' | bgzip > {output.modified_header} 2>> {log}

        zcat {output.modified_header} {input.data} | bgzip > {output.table_all} 2>> {log}
        """

rule index_1000Genomes:
    input:
        os.path.join(study, "1000Genomes_all.tsv.gz")
    output:
        os.path.join(study, "1000Genomes_all.tsv.gz.tbi")
    log:
        os.path.join(study, "logfiles/index_1000Genomes.log")
    shell: 
        """
        tabix -S 1 -s 1 -b 2 -e 2 {input} 2> {log}
        """   

