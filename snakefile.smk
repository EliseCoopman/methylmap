import os #This module provides a portable way of using operating system dependent functionality
import pandas as pd
from pathlib import Path


outdir = "/home/ecoopman/outputresults"
sample_info = pd.read_table("/home/ecoopman/meth_bam_files/sampleinfo_limited.txt").set_index("SAMPLE", drop=False)
sample_info.BAM = "/home/ecoopman/meth_bam_files/" + sample_info.BAM 
sample_info["FILE PRESENT"] = sample_info.BAM.apply(os.path.isfile)
sample_info = sample_info.loc[sample_info["FILE PRESENT"] == True]

sample_info_baminput = pd.read_table("/home/ecoopman/heatmap/baminput/bamfiles/bamfiles.txt").set_index("SAMPLE")
#sample_info_baminput = sample_info_baminput.rename(columns=sample_info_baminput.iloc[0])
sample_info_baminput.BAM = "/home/ecoopman/heatmap/baminput/bamfiles/" + sample_info_baminput.BAM
sample_info_baminput["FILE PRESENT"] = sample_info_baminput.BAM.apply(os.path.isfile)
sample_info_baminput = sample_info_baminput.loc[sample_info_baminput["FILE PRESENT"] == True]

rule all:
    input:
        os.path.join(outdir, "methylmap/methylmap.tsv"),
        os.path.join(outdir, "methylmap/methylmap.html"),


rule methylmap:
    input:
        inp = os.path.join("/home/ecoopman/heatmap/methfreqtable/methfreqtable_chr_limited.tsv.gz"),
        #inp = expand(os.path.join(outdir, "calculatemethylationfrequency/testfiles/{id}_phase{phase}.tsv.gz"), id=sample_info.index, phase = [1,2]),
        #tbi = expand(os.path.join(outdir, "calculatemethylationfrequency/testfiles/{id}_phase{phase}.tsv.gz.tbi"), id=sample_info.index, phase = [1,2]),
        #bam = expand(os.path.join("/home/ecoopman/heatmap/baminput/bamfiles/{id}.bam"), id=sample_info_baminput.index),
        #fasta = os.path.join("/home/ecoopman/heatmap/baminput/hg38.fa.gz"),
        gff = os.path.join("/home/ecoopman/ONT-meth-Elise/gff3/gencode_v40_annotation_sorted.gff3.gz"),
    #params:
        #names = ["name1", "name2", "name3", "name4","name5"],
        #groups = ["case", "control", "control", "case"]
    output:
       outtable = os.path.join(outdir, "methylmap/methylmap.tsv"),
       outfig = os.path.join(outdir, "methylmap/methylmap.html"),
    conda:
        "/home/ecoopman/ONT-meth-Elise/envs/heatmap.yml"
    threads: 8
    log:
        os.path.join(outdir, "logs/methylmap.log")
    shell:
        "python methylmap/methylmap.py --table {input.inp} --gff {input.gff} --outfig {output.outfig} --outtable {output.outtable} 2> {log}"
