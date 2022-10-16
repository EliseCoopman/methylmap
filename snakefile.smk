import os #This module provides a portable way of using operating system dependent functionality
import pandas as pd
from pathlib import Path


outdir = "/home/ecoopman/outputresults"
sample_info = pd.read_table("/home/ecoopman/meth_bam_files/sampleinfo_limited.txt").set_index("SAMPLE", drop=False)
sample_info.BAM = "/home/ecoopman/meth_bam_files/" + sample_info.BAM 
sample_info["FILE PRESENT"] = sample_info.BAM.apply(os.path.isfile)
sample_info = sample_info.loc[sample_info["FILE PRESENT"] == True]


rule all:
    input:
        os.path.join(outdir, "heatmap/heatmap.tsv"),
        os.path.join(outdir, "heatmap/heatmap.html"),


rule heatmap:
    input:
        inp = os.path.join("/home/ecoopman/methfreqtable_chr_limited.tsv.gz"),
        #inp = expand(os.path.join(outdir, "calculatemethylationfrequency/testfiles/{id}_phase{phase}.tsv.gz"), id=sample_info.index, phase = [1,2]),
        #tbi = expand(os.path.join(outdir, "calculatemethylationfrequency/testfiles/{id}_phase{phase}.tsv.gz.tbi"), id=sample_info.index, phase = [1,2]),
        gff = os.path.join("/home/ecoopman/ONT-meth-Elise/gff3/gencode_v40_annotation_sorted.gff3.gz"),
    params:
        names = ["name1", "name2", "name3", "name4"],
        groups = ["case", "control", "control", "case"]
    output:
       outtable = os.path.join(outdir, "heatmap/heatmap.tsv"),
       outfig = os.path.join(outdir, "heatmap/heatmap.html"),
    conda:
        "/home/ecoopman/ONT-meth-Elise/envs/heatmap.yml"
    threads: 8
    log:
        os.path.join(outdir, "logs/heatmap.log")
    shell:
         "python heatmap_test.py --table {input.inp} --window chr17:44345246-44353106 --gff {input.gff} --expand 10000 --groups {params.groups} --outfig {output.outfig} --outtable {output.outtable} 2> {log}"
