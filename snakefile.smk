import os #This module provides a portable way of using operating system dependent functionality
import pandas as pd
from pathlib import Path

outdir = "/home/ecoopman/outputresults"
sample_info = pd.read_table("/home/ecoopman/meth_bam_files/sampleinfo.txt").set_index("SAMPLE", drop=False)
sample_info.BAM = "/home/ecoopman/meth_bam_files/" + sample_info.BAM 
sample_info["FILE PRESENT"] = sample_info.BAM.apply(os.path.isfile)
sample_info = sample_info.loc[sample_info["FILE PRESENT"] == True]


rule all:
    input:
        os.path.join(outdir, "heatmap/heatmap.tsv"),
        os.path.join(outdir, "heatmap/heatmap.html"),

rule heatmap:
    input:
        files = expand(os.path.join(outdir, "calculatemethylationfrequency/{id}_phase{phase}.tsv.gz"), id=sample_info.index, phase = [1,2]),
        tbi = expand(os.path.join(outdir, "calculatemethylationfrequency/{id}_phase{phase}.tsv.gz.tbi"), id=sample_info.index, phase = [1,2]),
        gff = os.path.join("/home/ecoopman/ONT-meth-Elise/gff3/gencode_v40_annotation_sorted.gff3.gz"),
    output:
       outtable = os.path.join(outdir, "heatmap/heatmap.tsv"),
       outfig = os.path.join(outdir, "heatmap/heatmap.html"),
    conda:
        "/home/ecoopman/ONT-meth-Elise/envs/heatmap.yml"
    threads: 8
    log:
        os.path.join(outdir, "logs/heatmap.log")
    shell:
         "python heatmap_test.py --files {input.files} --window chr17:44345246-44353106 --expand 10000 --gff {input.gff} --outtable {output.outtable} --outfig {output.outfig} 2> {log}"


