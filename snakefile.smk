import os
import pandas as pd
from pathlib import Path


outdir = "/home/ecoopman/outputresults"
sample_info = pd.read_table("/home/ecoopman/meth_bam_files/samples_mayo.txt").set_index("SAMPLE", drop=False)
sample_info.BAM = "/home/ecoopman/meth_bam_files/" + sample_info.BAM 
sample_info["FILE PRESENT"] = sample_info.BAM.apply(os.path.isfile)
sample_info = sample_info.loc[sample_info["FILE PRESENT"] == True]

sample_info_baminput = pd.read_table("/home/ecoopman/heatmap/baminput/bamfiles/bamfiles.txt").set_index("SAMPLE")
sample_info_baminput.BAM = "/home/ecoopman/heatmap/baminput/bamfiles/" + sample_info_baminput.BAM
sample_info_baminput["FILE PRESENT"] = sample_info_baminput.BAM.apply(os.path.isfile)
sample_info_baminput = sample_info_baminput.loc[sample_info_baminput["FILE PRESENT"] == True]

rule all:
    input:
        os.path.join(outdir, "methylmap/methylmap.tsv"),
        os.path.join(outdir, "methylmap/methylmap.html"),


rule methylmap:
    input:
        #inp = os.path.join("/home/ecoopman/heatmap/methfreqtable/methfreqtable_chr_limited.tsv.gz"),
        inp = expand(os.path.join(outdir, "calculatemethylationfrequency/{id}_phase{phase}.tsv.gz"), id=sample_info.index, phase = [1,2]),
        tbi = expand(os.path.join(outdir, "calculatemethylationfrequency/{id}_phase{phase}.tsv.gz.tbi"), id=sample_info.index, phase = [1,2]),
        #bam = expand(os.path.join("/home/ecoopman/heatmap/baminput/bamfiles/{id}.bam"), id=sample_info_baminput.index),
        #fasta = os.path.join("/home/ecoopman/heatmap/baminput/hg38.fa.gz"),
        gff = os.path.join("/home/ecoopman/ONT-meth-Elise/gff3/gencode_v40_annotation_sorted.gff3.gz"),
    params:
        names = ["sample_01_1","sample_01_2","sample_02_1","sample_02_2","sample_03_1","sample_03_2","sample_04_1","sample_04_2","sample_05_1","sample_05_2", "sample_06_1","sample_06_2","sample_07_1","sample_07_2","sample_08_1","sample_08_2","sample_09_1","sample_09_2","sample_10_1","sample_10_2","sample_11_1","sample_11_2","sample_12_1","sample_12_2","sample_13_1","sample_13_2","sample_14_1","sample_14_2","sample_15_1","sample_15_2","sample_16_1","sample_16_2","sample_17_1","sample_17_2","sample_18_1","sample_18_2","sample_19_1","sample_19_2","sample_20_1","sample_20_2","sample_21_1","sample_21_2","sample_22_1","sample_22_2","sample_23_1","sample_23_2","sample_24_1","sample_24_2","sample_25_1","sample_25_2","sample_26_1","sample_26_2","sample_27_1","sample_27_2","sample_28_1","sample_28_2","sample_29_1","sample_29_2","sample_30_1","sample_30_2","sample_31_1","sample_31_2","sample_32_1","sample_32_2","sample_33_1","sample_33_2","sample_34_1","sample_34_2","sample_35_1","sample_35_2","sample_36_1","sample_36_2","sample_37_1","sample_37_2","sample_38_1","sample_38_2","sample_39_1","sample_39_2","sample_40_1","sample_40_2","sample_41_1","sample_41_2","sample_42_1","sample_42_2","sample_43_1","sample_43_2","sample_44_1","sample_44_2","sample_45_1","sample_45_2","sample_46_1","sample_46_2","sample_47_1","sample_47_2","sample_48_1","sample_48_2","sample_49_1","sample_49_2","sample_50_1","sample_50_2","sample_51_1","sample_51_2","sample_52_1","sample_52_2","sample_53_1","sample_53_2", "sample_54_1","sample_54_2","sample_55_1","sample_55_2","sample_56_1","sample_56_2","sample_57_1","sample_57_2","sample_58_1","sample_58_2","sample_59_1","sample_59_2","sample_60_1","sample_60_2","sample_61_1","sample_61_2","sample_62_1","sample_62_2","sample_63_1","sample_63_2","sample_64_1","sample_64_2","sample_65_1","sample_65_2","sample_66_1","sample_66_2","sample_67_1","sample_67_2","sample_68_1","sample_68_2","sample_69_1","sample_69_2","sample_70_1","sample_70_2","sample_71_1","sample_71_2","sample_72_1","sample_72_2","sample_73_1","sample_73_2","sample_74_1","sample_74_2","sample_75_1","sample_75_2","sample_76_1","sample_76_2","sample_77_1","sample_77_2","sample_78_1","sample_78_2","sample_79_1","sample_79_2","sample_80_1","sample_80_2","sample_81_1","sample_81_2","sample_82_1","sample_82_2","sample_83_1","sample_83_2","sample_84_1","sample_84_2","sample_85_1","sample_85_2","sample_86_1","sample_86_2","sample_87_1","sample_87_2"],
        #groups = ["case", "control", "control", "case"]
    output:
       outtable = os.path.join(outdir, "methylmap/methylmap.tsv"),
       outfig = os.path.join(outdir, "methylmap/methylmap.html"),
    conda:
        "/home/ecoopman/ONT-meth-Elise/envs/methylmap.yml"
    threads: 8
    log:
        os.path.join(outdir, "logs/methylmap.log")
    shell:
        "methylmap --files {input.inp} --window chr20:58839718-58911192 --simplify --names {params.names} --gff {input.gff} --outfig {output.outfig} --outtable {output.outtable} 2> {log}"
