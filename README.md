## METHYLMAP

Methylmap is a tool for visualization of modified nucleotide frequencies for large cohort sizes. The tools is available through bioconda and pypi, and can be installed using the following commands:
```
conda install -c bioconda methylmap
pip install methylmap
```

or through the methylmap web application at https://methylmap.bioinf.be. The methylmap web application visualizes nucleotide modification frequencies from the 1000Genomes ONT project, and allows users to upload their own data for visualization.
 
If this application is useful for your research, please cite:
https://www.biorxiv.org/content/10.1101/2022.11.28.518239v1 (methylmap)


https://www.medrxiv.org/content/10.1101/2024.03.05.24303792v1 (the underlying 1000Genomes ONT dataset)


### EXAMPLE

![GNAS methylmap](assets/1000Genomes_GNAS.png)  

### METHYLMAP WEB APPLICATION
#### INPUT POSSIBILITIES 

The methylmap web application supports the visualization of own modification frequencies data by uploading a tab separated .tsv file. The file should contain the following columns: "chrom", "position", "sample_1", "sample_2", ... "sample_n". Example:
```
chrom	position	sample_1	sample_2	sample_3	sample_4
chr1	100000.0	0.000	0.167	0.000	0.077
chr1	100000.5	0.000	0.000	0.100	0.000
chr1	100001.0	0.000	0.000	0.000	0.222
chr1	100002.0	0.000	0.000	0.000	0.000
chr1	100003.0	0.000	0.000	0.000	0.000
```

Such a table can be generated using the multiparsetable.py script, that supports the following input possibilities:
- BAM/CRAM files with MM and ML tags. 

- files from nanopolish (as processed by calculate_methylation_frequency.py). The methylation calls can additionally be phased using the available scripts in the "scripts" folder.

#### METHYLMAP COMMAND LINE TOOL
#### INPUT POSSIBILITIES
- BAM/CRAM files with MM and ML tags. Use --files input option.
- files from nanopolish (as processed by calculate_methylation_frequency.py). The methylation calls can additionally be phased using the available scripts in the "scripts" folder. Use --files input option.
- an own tab separtated table with nucleotide modification frequencies over all positions (methfreqtable), required header names are "chrom" (column with chromosome information) and "position" (columns with position information). Use --table input option. Example:
```
chrom	position	sample_1	sample_2	sample_3	sample_4
chr1	100000.0	0.000	0.167	0.000	0.077
chr1	100000.5	0.000	0.000	0.100	0.000
chr1	100001.0	0.000	0.000	0.000	0.222
chr1	100002.0	0.000	0.000	0.000	0.000
chr1	100003.0	0.000	0.000	0.000	0.000
```
- a tab separated file with an overview table containing all nanopolish or BAM/CRAM files and their sample name and experimental group (header requires "path", "name" and "group"). Use --table input option. Example:
```
path    name    group
/home/path_to_file/bamfile_sample_1.bam   samplename_1    case
/home/path_to_file/bamfile_sample_2.bam   samplename_2    control
/home/path_to_file/bamfile_sample_3.bam   samplename_3    control
/home/path_to_file/bamfile_sample_4.bam   samplename_4    case
```

#### USAGE
```
usage: methylmap [-h] [-f FILES [FILES ...] | -t TABLE] [-w WINDOW]
                 [-n [NAMES ...]] --gff GFF [--output OUTPUT]
                 [--groups [GROUPS ...]] [-s] [--fasta FASTA]
                 [--mod {m,h}] [--hapl] [--dendro]
                 [--threads THREADS] --db DB [--quiet] [--debug]
                 [--host HOST] [--port PORT] [-v]
```

### MORE INFORMATION

More information: https://www.biorxiv.org/content/10.1101/2022.11.28.518239v1