# METHYLMAP

## EXAMPLE

![GNAS methylmap](assets/1000Genomes_GNAS.png)  

Methylmap is a tool for visualization of modified nucleotide frequencies for large cohort sizes, and allows for quick and easy consulting of nucleotide methylation frequencies of individuals in the 1000Genomes ONT project.

You can visualize your own data through the methylmap command line tool (available through bioconda and pypi) and on the methylmap web application <https://methylmap.bioinf.be>.
Consulting the methylation frequencies of the 1000Genomes ONT project can be done through the methylmap web application.

Installation of the methylmap command line tool:

```bash
conda install -c bioconda methylmap
```

Or:

```bash
pip install methylmap
```

If this application is useful for your research, please cite:
[our preprint](https://www.biorxiv.org/content/10.1101/2022.11.28.518239v1) and [the underlying 1000 Genomes Project ONT dataset](https://www.medrxiv.org/content/10.1101/2024.03.05.24303792v1).

### METHYLMAP WEB APPLICATION

#### INPUT POSSIBILITIES

The methylmap web application supports the visualization of own modification frequencies data by uploading a tab separated .tsv/.tsv.gz file. The file should contain the following columns: "chrom", "position", "sample_1", "sample_2", ... "sample_n". Example:

```text
chrom position sample_1 sample_2 sample_3 sample_4
chr1 100000.0 0.000 0.167 0.000 0.077
chr1 100000.5 0.000 0.000 0.100 0.000
chr1 100001.0 0.000 0.000 0.000 0.222
chr1 100002.0 0.000 0.000 0.000 0.000
chr1 100003.0 0.000 0.000 0.000 0.000
```

The maximum file size is 100MB.


#### GENERATING A METHYLATION FREQUENCY TABLE WITH THE MULTIPARSETABLE.PY SCRIPT

The required input table can be generated using the multiparsetable.py script, that supports the following input possibilities:

- BAM/CRAM files with MM and ML tags.

- files from nanopolish (as processed by calculate_methylation_frequency.py). The methylation calls can additionally be phased using the available scripts in the "scripts" folder.

=> The multiparsetable.py script can be found in the "scripts" folder. Example:

```bash
python multiparsetable.py --files cramfileA.cram cramfileB.cram --fasta reference.fa --output methfreqtable.tsv --window chr20:58839718-58911192
python multiparsetable.py --files nanopolishfileA.tsv nanopolishfileB.tsv --output methfreqtable.tsv --window chr20:58839718-58911192 
```

#### ANNOTATION FILES

- currently available annotation files on the methylmap website are:

  - <https://www.gencodegenes.org/human/release_46.html>:   Release 46 (GRCh38.p14) - comprehensive gene annotation

  - <https://www.gencodegenes.org/mouse/>:                  Release M36 (GRCm39) - comprehensive gene annotation

  If you would like to use another annotation file, please upload your request through the Github Issues page.


### METHYLMAP COMMAND LINE TOOL

#### INPUT POSSIBILITIES

- BAM/CRAM files with MM and ML tags. Use --files input option and --fasta for the reference genome.
- files from nanopolish (as processed by calculate_methylation_frequency.py). The methylation calls can additionally be phased using the available scripts in the "scripts" folder. Use --files input option.
- a tab separated table with nucleotide modification frequencies over all positions (methfreqtable), required header names are "chrom" (column with chromosome information) and "position" (columns with position information). Use --table input option. Example:

```text
chrom position sample_1 sample_2 sample_3 sample_4
chr1 100000.0 0.000 0.167 0.000 0.077
chr1 100000.5 0.000 0.000 0.100 0.000
chr1 100001.0 0.000 0.000 0.000 0.222
chr1 100002.0 0.000 0.000 0.000 0.000
chr1 100003.0 0.000 0.000 0.000 0.000
```

- a tab separated file with an overview table containing all nanopolish or BAM/CRAM files and their sample name and experimental group (header requires "file", "name" and "group"). Use --tsv input option. When using BAM/CRAM files, please provide the reference genome with the --fasta argument. Example:

```text
file    name    group
/home/path_to_file/bamfile_sample_1.bam   samplename_1    case
/home/path_to_file/bamfile_sample_2.bam   samplename_2    control
/home/path_to_file/bamfile_sample_3.bam   samplename_3    control
/home/path_to_file/bamfile_sample_4.bam   samplename_4    case
```

#### EXAMPLE COMMAND LINE TOOL USAGE

Example command line tool usage:

```bash
methylmap --files cramfileA.cram cramfileB.cram --fasta reference.fa --gff annotation.gff3.gz --window chr20:58839718-58911192
methylmap --files nanopolishfileA.tsv nanopolishfileB.tsv --gff annotation.gff3.gz --window chr20:58839718-58911192 
methylmap --table methfreqtable.tsv  --gff annotation.gff3.gz --window chr20:58839718-58911192
methylmap --table overviewtable.tsv --fasta reference.fa --gff annotation.gff3.gz --window chr20:58839718-58911192                                        (--fasta argument required when files in overviewtable are BAM/CRAM files)
methylmap --files cramfileA.cram cramfileB.cram --fasta reference.fa --gff annotation.gff3.gz --window chr20:58839718-58911192 --names sampleA sampleB sampleC sampleD --groups case control case control
```

#### IMPORTANT INFORMATION

Important: Adding a GFF3 file is required, use the --gff argument.

- File should be bgzipped
- File should be sorted (use: zcat annotation.gff3.gz  | sort -k1,1V -k4,4n | bgzip > annotation_sorted.gff3.gz)
- File should be indexed (use: tabix -p gff annotation_sorted.gff3.gz)

Important: When using BAM/CRAM files as input, a reference genome with the --fasta argument is required.

- File should not be zipped/bgzipped.
- Indexed file should be present in the same directory as the fasta file.

Important: Individuals with 40% or more missing data are removed. If there are still missing values, they are estimated using the pandas interpolate method. Any remaining missing values after interpolation result in the removal of those individuals.

```text
usage: methylmap [-h] [-f FILES [FILES ...] | -t TABLE] [-w WINDOW] [-n [NAMES ...]] --gff GFF [--output OUTPUT] [--groups [GROUPS ...]] [-s] [--fasta FASTA]
                 [--mod {m,h}] [--hapl] [--dendro] [--threads THREADS] [--quiet] [--debug] [--host HOST] [--port PORT] [-v]

Plotting tool for population scale nucleotide modifications.

options:
  -h, --help            show this help message and exit
  -f FILES [FILES ...], --files FILES [FILES ...]
                        list with BAM/CRAM files or nanopolish (processed with calculate_methylation_frequency.py) files
  -t TABLE, --table TABLE
                        methfreqtable or overviewtable input
  -w WINDOW, --window WINDOW
                        region to visualise, format: chr:start-end (example: chr20:58839718-58911192)
  -n [NAMES ...], --names [NAMES ...]
                        list with sample names
  --gff GFF             add annotation track based on GFF3 file
  --output OUTPUT       TSV file to write the frequencies to.
  --groups [GROUPS ...]
                        list of experimental group for each sample
  -s, --simplify        simplify annotation track to show genes rather than transcripts
  --fasta FASTA         fasta reference file, required when input is BAM/CRAM files or overviewtable with BAM/CRAM files
  --mod {m,h}           modified base of interest when BAM/CRAM files as input. Options are: m, h, default = m
  --hapl                display modification frequencies in input BAM/CRAM file for each haplotype (alternating haplotypes in methylmap)
  --dendro              perform hierarchical clustering on the samples/haplotypes and visualize with dendrogram on sorted heatmap as output
  --threads THREADS     number of threads to use when processing BAM/CRAM files
  --quiet               suppress modkit output
  --debug               Run the app in debug mode
  --host HOST           Host IP used to serve the application
  --port PORT           Port used to serve the application
  -v, --version         print version and exit
```

### MORE INFORMATION

The app.py script is the main script for the methylmap web application. The methylmap.py script is the main script for the methylmap command line tool.

More information: <https://www.biorxiv.org/content/10.1101/2022.11.28.518239v1>
