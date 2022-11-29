### METHYLMAP

Methylmap is a tool for visualization of modified nucleotide frequencies for large cohort sizes. Supported input possibilities are:
- BAM/CRAM files with MM and ML tags
- files from nanopolish (as processed by calculate_methylation_frequency.py). The methylation calls can additionally be phased using the available scripts in the "scripts" folder.
- an own tab separtated table with nucleotide modification frequencies over all positions, required header names are "chrom" (column with chromosome information) and "position" (columns with position information). Example:
```
	chrom	position	sample_1	sample_2	sample_3	sample_4
0	chr1	100000.0	0.000	0.167	0.000	0.077
1	chr1	100000.5	0.000	0.000	0.100	0.000
2	chr1	100001.0	0.000	0.000	0.000	0.222
3	chr1	100002.0	0.000	0.000	0.000	0.000
4	chr1	100003.0	0.000	0.000	0.000	0.000
```
- a tab separated file with an overview table containing all nanopolish or BAM/CRAM files and their sample name and experimental group (header requires "path", "name" and "group"). Example:
```
        path    name    group
0     /home/path_to_file/bamfile_sample_1.bam   samplename_1    case
1     /home/path_to_file/bamfile_sample_2.bam   samplename_2    control
2     /home/path_to_file/bamfile_sample_3.bam   samplename_3    control
3     /home/path_to_file/bamfile_sample_4.bam   samplename_4    case
````

## INSTALLATION
pip install methylmap

## USAGE

```
usage: methylmap [-h] (-f FILES [FILES ...] | -t TABLE) [-w WINDOW] [-n [NAMES ...]] [--gff GFF] [--expand EXPAND] [--outtable OUTTABLE] [--outfig OUTFIG] [--groups [GROUPS ...]] [-s] [--fasta FASTA]
                 [--mod {5mC,5hmC,5fC,5caC,5hmU,5fU,5caU,6mA,5oxoG,Xao}] [-v]

Plotting tool for population scale nucleotide modifications.

options:
  -h, --help            show this help message and exit
  -f FILES [FILES ...], --files FILES [FILES ...]
                        Nanopolish calculate_methylation_frequency.py output or BAM/CRAM files.
  -t TABLE, --table TABLE
                        Methfrequencytable or overviewtable input.
  -w WINDOW, --window WINDOW
                        Region to visualise. Format: chr:start-end (Example: chr20:58839718-58911192)
  -n [NAMES ...], --names [NAMES ...]
                        List with sample names.
  --gff GFF, --gtf GFF  Add annotation track based on GTF/GFF file.
  --expand EXPAND       Number of base pairs to expand the window with in both directions.
  --outtable OUTTABLE   File to write the frequencies table to.
  --outfig OUTFIG       File to write output heatmap (in HTML format) to.
  --groups [GROUPS ...]
                        List of experimental group for each sample.
  -s, --simplify        Simplify annotation track to show genes rather than transcripts.
  --fasta FASTA         Fasta reference file, required when input is BAM/CRAM files or overviewtable with BAM/CRAM files.
  --mod {5mC,5hmC,5fC,5caC,5hmU,5fU,5caU,6mA,5oxoG,Xao}
                        Modified base of interest when BAM/CRAM files as input. Options are: 5mC, 5hmC, 5fC, 5caC, 5hmU, 5fU, 5caU, 6mA, 5oxoG, Xao. Default = 5mC
  -v, --version         Print version and exit.
  ```

  