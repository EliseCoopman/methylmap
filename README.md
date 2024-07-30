## METHYLMAP

Repository for the methylmap web application at https://methylmap.bioinf.be. 
The methylmap web application visualizes nucleotide modification frequencies from the 1000Genomes ONT project, and allows users to upload their own data for visualization.
 
If this application is useful for your research, please cite:
https://www.biorxiv.org/content/10.1101/2022.11.28.518239v1 (methylmap)


https://www.medrxiv.org/content/10.1101/2024.03.05.24303792v1 (the underlying 1000Genomes ONT dataset)


### EXAMPLE

![GNAS methylmap](example/20221213182515.png)  

### INPUT POSSIBILITIES

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


### MORE INFORMATION

More information: https://www.biorxiv.org/content/10.1101/2022.11.28.518239v1