import shlex
import sys
import gzip
import logging
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path


class Region(object):
    def __init__(self, region, expand=False):
        try:
            self.chromosome, interval = region.replace(",", "").split(":")
            try:
                #see if just integer chromosomes are used
                self.chromosome = int(self.chromosome)
            except ValueError:
                pass
            self.begin, self.end = [int(i) for i in interval.split("-")]
        except ValueError:
            sys.exit(
                "\n\nERROR: Window (-w/--window) inproperly formatted, "
                "an example of accepted formats is:\n'chr5:150200605-150423790'\n\n"
            )
        if expand:
            self.begin = self.begin - int(expand)
            self.end = self.end + int(expand)
        self.start = self.begin #start as alias for begin
        self.size = self.end - self.begin
        if self.size < 0:
            sys.exit(
                "\n\nERROR: Window (-w/--window) inproperly formatted, "
                "begin of the interval has to be smaller than end\n\n"
            )
        self.string = f"{self.chromosome}_{self.begin}_{self.end}"
        self.fmt = f"{self.chromosome}:{self.begin}-{self.end}"


def read_mods(files, table, names, window, groups, gff, fasta, mod, dendro):
    """
    Deciding of input file(s) type and processing them.
    """
    if files:
        file_type = file_sniffer(files[0])
    elif table:
        file_type = file_sniffer(table)
    try:
        if file_type == "nanopolish_calc_meth_freq":
            return parse_nanopolish(files, table, names, window, groups, dendro)
        elif file_type == "overviewtable_nanopolishfiles":
            return parse_nanopolish(files, table, names, window, groups, dendro)
        elif file_type == "methfrequencytable":
            return parse_methfrequencytable(table, names, window, groups, gff, dendro)
        elif file_type in ["cram", "bam"]:
            rc = subprocess.call(['which', 'modbam2bed'])
            if not rc == 0:
                sys.exit(f"\n\n\nIs modbam2bed installed? Installation: mamba install -c epi2melabs modbam2bed")
            else:
                return parse_bam(files, table, names, window, groups, fasta, mod, dendro)
        elif file_type in ["overviewtable_bam", "overviewtable_cram"]:
            rc = subprocess.call(['which', 'modbam2bed'])
            if not rc == 0:
                sys.exit(f"\n\n\nIs modbam2bed installed? Instalation: mamba install -c epi2melabs modbam2bed")
            else:
                return parse_bam(files, table, names, window, groups, fasta, mod, dendro)
    except Exception as e:
        logging.error("Error processing input file(s).")
        logging.error(e, exc_info=True)
        sys.stderr.write("\n\n\nError processing input file(s)!\n")
        raise


def parse_overviewtable(table):
    """
    Parses overviewtable, sorts files by group and makes list of all files and sample names.
    """
    overviewtable = pd.read_table(table).sort_values(["group", "name"])
    files = overviewtable["path"].tolist()
    names = overviewtable["name"].tolist()
    return files, names


def parse_nanopolish(files, table, names, window, groups, dendro):
    """
    Converts a file from nanopolish to a pandas dataframe
    input can be from calculate_methylation_frequency or overviewtable with these files
    which will return a dataframe with 'chromosome', 'pos', 'methylated_frequency'.
    """
    if table:
        logging.info("Extract files and names from overviewtable")
        files, names = parse_overviewtable(table)

    dfs = []
    for file, name in zip(files, names):
        if window:
            if not Path(file + ".tbi").is_file():
                logging.info(
                    "Make tabix file of input files for fast selection of window of interest"
                )
                try:
                    make_tabix = subprocess.Popen(
                        shlex.split(f"tabix -S1 -s1 -b2 -e3 {file}"), stderr=subprocess.PIPE
                    )
                except FileNotFoundError as e:
                    logging.error("Error when making a .tbi file.")
                    logging.error(e, exc_info=True)
                    sys.stderr.write("\n\nERROR when making a .tbi file.\n")
                    sys.stderr.write("Is tabix installed and on the PATH?.")
                    sys.stderr.write(f"\n\n\nDetailed error:\n\n{make_tabix.stderr.read()}\n")
                    raise
                if make_tabix.returncode:
                    sys.exit(f"\n\n\nReceived tabix error\n\n{make_tabix.stderr.read()}")

            try:
                logging.info(f"Reading {file} using a tabix stream.")
                tabix_stream = subprocess.Popen(
                    shlex.split(f"tabix {file} {window.fmt}"),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                )
            except FileNotFoundError as e:
                logging.error("Error when opening a tabix stream.")
                logging.error(e, exc_info=True)
                sys.stderr.write("\n\nERROR when opening a tabix stream.\n")
                sys.stderr.write("Is tabix installed and on the PATH?")
                sys.stderr.write(f"\n\n\nDetailed error: {tabix_stream.stderr.read()}\n")
                raise
            if tabix_stream.returncode:
                sys.exit(f"\n\n\nReceived tabix error\n\n{tabix_stream.stderr.read()}")

            header = gzip.open(file, "rt").readline().rstrip().split("\t")
            table = pd.read_csv(tabix_stream.stdout, sep="\t", header=None, names=header)
            logging.info("Read the file in a dataframe.")
        else:
            table = pd.read_csv(file, sep="\t")
            logging.info("Read the file in a dataframe.")
        table.drop(
            ["group_sequence", "called_sites_methylated", "num_motifs_in_group", "called_sites"],
            axis=1,
            inplace=True,
        )
        table["position"] = (table["start"] + table["end"]) / 2
        table = table.set_index("position").drop(["chromosome", "start", "end"], axis=1)
        dfs.append(table.rename(columns={"methylated_frequency": name}))
    methfrequencytable = dfs[0].join(dfs[1:], how="outer")

    if len(methfrequencytable) == 0:
        logging.error(
            f"WARNING: length of methylation frequency table is zero. Do the input files contain data?"
        )
        sys.exit(
            f"WARNING: length of methylation frequency table is zero. Do the input files contain data?"
        )

    methfreqtable = methfrequencytable.sort_values("position", ascending=True)
    # output is a meth frequency table with position as index and for each sample a column with all the methylation frequencies

    if files:
        if groups:
            if dendro:
                logging.warning(
                    "Columns will not be sorted based on --group input since hierarchical clustering with --dendro is requested."
                )
                sys.stderr.write(
                    "Columns will not be sorted based on --group input since hierarchical clustering with --dendro is requested."
                )
            else:
                logging.info("Sort columns of methfrequencytable based on group")
                headerlist = list(methfreqtable.columns.values)
                if len(headerlist) == len(groups):
                    res = zip(headerlist, groups)
                    output = sorted(list(res), key=lambda x: x[1])
                    orderedlist = [i[0] for i in output]
                    methfreqtable = methfreqtable.reindex(columns=orderedlist)
                else:
                    sys.exit(
                        f"ERROR when matching --groups with samples, is length of --groups list ({len(groups)}) matching with number of sample files?"
                    )

    return methfreqtable, window


def parse_methfrequencytable(table, names, window, groups, gff, dendro):
    """
    Parsing methfrequencytable input.
    """
    df = pd.read_table(table).sort_values("position", ascending=True)

    if window:
        logging.info("Select window out of methfreqtable")
        df = df[df["position"].between(window.begin, window.end)]

    else:
        if gff:
            logging.info(
                "If no window given and annotation requested, take window out of methfreqtable."
            )
            if len(df["chrom"].unique()) > 1:
                sys.exit(
                    "\n\nError when extracting window out of methfreqtable positions. Chrom column can not contain more than one chromosome \n"
                )
            chrom = df.iloc[0, df.columns.get_loc("chrom")]
            if chrom.startswith("chr"):  # if in format "chr1"
                chrom = chrom
            else:  # if in format "1"
                chrom = "chr" + chrom
            numberofpositions = len(df) - 1
            begin = float(df["position"].iat[0])
            end = float(df["position"].iat[numberofpositions])
            window = Region(f"{chrom}:{round(begin)}-{round(end)}")

    df.drop(["chrom"], axis=1, inplace=True)
    df.set_index("position", inplace=True)

    if len(df) == 0:
        logging.error(
            f"WARNING: length of methylation frequency table is zero. Does the input table contain data?"
        )
        sys.exit(
            f"WARNING: length of methylation frequency table is zero. Does the input table contain data?"
        )

    if names:
        logging.info("Changes sample names to names from --names input")
        df.columns = names

    if groups:
        logging.info("Sort columns of methfrequencytable based on group")
        if dendro:
                logging.warning(
                    "Columns will not be sorted based on --group input since hierarchical clustering with --dendro is requested."
                )
                sys.stderr.write(
                    "Columns will not be sorted based on --group input since hierarchical clustering with --dendro is requested."
                )
        else:
            headerlist = list(df.columns.values)
            if len(headerlist) == len(groups):
                res = zip(headerlist, groups)
                output = sorted(list(res), key=lambda x: x[-1])
                orderedlist = [i[0] for i in output]
                df = df.reindex(columns=orderedlist)
            else:
                logging.error(
                    f"Error when matching --groups with samples. Is length of the --groups list ({len(groups)}) matching with number of samples in table?"
                )
                sys.exit(
                    f"Error when matching --groups with samples. Is length of the --groups list ({len(groups)}) matching with number of samples in table?"
                )

    return df, window


def parse_bam(files, table, names, window, groups, fasta, mod, dendro):
    """
    Converts a bam/cram file to a pandas dataframe.
    """
    if not fasta:
        logging.info("Stop script when no --fasta input")
        sys.exit(
            f"ERROR when parsing bam/cram file, can not find fasta file. Is fasta file given with --fasta argument?"
        )

    if table:
        logging.info("Extract files and names from overviewtable")
        files, names = parse_overviewtable(table)

    dfs = []
    for file, name in zip(files, names):
        try:
            logging.info(
                "Extract modification frequencies from bam/cram files using modbam2bed tool"
            )
            modbam_stream = subprocess.Popen(
                shlex.split(
                    f"modbam2bed --mod_base={mod} --region={window.fmt} --cpg {fasta} {file}"
                ),
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except FileNotFoundError as e:
            logging.error(e, exc_info=True)
            sys.stderr.write("\n\nError when making bedfile of bam/cram file with modbam2bed.\n")
            sys.stderr.write(
                "Is modbam2bed (conda install -c epi2melabs modbam2bed) installed and on the PATH?"
            )
            sys.stderr.write(f"\n\n\nDetailed error: {modbam_stream.stderr.read()}\n")
            raise
        if modbam_stream.returncode:
            sys.exit(f"Received modbam2bed error:\n{modbam_stream.stderr.read()}\n")
        headerlist = [
            "chromosome",
            "start",
            "position",
            "modification",
            "score",
            "strand",
            "start2",
            "end2",
            "column1",
            "read_coverage",
            "methylated_frequency",
        ]

        table = pd.read_csv(modbam_stream.stdout, sep="\t", header=None,names=headerlist, usecols=["position","methylated_frequency"])
        logging.info("Read the file in a dataframe.")
        table = table.set_index("position")
        dfs.append(table.rename(columns={"methylated_frequency": name}))

    methfrequencytable = dfs[0].join(dfs[1:], how="outer")

    if len(methfrequencytable) == 0:
        logging.error(
            f"WARNING: length of methylation frequency table is zero. Do the input files contain data?"
        )
        sys.exit(
            f"WARNING: length of methylation frequency table is zero. Do the input files contain data?"
        )

    methfreqtable = methfrequencytable.sort_values("position", ascending=True)

    if files:
        logging.info("Sort columns of methfrequencytable based on group")
        if groups:
            if dendro:
                logging.warning(
                    "Columns will not be sorted based on --group input since hierarchical clustering with --dendro is requested."
                )
                sys.stderr.write(
                    "Columns will not be sorted based on --group input since hierarchical clustering with --dendro is requested."
                )
            else:
                headerlist = list(methfreqtable.columns.values)
                if len(headerlist) == len(groups):
                    res = zip(headerlist, groups)
                    output = sorted(list(res), key=lambda x: x[1])
                    orderedlist = [i[0] for i in output]
                    methfreqtable = methfreqtable.reindex(columns=orderedlist)
                else:
                    sys.exit(
                        f"ERROR when matching --groups with samples, is length of --groups list ({len(groups)}) matching with number of sample files?")
                    
    return methfreqtable, window


def file_sniffer(filename):
    """
    Takes in a filename and tries to guess the input file type.
    """
    if not Path(filename).is_file():
        sys.exit(f"\n\nERROR: File {filename} does not exist, please check the path!\n")
    if is_bam_file(filename):  # input BAM
        return "bam"
    if is_cram_file(filename):  # input CRAM
        return "cram"
    # input: calculate_methylation_frequency.py output (.tsv or .tsv.gz) OR own methfreqtable (.tsv or .tsv.gz)
    if is_gz_file(filename):
        import gzip

        header = gzip.open(filename, "rt").readline()
    else:
        header = open(filename, "r").readline()

    if "methylated_frequency" in header:
        # calculate_methylation_frequency.py output as input
        return "nanopolish_calc_meth_freq"
    if "path" in header:  # overviewtable
        df = pd.read_table(filename)
        if df["path"].iloc[0].endswith(".tsv"):  ###only checks first file in overviewtable
            return "overviewtable_nanopolishfiles"
        elif df["path"].iloc[0].endswith(".tsv.gz"):  ###only checks first file in overviewtable
            return "overviewtable_nanopolishfiles"
        elif df["path"].iloc[0].endswith(".bam"):  ####only checks first file in overviewtable
            return "overviewtable_bam"
        elif df["path"].iloc[0].endswith(".cram"):  ####only checks first file in overviewtable
            return "overviewtable_cram"
    if header.startswith("chrom"):
        # own methfrequencytable as input: needs first column header to be "chrom"
        return "methfrequencytable"
    sys.exit(f"\n\n\nInput file {filename} not recognized!\n")


def is_gz_file(filepath):
    import binascii

    with open(filepath, "rb") as test_f:
        return binascii.hexlify(test_f.read(2)) == b"1f8b"


def is_cram_file(filepath):
    with open(filepath, "rb") as test_f:
        return test_f.read(4) == b"CRAM"


def is_bam_file(filepath):
    import gzip

    try:
        with gzip.open(filepath) as test_f:
            return test_f.read(3) == b"BAM"
    except OSError:
        return False
