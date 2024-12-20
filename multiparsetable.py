import os
import sys
import gzip
import shlex
from pathlib import Path
import subprocess
import numpy as np
import logging
import itertools
import tempfile
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ThreadPoolExecutor
from argparse import ArgumentParser


def get_args():
    parser = ArgumentParser(
        description="Plotting tool for population scale nucleotide modifications."
    )
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument(
        "-f",
        "--files",
        nargs="+",
        help="list of CRAM or BAM files",
    )
    action.add_argument("--tsv", help="TSV file with file name, sample name and group")
    parser.add_argument(
        "-w",
        "--window",
        help="region to visualise, format: chr:start-end (example: chr20:58839718-58911192)",
    )
    parser.add_argument(
        "-n", "--names", nargs="*", default=[], help="list with sample names"
    )
    parser.add_argument(
        "--expand",
        help="number of base pairs to expand the window with in both directions",
        type=int,
        default=0,
    )
    parser.add_argument(
        "-o", "--output", help="TSV file to write the frequencies to", required=True
    )
    parser.add_argument(
        "--groups", nargs="*", help="list of experimental group for each sample"
    )
    parser.add_argument(
        "--fasta",
        help="fasta reference file, required when input is BAM/CRAM files or overviewtable with BAM/CRAM files",
    )
    parser.add_argument(
        "--mod",
        help="modified base of interest when BAM/CRAM files as input. Options are: m, h, default = m",
        default="m",
        choices=["m", "h"],  # methylation or hydroxymethylation
    )
    parser.add_argument(
        "--hapl",
        action="store_true",
        help="display modification frequencies in input BAM/CRAM file for each haplotype (alternating haplotypes in methylmap)",
    )
    parser.add_argument(
        "--threads",
        help="number of threads to use when processing BAM/CRAM files",
        type=int,
        default=12,
    )
    parser.add_argument("--quiet", action="store_true", help="suppress modkit output")
    args = parser.parse_args()
    if args.files:
        if len(args.names) == 0:
            for b in args.files:
                c = b.rsplit(".", 2)[0]
                args.names.append(os.path.basename(c))
        else:
            if len(args.files) != len(args.names):
                sys.exit(
                    f"ERROR: expecting the same number of input files [{len(args.files)}] and names [{len(args.names)}]"
                )
    if args.files or args.tsv:
        if not args.window:
            sys.exit("ERROR: please provide a genomic region with --window")
    if args.tsv:
        if not Path(args.tsv).is_file():
            sys.exit(f"ERROR: file {args.tsv} does not exist, please check the path!")
        tsv_df = pd.read_csv(args.tsv, sep="\t")
        if not all(col in tsv_df.columns for col in ["file", "name"]):
            sys.exit("ERROR: provide a TSV file with minimally columns 'file' and 'name'")
        args.files = tsv_df["file"].tolist()
        args.names = tsv_df["name"].tolist()
        if "group" in tsv_df.columns:
            args.groups = tsv_df["group"].tolist()
    if args.files:
        first_file = args.files[0]
        if first_file.endswith(".bam") or first_file.endswith(".cram"):
            if not args.fasta:
                sys.exit("ERROR: please provide a reference fasta file with --fasta")
    return args


def main():
    args = get_args()
    window = Region(args.window, args.expand)
    file_type = file_sniffer(args.files[0])
    try:
        if file_type == "nanopolish_calc_meth_freq":
            overviewtable = parse_nanopolish(args, window)
        elif file_type in ["cram", "bam"]:
            check_modkit()
            overviewtable = parse_bam(args, window)
        if args.output:
            overviewtable.to_csv(args.output, sep="\t", na_rep='NA', header=True)
        else:
            print(overviewtable.to_csv(sep="\t", na_rep='NA', header=True))
    except Exception as e:
        logging.error("Error processing input file(s).")
        logging.error(e, exc_info=True)
        sys.stderr.write("\n\n\nError processing input file(s)!\n")
        raise


def check_modkit():
    rc = subprocess.call(["which", "modkit"], stdout=subprocess.DEVNULL)
    if rc != 0:
        sys.exit("\n\nIs modkit installed? See github.com/nanoporetech/modkit")


def parse_bam(args, window):
    args_list = [
        (file, name, args, window) for file, name in zip(args.files, args.names)
    ]
    # Process files concurrently with a maximum of 12 threads
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        results = list(
            tqdm(executor.map(process_single_file, args_list), total=len(args_list))
        )

    dfs = list(itertools.chain(*results))

    methfrequencytable = dfs[0].join(dfs[1:], how="outer")

    if len(methfrequencytable) == 0:
        err = "WARNING: length of methylation frequency table is zero. Do the input files contain data?"
        logging.error(err)
        sys.exit(err)

    if args.groups:
        logging.info("Sort columns of methfrequencytable based on group")
        if args.hapl:
            args.groups = list(itertools.chain(*zip(args.groups, args.groups)))
        headerlist = list(methfrequencytable.columns.values)
        if len(headerlist) == len(args.groups):
            res = zip(headerlist, args.groups)
            output = sorted(list(res), key=lambda x: x[1])
            orderedlist = [i[0] for i in output]
            methfrequencytable = methfrequencytable.reindex(columns=orderedlist)
        else:
            sys.exit(
                f"ERROR when matching --groups with samples, is length of --groups list ({len(args.groups)}) matching with number of sample files?"
            )

    return methfrequencytable


def process_single_file(function_args):
    input, name, args, window = function_args

    with tempfile.TemporaryDirectory() as temp_dir:
        file_basename = os.path.basename(input)
        logging.info("Extract mod frequencies in bam/cram using modkit")
        # when processing haplotypes, <output> is a directory
        output = (
            temp_dir if args.hapl else os.path.join(temp_dir, f"{file_basename}.bed")
        )
        logging.info(f"Processing file {input}")
        stderr = run_modkit_pileup(
            input=input,
            output=output,
            ignore="m" if args.mod == "h" else "h",
            window=window,
            fasta=args.fasta,
            log_file=os.path.join(temp_dir, "modkit.log"),
            haplotype=args.hapl,
        )
        if not args.quiet:
            tqdm.write(stderr, file=sys.stderr)
        logging.info(f"Read {input} in a dataframe per haplotype.")
        if args.hapl:
            return [
                process_modkit_tsv(
                    f"{output}/H_{haplotype}.bed", name=f"{name}_{haplotype}"
                )
                for haplotype in [1, 2]
            ]
        else:
            return [process_modkit_tsv(output, name=name)]


def run_modkit_pileup(input, output, ignore, window, fasta, log_file, haplotype=False):
    cmd = f"modkit pileup {input} {output} --ignore {ignore} --region={window.fmt} --cpg --ref={fasta} --log-filepath {log_file} --only-tabs"
    if haplotype:
        cmd += " --partition-tag HP --prefix H"
    try:
        modkit = subprocess.Popen(
            shlex.split(cmd),
            stderr=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
        )
        stderr = modkit.stderr.read()
    except Exception as e:
        logging.error(e, exc_info=True)
        sys.stderr.write("\n\nError parsing bam/cram file with modkit.\n")
        sys.stderr.write(f"\n\n\nDetailed error: {modkit.stderr.read()}\n")
        raise
    if modkit.returncode:
        sys.exit(f"Received modkit error:\n{modkit.stderr.read()}\n")
    return stderr


def process_modkit_tsv(filename, name):
    headerlist = [
        "chrom",
        "position",
        "endposition",
        "modified_base_code",
        "score",
        "strand",
        "startposition2",
        "endposition2",
        "color",
        "Nvalid_cov",
        "fractionmodified",
        "Nmod",
        "Ncanonical",
        "Nother_mod",
        "Ndelete",
        "Nfail",
        "Ndiff",
        "Nnocall",
    ]
    if not Path(filename).is_file():
        sys.exit(f"\n\nERROR: File {filename} for {name} does not exist, please check the path and log!\n")
    logging.info(f"Reading the modkit file for {name} in a dataframe.")
    df = pd.read_table(
        filename,
        sep="\t",
        header=None,
        names=headerlist,
        usecols=["chrom", "position", "fractionmodified"],
    )

    return df.set_index(["chrom", "position"]).rename(
        columns={"fractionmodified": name}
    )


def parse_nanopolish(args, window):
    """
    Converts a file from nanopolish to a pandas dataframe
    input can be from calculate_methylation_frequency.py or overviewtable with these files
    which will return a dataframe with 'chromosome', 'pos', 'methylated_frequency'.
    """
    dfs = []
    for file, name in zip(args.files, args.names):
        if window:
            if not Path(file + ".tbi").is_file():
                logging.info(
                    "Make tabix file of input files for fast selection of window of interest"
                )
                try:
                    make_tabix = subprocess.Popen(
                        shlex.split(f"tabix -S1 -s1 -b2 -e3 {file}"),
                        stderr=subprocess.PIPE,
                    )
                except FileNotFoundError as e:
                    logging.error("Error when making a .tbi file.")
                    logging.error(e, exc_info=True)
                    sys.stderr.write("\n\nERROR when making a .tbi file.\n")
                    sys.stderr.write("Is tabix installed and on the PATH?\n.")
                    sys.stderr.write(
                        f"\n\n\nDetailed error:\n\n{make_tabix.stderr.read()}\n"
                    )
                    raise
                if make_tabix.returncode:
                    sys.exit(
                        f"\n\n\nReceived tabix error\n\n{make_tabix.stderr.read()}"
                    )

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
                sys.stderr.write("Is tabix installed and on the PATH?\n")
                sys.stderr.write(
                    f"\n\n\nDetailed error: {tabix_stream.stderr.read()}\n"
                )
                raise
            if tabix_stream.returncode:
                sys.exit(f"\n\n\nReceived tabix error\n\n{tabix_stream.stderr.read()}")

            header = gzip.open(file, "rt").readline().rstrip().split("\t")
            table = pd.read_csv(
                tabix_stream.stdout, sep="\t", header=None, names=header
            )
            logging.info("Read the file in a dataframe.")
        else:
            table = pd.read_csv(file, sep="\t")
        logging.info("Read the file in a dataframe.")
        table.drop(
            [
                "group_sequence",
                "called_sites_methylated",
                "num_motifs_in_group",
                "called_sites",
            ],
            axis=1,
            inplace=True,
        )
        table["midpoint"] = (table["start"] + table["end"]) / 2
        table["position"] = table["chromosome"] + ":" + table["midpoint"].astype(str)
        table = table.set_index("position").drop(
            ["chromosome", "start", "end", "midpoint"], axis=1
        )
        dfs.append(table.rename(columns={"methylated_frequency": name}))
    modfrequencytable = dfs[0].join(dfs[1:], how="outer")
    modfrequencytable.reset_index(inplace=True)
    modfrequencytable["chrom"] = modfrequencytable["position"].str.split(
        ":", expand=True
    )[0]
    modfrequencytable["position"] = modfrequencytable["position"].str.split(
        ":", expand=True
    )[1]
    modfrequencytable.set_index(["chrom", "position"], inplace=True)
    modfrequencytable.reset_index(inplace=True)
    if len(modfrequencytable) == 0:
        err = "WARNING: length of modification frequency table is zero. Do the input files contain data?"
        logging.error(err)
        sys.exit(err)

    if args.files:
        if args.groups:
            if args.dendro:
                err = "Columns will not be sorted based on --group input since hierarchical clustering with --dendro is requested."
                logging.warning(err)
                sys.stderr.write(err)
            else:
                logging.info("Sort columns of modfrequencytable based on group")
                headerlist = list(modfrequencytable.columns.values)
                if len(headerlist) == len(args.groups):
                    res = zip(headerlist, args.groups)
                    output = sorted(list(res), key=lambda x: x[1])
                    orderedlist = [i[0] for i in output]
                    modfrequencytable = modfrequencytable.reindex(columns=orderedlist)
                else:
                    sys.exit(
                        f"ERROR when matching --groups with samples, is length of --groups list ({len(args.groups)}) matching with number of sample files?"
                    )
    modfrequencytable.set_index(["chrom"], inplace=True)
    modfrequencytable["position"] = (
        modfrequencytable["position"].astype(float).astype(int)
    )
    return modfrequencytable


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
    # input: calculate_methylation_frequency.py output (.tsv or .tsv.gz) OR own modfreqtable (.tsv or .tsv.gz)
    if is_gz_file(filename):
        import gzip

        header = gzip.open(filename, "rt").readline()
    else:
        header = open(filename, "r").readline()

    if "methylated_frequency" in header:
        # calculate_methylation_frequency.py output as input
        return "nanopolish_calc_meth_freq"
    sys.exit(f"\n\n\nInput file {filename} not recognized!\n")


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


def is_gz_file(filepath):
    import binascii

    with open(filepath, "rb") as test_f:
        return binascii.hexlify(test_f.read(2)) == b"1f8b"


class Region(object):
    def __init__(self, region=None, expand=False, tup=None):
        if tup:
            self.chromosome, self.begin, self.end = tup
        else:
            try:
                self.chromosome, interval = region.replace(",", "").split(":")
                try:
                    # see if just integer chromosomes are used
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
        self.start = self.begin  # start as alias for begin
        self.size = self.end - self.begin
        if self.size < 0:
            sys.exit(
                "\n\nERROR: Window (-w/--window) inproperly formatted, "
                "begin of the interval has to be smaller than end\n\n"
            )
        self.string = f"{self.chromosome}_{self.begin}_{self.end}"
        self.fmt = f"{self.chromosome}:{self.begin}-{self.end}"

    def __mul__(self, other):
        new_half_size = int((self.size * other) / 2)
        middlepoint = (self.begin + self.end) / 2
        return Region(
            tup=(
                self.chromosome,
                int(middlepoint - new_half_size),
                int(middlepoint + new_half_size),
            )
        )

    def __truediv__(self, other):
        new_half_size = int((self.size / other) / 2)
        middlepoint = (self.begin + self.end) / 2
        return Region(
            tup=(
                self.chromosome,
                int(middlepoint - new_half_size),
                int(middlepoint + new_half_size),
            )
        )


if __name__ == "__main__":
    main()
