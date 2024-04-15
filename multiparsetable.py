import os
import sys
import shlex
from pathlib import Path
import subprocess
import numpy as np
import logging
import itertools
import tempfile
import pandas as pd
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
        help="list with BAM/CRAM files or nanopolish (processed with calculate_methylation_frequency.py) files",
    )
    parser.add_argument(
        "-w",
        "--window",
        help="region to visualise, format: chr:start-end (example: chr20:58839718-58911192)",
    )
    parser.add_argument("--outtable", help="File to write the frequencies table to.")
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
        "--outtable", help="file to write the frequencies table to in tsv format"
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
        choices=["m", "h"],  # methylation  # hydroxymethylation
    )
    parser.add_argument(
        "--hapl",
        action="store_true",
        help="display modification frequencies in input BAM/CRAM file for each haplotype (alternating haplotypes in methylmap)",
    )
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
    return args


def main():
    args = get_args()
    overviewtable = read_mods(
        args.files,
        args.names,
        args.window,
        args.groups,
        args.fasta,
        args.mod,
        args.hapl,
        args.expand,
    )
    overviewtable.to_csv(args.outtable, sep="\t", na_rep=np.NaN, header=True)


def read_mods(files, names, window, groups, fasta, mod, hapl, expand=False):
    window = Region(window, expand)
    """
    Deciding of input file(s) type and processing them.
    """
    if files:
        file_type = file_sniffer(files[0])
    try:
        if file_type in ["cram", "bam"]:
            rc = subprocess.call(["which", "modkit"])
            if not rc == 0:
                sys.exit(
                    "\n\n\nIs modkit installed? Installation: see https://github.com/nanoporetech/modkit"
                )
            else:
                return parse_bam(files, names, window, groups, fasta, mod, hapl)
    except Exception as e:
        logging.error("Error processing input file(s).")
        logging.error(e, exc_info=True)
        sys.stderr.write("\n\n\nError processing input file(s)!\n")
        raise


def parse_bam(files, names, window, groups, fasta, mod, hapl):
    if not fasta:
        logging.info("Stop script when no --fasta input")
        sys.exit(
            "ERROR when parsing bam/cram file, can not find fasta file. Is fasta file given with --fasta argument?"
        )
    with ThreadPoolExecutor(max_workers=12) as executor:
        args_list = [
            (file, name, fasta, mod, window, hapl) for file, name in zip(files, names)
        ]

        # Process files concurrently with a maximum of 12 threads
        results = list(executor.map(process_single_file, args_list))

    dfs = list(itertools.chain(*results))

    methfrequencytable = dfs[0].join(dfs[1:], how="outer")

    if len(methfrequencytable) == 0:
        logging.error(
            "WARNING: length of methylation frequency table is zero. Do the input files contain data?"
        )
        sys.exit(
            "WARNING: length of methylation frequency table is zero. Do the input files contain data?"
        )

    methfreqtable = methfrequencytable.sort_values("position", ascending=True)

    if groups:
        logging.info("Sort columns of methfrequencytable based on group")
        if hapl:
            groupshapl = list(itertools.chain(*zip(groups, groups)))
            groups = groupshapl
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

    return methfreqtable


def process_single_file(args):
    file, name, fasta, mod, window, hapl = args
    dfs = []
    dfs_1 = []
    dfs_2 = []
    headerlist = [
        "chrom",
        "startposition",
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

    with tempfile.TemporaryDirectory() as temp_dir:
        file_basename = os.path.basename(file)
        file_temp_path_hapl = f"{temp_dir}/{file_basename}"
        file_temp_path_nohapl = f"{temp_dir}/{file_basename}/{file_basename}.bed"
        log_temp_path = f"{temp_dir}/{file_basename}/logging.log"

        if hapl:
            if mod == "m":
                try:
                    logging.info(
                        "Extract modification frequencies from haplotype 1 in bam/cram files using modkit tool"
                    )
                    modkit_stream = subprocess.Popen(
                        shlex.split(
                            f"modkit pilup {file} {file_temp_path_hapl} --ignore h --region={window.fmt} --cpg --ref={fasta} --log-filepath {log_temp_path} --partition-tag HP --prefix H --only-tabs"
                        ),
                        stderr=subprocess.PIPE,
                    )
                    logging.info(
                        "Extract modification frequencies haplotypes in bam/cram files using modkit tool"
                    )
                except FileNotFoundError as e:
                    logging.error(e, exc_info=True)
                    sys.stderr.write(
                        "\n\nError when making bedfile of bam/cram file with modkit.\n"
                    )
                    sys.stderr.write("Is modkit installed and on the PATH?")
                    sys.stderr.write(
                        f"\n\n\nDetailed error: {modkit_stream.stderr.read()}\n"
                    )
                    raise
                if modkit_stream.returncode:
                    sys.exit(f"Received modkit error:\n{modkit_stream.stderr.read()}\n")
            if mod == "h":
                try:
                    logging.info(
                        "Extract modification frequencies from haplotype 1 in bam/cram files using modkit tool"
                    )
                    modkit_stream = subprocess.Popen(
                        shlex.split(
                            f"modkit pilup {file} {file_temp_path_hapl} --ignore m --region={window.fmt} --cpg --ref={fasta} --log-filepath {log_temp_path} --partition-tag HP --prefix H --only-tabs"
                        ),
                        stderr=subprocess.PIPE,
                    )
                    logging.info(
                        "Extract modification frequencies haplotypes in bam/cram files using modkit tool"
                    )
                except FileNotFoundError as e:
                    logging.error(e, exc_info=True)
                    sys.stderr.write(
                        "\n\nError when making bedfile of bam/cram file with modkit.\n"
                    )
                    sys.stderr.write("Is modkit installed and on the PATH?")
                    sys.stderr.write(
                        f"\n\n\nDetailed error: {modkit_stream.stderr.read()}\n"
                    )
                    raise
                if modkit_stream.returncode:
                    sys.exit(f"Received modkit error:\n{modkit_stream.stderr.read()}\n")

            logging.info("Read the file in a dataframe per haplotype.")

            df_H1 = pd.read_table(
                f"{file_temp_path_hapl}/H_1.bed",
                sep="\t",
                header=None,
                names=headerlist,
                usecols=["startposition", "endposition", "fractionmodified"],
            )

            df_H2 = pd.read_table(
                f"{file_temp_path_hapl}/H_2.bed",
                sep="\t",
                header=None,
                names=headerlist,
                usecols=["startposition", "endposition", "fractionmodified"],
            )

            df_H1["position"] = (df_H1["startposition"] + df_H1["endposition"]) / 2
            df_H1 = df_H1.set_index("position").drop(
                ["startposition", "endposition"], axis=1, inplace=True
            )

            df_H2["position"] = (df_H2["startposition"] + df_H2["endposition"]) / 2
            df_H2 = df_H2.set_index("position").drop(
                ["startposition", "endposition"], axis=1, inplace=True
            )

            dfs_1.append(df_H1.rename(columns={"methylated_frequency": f"{name}_1"}))
            dfs_2.append(df_H2.rename(columns={"methylated_frequency": f"{name}_2"}))

            import itertools

            dfs = list(itertools.chain(*zip(dfs_1, dfs_2)))
            return dfs

        if not hapl:
            if mod == "m":
                try:
                    logging.info(
                        "Extract modification frequencies from bam/cram files using modkit tool"
                    )
                    modkit_stream = subprocess.Popen(
                        shlex.split(
                            f"modkit pilup {file} {file_temp_path_nohapl} --ignore h --region={window.fmt} --cpg --ref={fasta} --log-filepath {log_temp_path} --only-tabs"
                        ),
                        stderr=subprocess.PIPE,
                    )
                except FileNotFoundError as e:
                    logging.error(e, exc_info=True)
                    sys.stderr.write(
                        "\n\nError when making bedfile of bam/cram file with modkit.\n"
                    )
                    sys.stderr.write("Is modkit installed and on the PATH?")
                    sys.stderr.write(
                        f"\n\n\nDetailed error: {modkit_stream.stderr.read()}\n"
                    )
                    raise
                if modkit_stream.returncode:
                    sys.exit(f"Received modkit error:\n{modkit_stream.stderr.read()}\n")
            if mod == "h":
                try:
                    logging.info(
                        "Extract modification frequencies from bam/cram files using modkit tool"
                    )
                    modkit_stream = subprocess.Popen(
                        shlex.split(
                            f"modkit pilup {file} {file_temp_path_nohapl} --ignore m --region={window.fmt} --cpg --ref={fasta} --log-filepath {log_temp_path} --only-tabs"
                        ),
                        stderr=subprocess.PIPE,
                    )
                except FileNotFoundError as e:
                    logging.error(e, exc_info=True)
                    sys.stderr.write(
                        "\n\nError when making bedfile of bam/cram file with modkit.\n"
                    )
                    sys.stderr.write("Is modkit installed and on the PATH?")
                    sys.stderr.write(
                        f"\n\n\nDetailed error: {modkit_stream.stderr.read()}\n"
                    )
                    raise
                if modkit_stream.returncode:
                    sys.exit(f"Received modkit error:\n{modkit_stream.stderr.read()}\n")

                logging.info("Read the file in a dataframe.")
            df = pd.read_table(
                file_temp_path_nohapl,
                sep="\t",
                header=None,
                names=headerlist,
                usecols=["startposition", "endposition", "fractionmodified"],
            )
            df["position"] = (df["startposition"] + df["endposition"]) / 2
            df = df.set_index("position").drop(
                ["startposition", "endposition"], axis=1, inplace=True
            )

            dfs.append(df.rename(columns={"methylated_frequency": name}))
            return dfs


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
    else:
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
