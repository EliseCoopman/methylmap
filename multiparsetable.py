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
        "--output", help="TSV file to write the frequencies to [default: stdout]"
    )
    parser.add_argument(
        "--groups", nargs="*", help="list of experimental group for each sample"
    )
    parser.add_argument(
        "--fasta",
        help="fasta reference file, required when input is BAM/CRAM files or overviewtable with BAM/CRAM files",
        required=True,
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
    return args


def main():
    args = get_args()
    window = Region(args.window, args.expand)
    file_sniffer(args.files[0])
    check_modkit()
    try:
        overviewtable = parse_bam(args, window)
        if args.output:
            overviewtable.to_csv(args.output, sep="\t", na_rep=np.NaN, header=True)
        else:
            print(overviewtable.to_csv(sep="\t", na_rep=np.NaN, header=True))
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

    methfreqtable = methfrequencytable.sort_values("position", ascending=True)

    if args.groups:
        logging.info("Sort columns of methfrequencytable based on group")
        if args.hapl:
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


def process_single_file(function_args):
    input, name, args, window = function_args

    with tempfile.TemporaryDirectory() as temp_dir:
        file_basename = os.path.basename(input)
        logging.info("Extract mod frequencies in bam/cram using modkit")
        # when processing haplotypes, <output> is a directory
        output = (
            temp_dir if args.hapl else os.path.join(temp_dir, f"{file_basename}.bed")
        )

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
        logging.info("Read the file in a dataframe per haplotype.")
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

    logging.info("Read the modkit file in a dataframe.")
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


def file_sniffer(filename):
    """
    Takes in a filename and tries to guess the input file type.
    """
    if not Path(filename).is_file():
        sys.exit(f"\n\nERROR: File {filename} does not exist, please check the path!\n")
    if not (is_bam_file(filename) or is_cram_file(filename)):
        sys.exit(f"\n\n\nInput file type {filename} not recognized!\n")


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
