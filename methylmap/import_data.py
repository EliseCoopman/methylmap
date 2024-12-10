from methylmap.region import Region

import os
import io
import sys
import gzip
import shlex
import base64
import logging
import tempfile
import itertools
import subprocess
import pandas as pd
from tqdm import tqdm
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor


def read_mods(args, window, upload_data, filename, last_modified):
    """
    Deciding of input file(s) type and processing them.
    """
    if args == False:
        if upload_data:
            file_type = "modfrequencytable"
            return parse_modfrequencytable(
                args, window, upload_data, filename, last_modified
            )
    if not args.files and not args.table:
        return
    elif args.files:
        file_type = file_sniffer(args.files[0])
    elif args.table:
        file_type = file_sniffer(args.table)
    try:
        if file_type == "nanopolish_calc_meth_freq":
            return parse_nanopolish(args, window)
        elif file_type == "overviewtable_nanopolishfiles":
            return parse_nanopolish(args, window)
        elif file_type == "modfrequencytable":
            return parse_modfrequencytable(args, window)
        elif file_type in ["cram", "bam"]:
            check_modkit()
            rc = subprocess.call(["which", "modkit"])
            return parse_bam(args, window)
        elif file_type in ["overviewtable_bam", "overviewtable_cram"]:
            check_modkit()
            return parse_bam(args, window)
    except Exception as e:
        logging.error("Error processing input file(s).")
        logging.error(e, exc_info=True)
        sys.stderr.write("\n\n\nError processing input file(s)!\n")
        raise


def check_modkit():
    rc = subprocess.call(["which", "modkit"], stdout=subprocess.DEVNULL)
    if rc != 0:
        sys.exit("\n\nIs modkit installed? See github.com/nanoporetech/modkit")


def parse_overviewtable(table):
    """
    Parses overviewtable, sorts files by group and makes list of all files and sample names.
    """
    overviewtable = pd.read_table(table).sort_values(["group", "name"])
    files = overviewtable["path"].tolist()
    names = overviewtable["name"].tolist()
    return files, names


def parse_nanopolish(args, window):
    """
    Converts a file from nanopolish to a pandas dataframe
    input can be from calculate_methylation_frequency.py or overviewtable with these files
    which will return a dataframe with 'chromosome', 'pos', 'methylated_frequency'.
    """
    if args.table:
        logging.info("Extract files and names from overviewtable")
        files, names = parse_overviewtable(args.table)
    else:
        files, names = args.files, args.names

    dfs = []
    for file, name in zip(files, names):
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
        table["position"] = (table["start"] + table["end"]) / 2
        table = table.set_index(["chromosome", "position"]).drop(
            ["start", "end"], axis=1
        )
        dfs.append(table.rename(columns={"methylated_frequency": name}))
    modfrequencytable = dfs[0].join(dfs[1:], how="outer")

    modfrequencytable.reset_index(inplace=True)
    modfrequencytable.rename(columns={"chromosome": "chrom"}, inplace=True)
    if len(modfrequencytable) == 0:
        err = "WARNING: length of modification frequency table is zero. Do the input files contain data?"
        logging.error(err)
        sys.exit(err)

    modfreqtable = modfrequencytable.sort_values("position", ascending=True)
    # output is a modification frequency table with position as index and for each sample a column with all modification frequencies
    modfreqtable.set_index(["chrom", "position"], inplace=True)
    if args.files:
        if args.groups:
            if args.dendro:
                err = "Columns will not be sorted based on --group input since hierarchical clustering with --dendro is requested."
                logging.warning(err)
                sys.stderr.write(err)
            else:
                logging.info("Sort columns of modfrequencytable based on group")
                headerlist = list(modfreqtable.columns.values)
                if len(headerlist) == len(args.groups):
                    res = zip(headerlist, args.groups)
                    output = sorted(list(res), key=lambda x: x[1])
                    orderedlist = [i[0] for i in output]
                    modfreqtable = modfreqtable.reindex(columns=orderedlist)
                else:
                    sys.exit(
                        f"ERROR when matching --groups with samples, is length of --groups list ({len(args.groups)}) matching with number of sample files?"
                    )

    return modfreqtable


def parse_modfrequencytable(
    args, window, upload_data=False, filename=None, last_modified=None
):
    """
    Parsing modfrequencytable input.
    """
    if args is False:
        content_type, content_string = upload_data.split(",")
        decoded = base64.b64decode(content_string)
        try:
            if "tsv" in filename:
                if window is not None:
                    df = pd.read_csv(
                        io.StringIO(decoded.decode("utf-8")), sep="\t"
                    ).sort_values("position", ascending=True)
                    df = df[df["chrom"] == window.chromosome]
                    df = df[df["position"].between(window.begin, window.end)]
                    df.drop(["chrom"], axis=1, inplace=True)
                    df.set_index(["position"], inplace=True)
                if window is None:
                    df = pd.read_csv(
                        io.StringIO(decoded.decode("utf-8")), sep="\t", nrows=500
                    ).sort_values("position", ascending=True)
                    first_chromosome = df["chrom"].unique()[0]
                    df = df[df["chrom"] == first_chromosome]
            else:
                return
        except Exception as e:
            print(e)
            return

    else:
        df = pd.read_table(args.table).sort_values("position", ascending=True)
        if window:
            df = df[df["chrom"] == window.chromosome]
            logging.info("Select window out of modfreqtable")
            df = df[df["position"].between(window.begin, window.end)]
        else:
            if args.gff:
                logging.info(
                    "If no window given and annotation requested, take window out of modfreqtable."
                )
                if len(df["chrom"].unique()) > 1:
                    sys.exit(
                        "\n\nError when extracting window out of modfreqtable positions. Chrom column can not contain more than one chromosome \n"
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

        df.set_index(["chrom", "position"], inplace=True)
        if len(df) == 0:
            err = "WARNING: length of modification frequency table is zero. Does the input table contain data?"
            logging.error(err)
            sys.exit(err)

        if args.names:
            logging.info("Changes sample names to names from --names input")
            df.columns = args.names

        if args.groups:
            logging.info("Sort columns of modfrequencytable based on group")
            if args.dendro:
                err = "Columns will not be sorted based on --group input since hierarchical clustering with --dendro is requested."
                logging.warning(err)
                sys.stderr.write(err)
            else:
                headerlist = list(df.columns.values)
                if len(headerlist) == len(args.groups):
                    res = zip(headerlist, args.groups)
                    output = sorted(list(res), key=lambda x: x[-1])
                    orderedlist = [i[0] for i in output]
                    df = df.reindex(columns=orderedlist)
                else:
                    err = f"Error when matching --groups with samples. Is length of the --groups list ({len(args.groups)}) matching with number of samples in table?"
                    logging.error(err)
                    sys.exit(err)
    return df


def process_single_file(function_args):
    input, name, window, args = function_args
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
            log_file=os.path.join("/home/ecoopman/modkit.log"),
            haplotype=args.hapl,
        )
        if not args.quiet:
            tqdm.write(stderr, file=sys.stderr)
        if args.hapl:
            logging.info("Read the file in a dataframe per haplotype.")
            return [
                process_modkit_tsv(
                    f"{output}/H_{haplotype}.bed", name=f"{name}_{haplotype}"
                )
                for haplotype in [1, 2]
            ]
        else:
            return [process_modkit_tsv(output, name=name)]


def run_modkit_pileup(input, output, ignore, window, fasta, log_file, haplotype=False):
    cmd = f"modkit pileup {input} {output} --ignore {ignore} --region={window.fmt} --ref={fasta} --log-filepath {log_file} --only-tabs --suppress-progress"
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


def parse_bam(args, window):
    if not args.fasta:
        logging.info("Error: No --fasta input, this is required.")
    if args.files:
        args_list = [
            (file, name, window, args) for file, name in zip(args.files, args.names)
        ]
    if args.table:
        logging.info("Extract files and names from overviewtable")
        files, names = parse_overviewtable(args.table)
        args_list = [
            (file, name, window, args) for file, name in zip(files, names)
        ]

    # Process files concurrently with a maximum of 12 threads
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        results = list(
            tqdm(executor.map(process_single_file, args_list), total=len(args_list))
        )

    dfs = list(itertools.chain(*results))

    modfrequencytable = dfs[0].join(dfs[1:], how="outer")

    if len(modfrequencytable) == 0:
        err = "WARNING: length of modification frequency table is zero. Do the input files contain data?"
        logging.error(err)
        sys.exit(err)

    modfreqtable = modfrequencytable.sort_values("position", ascending=True)

    if args.files:
        if args.groups:
            logging.info("Sort columns of modfrequencytable based on group")
            if args.dendro:
                err = "Columns will not be sorted based on --group input since hierarchical clustering with --dendro is requested."
                logging.warning(err)
                sys.stderr.write(err)
            else:
                if args.hapl:
                    groupshapl = list(itertools.chain(*zip(groups, groups)))
                    groups = groupshapl
                else:
                    groups = args.groups
                headerlist = list(modfreqtable.columns.values)
                if len(headerlist) == len(groups):
                    res = zip(headerlist, groups)
                    output = sorted(list(res), key=lambda x: x[1])
                    orderedlist = [i[0] for i in output]
                    modfreqtable = modfreqtable.reindex(columns=orderedlist)
                else:
                    sys.exit(
                        f"ERROR when matching --groups with samples, is length of --groups list ({len(groups)}) matching with number of sample files?"
                    )

    return modfreqtable


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
    if "path" in header:  # overviewtable
        df = pd.read_table(filename)
        if (
            df["path"].iloc[0].endswith(".tsv")
        ):  ###only checks first file in overviewtable
            return "overviewtable_nanopolishfiles"
        elif (
            df["path"].iloc[0].endswith(".tsv.gz")
        ):  ###only checks first file in overviewtable
            return "overviewtable_nanopolishfiles"
        elif (
            df["path"].iloc[0].endswith(".bam")
        ):  ####only checks first file in overviewtable
            return "overviewtable_bam"
        elif (
            df["path"].iloc[0].endswith(".cram")
        ):  ####only checks first file in overviewtable
            return "overviewtable_cram"
    if header.startswith("chrom"):
        # own modfrequencytable as input: needs first column header to be "chrom"
        return "modfrequencytable"
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
