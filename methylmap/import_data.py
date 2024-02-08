import shlex
import sys
import os
import gzip
import logging
import subprocess
import numpy as np
import pandas as pd
from pathlib import Path
from methylmap.region import Region
from concurrent.futures import ThreadPoolExecutor
import itertools
import tempfile


def read_mods(files, table, names, window, groups, gff, fasta, mod, hapl, dendro):
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
            rc = subprocess.call(["which", "modkit"])
            if not rc == 0:
                sys.exit(
                    "\n\n\nIs modkit installed? Installation: see https://github.com/nanoporetech/modkit"
                )
            else:
                return parse_bam(
                    files, table, names, window, groups, fasta, mod, hapl, dendro
                )
        elif file_type in ["overviewtable_bam", "overviewtable_cram"]:
            rc = subprocess.call(["which", "modkit"])
            if not rc == 0:
                sys.exit(
                    "\n\n\nIs modkit installed? Installation: see https://github.com/nanoporetech/modkit"
                )
            else:
                return parse_bam(
                    files, table, names, window, groups, fasta, mod, hapl, dendro
                )
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
                        shlex.split(f"tabix -S1 -s1 -b2 -e3 {file}"),
                        stderr=subprocess.PIPE,
                    )
                except FileNotFoundError as e:
                    logging.error("Error when making a .tbi file.")
                    logging.error(e, exc_info=True)
                    sys.stderr.write("\n\nERROR when making a .tbi file.\n")
                    sys.stderr.write("Is tabix installed and on the PATH?.")
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
                sys.stderr.write("Is tabix installed and on the PATH?")
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
        table = table.set_index("position").drop(["chromosome", "start", "end"], axis=1)
        dfs.append(table.rename(columns={"methylated_frequency": name}))
    methfrequencytable = dfs[0].join(dfs[1:], how="outer")

    if len(methfrequencytable) == 0:
        logging.error(
            "WARNING: length of methylation frequency table is zero. Do the input files contain data?"
        )
        sys.exit(
            "WARNING: length of methylation frequency table is zero. Do the input files contain data?"
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
            "WARNING: length of methylation frequency table is zero. Does the input table contain data?"
        )
        sys.exit(
            "WARNING: length of methylation frequency table is zero. Does the input table contain data?"
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
        file_temp_path_hapl = f"{temp_dir}"
        file_temp_path_nohapl = f"{temp_dir}/{file_basename}/{file_basename}.bed"
        log_temp_path = f"{temp_dir}/{file_basename}/logging.log"
        print("Temporary directory created at:", temp_dir)
        assert os.path.isdir(temp_dir)  
        print(os.listdir(file_temp_path_hapl))

        if hapl:
            print("Temporary directory created at:", temp_dir)
            ignore = "h" if mod == "m" else "m"
            try:
                logging.info(
                    f"Processing file {file}: Extract modification frequencies per haplotype in bam/cram files using modkit tool"
                )
                modkit_cmd = f"modkit pileup {file} {temp_dir} --ignore {ignore} --region={window.fmt} --log-filepath {log_temp_path} --partition-tag HP --prefix H --only-tabs"
                print(os.listdir(file_temp_path_hapl))
                print(modkit_cmd)
                modkit_stream = subprocess.Popen(
                    shlex.split(modkit_cmd),
                    stderr=subprocess.PIPE
                )
                logging.info(
                    f"Extract modification frequencies haplotypes in bam/cram files using modkit tool for file {file}."
                )
                print("test1")  
                assert os.path.isdir(file_temp_path_hapl)
                print(os.listdir(file_temp_path_hapl))
                print("test2")

            except FileNotFoundError as e:
                logging.error(e, exc_info=True)
                sys.stderr.write(
                    f"\n\nError when making bedfile of bam/cram file with modkit for file {file}.\n"
                )
                sys.stderr.write(
                    "Is modkit installed and on the PATH? Is temporary directory accessible?"
                )
                sys.stderr.write(
                    f"\n\n\nDetailed error: {modkit_stream.stderr.read()}\n"
                )
                raise
            if modkit_stream.returncode:
                sys.exit(f"Received modkit error:\n{modkit_stream.stderr.read()}\n")

            logging.info(
                f"Read the file in a dataframe per haplotype for sample {file}."
            )

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
            ignore = "h" if mod == "m" else "m"
            try:
                logging.info(
                    "Extract modification frequencies from bam/cram files using modkit tool"
                )
                modkit_stream = subprocess.Popen(
                    shlex.split(
                        f"modkit pilup {file} {file_temp_path_nohapl} --ignore {ignore} --region={window.fmt} --ref={fasta} --log-filepath {log_temp_path} --only-tabs"
                    ),
                    stderr=subprocess.PIPE,
                )
            except FileNotFoundError as e:
                logging.error(e, exc_info=True)
                sys.stderr.write(
                    "\n\nError when making bedfile of bam/cram file with modkit.\n"
                )
                sys.stderr.write(
                    "Is modkit installed and on the PATH? IS temporary directory accessible?"
                )
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


def parse_bam(files, table, names, window, groups, fasta, mod, hapl, dendro):
    if not fasta:
        logging.info("Stop script when no --fasta input")
        sys.exit(
            "ERROR when parsing bam/cram file, can not find fasta file. Is fasta file given with --fasta argument?"
        )
    if table:
        logging.info("Extract files and names from overviewtable")
        files, names = parse_overviewtable(table)

    with ThreadPoolExecutor(max_workers=1) as executor:
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

    if files:
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
