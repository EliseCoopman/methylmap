import sys
import gzip
import shlex
import logging
import subprocess
import pandas as pd
from pathlib import Path


def process_1000Genomes(file, window):
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
        sys.stderr.write(f"\n\n\nDetailed error: {tabix_stream.stderr.read()}\n")
        raise
    if tabix_stream.returncode:
        sys.exit(f"\n\n\nReceived tabix error\n\n{tabix_stream.stderr.read()}")

    header = gzip.open(file, "rt").readline().rstrip().split("\t")
    modfreqtable_1000Genomes = pd.read_csv(
        tabix_stream.stdout, sep="\t", header=None, names=header
    ).set_index(["chrom", "position"])
    return modfreqtable_1000Genomes
