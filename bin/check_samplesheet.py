#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FQ_FORMATS = (
        ".fq.gz",
        ".fastq.gz",
        ".fq",
        ".fastq",
    )

    VALID_FA_FORMATS = (
        ".fa.gz",
        ".fasta.gz",
        ".fa",
        ".fasta",
    )

    def __init__(
        self,
        sample_col="sample",
        first_col="fastq",
        second_col="fasta",
        third_col="contig",
        fourth_col="target_size",
        fifth_col="replicon",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the FASTQ file path
                (default "fastq").
            second_col (str): The name of the column that contains the FASTA file
                (default "fasta").
            third_col (str): The name of the column that contains the contig name
                (default "contig").
            fourth_col (str): The name of the column that contains the contig's target
                size (default "target_size").
            fifth_col (str): The name of the column that contains the replicon name
                (default "replicon").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._second_col = second_col
        self._third_col = third_col
        self._fourth_col = fourth_col
        self._fifth_col = fifth_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_first(row)
        self._validate_second(row)
        self._validate_third(row)
        self._validate_fourth(row)
        self._validate_fifth(row)
        self._seen.add((row[self._sample_col], row[self._third_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the FASTQ entry is non-empty and has the right format."""
        if len(row[self._first_col]) <= 0:
            raise AssertionError("The FASTQ file is required.")
        self._validate_fastq_format(row[self._first_col])

    def _validate_second(self, row):
        """Assert that the FASTA entry is non-empty and has the right format."""
        if len(row[self._second_col]) <= 0:
            raise AssertionError("The FASTA file is required.")
        self._validate_fasta_format(row[self._second_col])

    def _validate_third(self, row):
        """Assert that the contig name is non-empty."""
        if len(row[self._third_col]) <= 0:
            raise AssertionError("The contig name is required.")

    def _validate_fourth(self, row):
        """Assert that the target size is non-empty and is an integer greater than 0."""
        if len(row[self._fourth_col]) <= 0:
            raise AssertionError("The target size is required.")
        elif not row[self._fourth_col].isdigit():
            raise AssertionError("The target size should be an integer greater than 0")
        elif int(row[self._fourth_col]) <= 0:
            raise AssertionError("The target size should be an integer greater than 0")
            
    def _validate_fifth(self, row):
        """Assert that the replicon name is non-empty."""
        if len(row[self._fifth_col]) <= 0:
            raise AssertionError("The replicon name is required.")

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has one of the expected FASTQ extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_FQ_FORMATS):
            raise AssertionError(
                f"The FASTQ file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_FQ_FORMATS)}"
            )

    def _validate_fasta_format(self, filename):
        """Assert that a given filename has one of the expected FASTA extensions."""
        if not any(filename.endswith(extension) for extension in self.VALID_FA_FORMATS):
            raise AssertionError(
                f"The FASTA file has an unrecognized extension: {filename}\n"
                f"It should be one of: {', '.join(self.VALID_FA_FORMATS)}"
            )

    def validate_unique_samples(self):
        """Assert that the combination of sample name and contig name is unique."""
        if len(self._seen) != len(self.modified):
            raise AssertionError("The pair of sample name and contig name must be unique.")
        for row in self.modified:
            sample = row[self._sample_col]
            contig = row[self._third_col]
            row[self._sample_col] = f"{sample}_{contig}"


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by nf-core pipelines.

    Validate the general shape of the table, expected columns, and each row.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure,

            sample,fastq,fasta,contig,target_size,replicon
            SAMPLE1,SAMPLE1.fastq.gz,stylo/SAMPLE1/flye/00-assembly/draft_assembly.fasta,contig_2,4448,ColpHAD28
            SAMPLE1,SAMPLE1.fastq.gz,stylo/SAMPLE1/flye/00-assembly/draft_assembly.fasta,contig_3,33063,IncX1
            SAMPLE1,SAMPLE1.fastq.gz,stylo/SAMPLE1/flye/00-assembly/draft_assembly.fasta,contig_4,4218,Col8282
            SAMPLE2,SAMPLE2.fastq.gz,stylo/SAMPLE2/flye/00-assembly/draft_assembly.fasta,contig_3,4448,ColpHAD28
            SAMPLE2,SAMPLE2.fastq.gz,stylo/SAMPLE2/flye/00-assembly/draft_assembly.fasta,contig_5,2096,ColpVC

    """
    required_columns = {"sample", "fastq", "fasta", "contig", "target_size", "replicon"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
        checker.validate_unique_samples()
    header = list(reader.fieldnames)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
