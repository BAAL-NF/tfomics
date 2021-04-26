"""Helper class for processing output from AlleleSeq"""
from enum import Enum, unique

import numpy
import pandas

from tfomics.data_structures.reference_genome import ReferenceGenome


class AlleleSeqData:
    """Helper object holding ChIP-seq data for a cell line and associated FDR estimates"""

    _field_separator = r"\s*\t"

    @unique
    class CountColumns(Enum):
        """Columns required in the AlleleSeq data files"""

        P_VALUE = "SymPval"
        REFERENCE_SNP = "ref"
        WINNING_ALLELE = "winning"
        MATERNAL_ALLELE = "mat_all"
        PATERNAL_ALLELE = "pat_all"

    @unique
    class FdrColumns(Enum):
        """Columns required in FDR data files"""

        P_VALUE = "pval"
        FDR = "FDR"

    def __init__(self, name, count_file, fdr_file):
        self.candidates = None
        self.name = name

        try:
            self.count = pandas.read_csv(
                count_file,
                AlleleSeqData._field_separator,
                dtype={"snppos": numpy.int64},
                engine="python",
                index_col=["chrm", "snppos"],
            )

            # The first line is a header comment, the last four lines have summary data.
            # We skip both of these and read in the summary data separately.
            self.fdr = pandas.read_csv(
                fdr_file,
                AlleleSeqData._field_separator,
                skipfooter=4,
                skiprows=1,
                engine="python",
            )

            with open(fdr_file) as file:
                # The first 30 lines have been read in above, so ignore them.
                for _ in range(30):
                    file.readline()

                for line in file.readlines():
                    parsed_line = [i for i in filter(None, line.strip().split(" "))]
                    if parsed_line[0] == "target":
                        self.target = parsed_line[1]
                    elif parsed_line[0] == "before":
                        self.before_fdrs = parsed_line[1:]
                    elif parsed_line[0] == "after":
                        self.after_fdrs = parsed_line[1:]

            count_columns = self.count.columns.values.tolist()
            for column in AlleleSeqData.CountColumns:
                assert (
                    column.value in count_columns
                ), "AlleleSeq data missing {} column".format(column.value)

            fdr_columns = self.fdr.columns.values.tolist()
            for column in AlleleSeqData.FdrColumns:
                assert column.value in fdr_columns, "FDR data missing {} column".format(
                    column.value
                )

        except FileNotFoundError as not_found:
            print(not_found)
            exit(1)
        except AssertionError as missing_column:
            print("{}, did you import the right file?".format(missing_column))
            exit(1)

    def get_pval(self, fdr):
        """Find P-value to use when selecting entries to ensure a
        False Discovery Rate <= fdr
        """
        candidates = self.fdr.query(
            "{} <= {}".format(AlleleSeqData.FdrColumns.FDR.value, fdr)
        )["pval"]

        if candidates.empty:
            print(
                "Warning: sample {} has no data points with FDR <= {}".format(
                    self.name, fdr
                )
            )
            return 0

        # Return the highest p value found with FDR < fdr
        return candidates.max()

    def get_candidates(self, fdr):
        """Select all SNPs with a small enough P value to achieve an FDR <= fdr
        and store them in self.candidates
        """

        # pval is used in the query below, but pylint doesn't understand the @-reference.
        pval = self.get_pval(fdr)  # pylint: disable=unused-variable

        # Occasionally a row is marked as weird because there are
        # reads that don't match the paternal or maternal allele.
        candidates = self.count.query(
            '({} <= @pval) and (winning != "?")'.format(
                AlleleSeqData.CountColumns.P_VALUE.value
            )
        )

        return candidates

    def get_row(self, chromosome, location):
        """Fetch a SNP at ``chromsome``:``loc``"""
        try:
            return self.count.loc[chromosome, location]
        except KeyError:
            return None
        except TypeError:
            return None

    def get_het_snp(self, chromosome, location, reference_snp=None):
        """Look up a SNP and return the nucleotide of the preferentially bound allele"""

        snp = self.get_row(chromosome, location)

        if reference_snp:
            expected_reference = snp[AlleleSeqData.CountColumns.REFERENCE_SNP.value]
            assert (
                expected_reference == reference_snp
            ), "Error: reference SNP mismatch at {}:{}, expected {}, found {}".format(
                chromosome, location, expected_reference, reference_snp
            )
        return AlleleSeqData.__get_winning_snp(snp)

    def get_winning_allele(self, chromosome, locations, pick_min=False):
        """Pool the results for all locations and return the overall preferentially bound allele"""
        rows = list(
            map(lambda location: self.get_row(chromosome, location.POS), locations)
        )

        if any(elem is None for elem in rows):
            return None

        maternal = 0
        paternal = 0
        for row in rows:
            maternal_allele = "c{}".format(
                row[AlleleSeqData.CountColumns.MATERNAL_ALLELE.value]
            )
            paternal_allele = "c{}".format(
                row[AlleleSeqData.CountColumns.PATERNAL_ALLELE.value]
            )
            maternal += row[maternal_allele]
            paternal += row[paternal_allele]

        row_dataframe = pandas.DataFrame(rows)

        if (maternal > paternal) ^ pick_min:
            return row_dataframe[
                [
                    AlleleSeqData.CountColumns.MATERNAL_ALLELE.value,
                    AlleleSeqData.CountColumns.REFERENCE_SNP.value,
                ]
            ]
        else:
            return row_dataframe[
                [
                    AlleleSeqData.CountColumns.PATERNAL_ALLELE.value,
                    AlleleSeqData.CountColumns.REFERENCE_SNP.value,
                ]
            ]

    @classmethod
    def set_reference_genome(cls, genome):
        """Assign a reference genome to the class, for use in creating sequences"""
        cls.genome = genome

    @classmethod
    def __create_sequences(cls, candidates):
        """Fetch sequences from the reference genome, modify them with the SNP as appropriate,
        and return a dataframe with the sequences inserted
        """
        sequences = []
        for (chromosome, position), row in candidates.iterrows():
            # Fetch the sequence from the reference genome
            sequence = cls.__get_sequence_from_reference_genome(
                chromosome, position, row
            )
            winning_snp = cls.__get_winning_snp(row)

            if winning_snp:
                # Mutate the sequence. Python strings are immutable so we create a copy
                sequence = (
                    sequence[: ReferenceGenome.offset]
                    + winning_snp
                    + sequence[ReferenceGenome.offset + 1 :]
                )
                sequences.append(sequence)

        return candidates.assign(sequence=sequences)

    @classmethod
    def __get_sequence_from_reference_genome(cls, chromosome, position, row):
        """Fetch a sequence from the reference genome and confirm
        it matches the expected value"""

        reference_base = row[cls.CountColumns.REFERENCE_SNP.value]

        # Fetch the sequence from the reference genome
        sequence = cls.genome.get_peak(
            chromosome, position, expected_base=reference_base
        )

        return sequence

    @staticmethod
    def __get_winning_snp(row):
        """Determine whether the TF bound more to the paternal or maternal allele and
        return the nucleotide that the TF bound to preferentially"""
        if row[AlleleSeqData.CountColumns.WINNING_ALLELE.value] == "P":
            return row[AlleleSeqData.CountColumns.PATERNAL_ALLELE.value]

        if row[AlleleSeqData.CountColumns.WINNING_ALLELE.value] == "M":
            return row[AlleleSeqData.CountColumns.MATERNAL_ALLELE.value]

        return None
