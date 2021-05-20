"""Helper class containing the HG38 reference genome"""
# Disable no-member error to work around https://github.com/pysam-developers/pysam/issues/819
# pylint: disable=no-member
import os
import pysam


class ReferenceGenome:
    """Holds an index to the reference genome and associated lookup methods"""

    offset = 100

    def __init__(self, filename):
        """Load reference genome

        Positional arguments:
        filename -- genome in fasta format
        """
        try:
            if not os.path.exists("{}.fai".format(filename)):
                pysam.faidx(filename)
            self._genome = pysam.FastaFile(filename)
        except pysam.utils.SamtoolsError as samtools_error:
            print(samtools_error)
            exit(1)

    def get_peak(self, reference, peak_position, expected_base=None):
        """Fetch a region of the genome centered on ``peak_position``, surrounded by
        an offset of base pairs on either side.

        Note that ``peak_position`` is one-indexed
        """

        start, end = self.get_peak_coords(peak_position)
        sequence = self.get_region(reference, start, end)

        # If we expect the peak position to have a particular base, assert this.
        if expected_base:
            base_location = min(self.offset, peak_position - 1)
            found_base = sequence[base_location].upper()
            assert (
                found_base == expected_base
            ), "Reference genome doesn't match expectations at {}:{}, expected {} found {}".format(
                reference, peak_position, expected_base, found_base
            )

        return sequence

    @classmethod
    def get_peak_coords(cls, peak_position):
        """Get the start and end coordinates for a region centered on ``peak_position``,
        with ``offset`` base pairs on either side.

        Input:
        ``peak_position``: 1-indexed position of peak in genome

        Output:
        0-indexed tuple of start and end position.
        """
        start = max(peak_position - cls.offset - 1, 0)
        end = max(peak_position + cls.offset, 2 * cls.offset + 1)
        return (start, end)

    def get_region(self, reference, start, end):
        """Fetch a region of the genome

        Positional arguments:
        reference -- region of the genome to access (i.e. chromosome number)
        start -- start position within region
        end -- end position within region
        """

        return self._genome.fetch(reference=reference, start=start, end=end).upper()
