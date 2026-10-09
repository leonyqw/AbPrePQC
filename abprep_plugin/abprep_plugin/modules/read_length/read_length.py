import logging
from typing import Counter, Dict, List
import plotly.graph_objects as go
import numpy as np

# from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    MultiQC with read length outputs
    """

    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name="Read length",
            anchor="read_length",
            info="Comparison of read lengths across samples. Read lengths are obtained from the samtools stats output.",
        )

        # Find and load any read_length reports
        read_len_data: Dict[str, List[int]] = dict()

        for f in self.find_log_files("read_length"):
            s_name = f["s_name"]

            read_len_data[s_name] = self.parse_read_length(f["f"])

            if s_name in read_len_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f)

        # Report if no samples found
        if len(read_len_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(read_len_data)} reports")

        # Remove empty samples
        read_len_data = {s_name: vals for s_name, vals in read_len_data.items() if vals}

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        self.read_length_plots(read_len_data)

    def parse_read_length(self, f) -> List[int]:
        """Parse read length files"""

        parsed_data = [int(x) for x in f.split("\n") if x.strip()]

        return parsed_data

    def read_length_plots(self, read_len_data):
        """
        Read length distribution per sample.
        Weighted histogram plot adapted from NanoComp.
        """
        # Drop zero-length reads and samples left empty
        lengths = {s_name: np.asarray(vals, dtype=np.int64) for s_name, vals in read_len_data.items()}
        lengths = {s_name: vals[vals > 0] for s_name, vals in lengths.items()}

        # Generate bins
        all_lengths = np.concatenate(list(lengths.values()))
        bins = min(max(round(all_lengths.max() / 500), 10), 500)
        edges = np.histogram_bin_edges(all_lengths, bins=bins)
        midpoints = [round(float(m), 1) for m in (edges[:-1] + edges[1:]) / 2]  # bin midpoints

        # Calculate read counts and bases per sample per bin
        reads_by_sample, bases_by_sample = {}, {}
        for s_name, vals in lengths.items():
            counts, _ = np.histogram(vals, bins=edges)
            bases, _ = np.histogram(vals, bins=edges, weights=vals)
            reads_by_sample[s_name] = dict(zip(midpoints, counts.astype(int).tolist()))
            bases_by_sample[s_name] = dict(zip(midpoints, bases.astype(int).tolist()))

        # Configuration for the line plots
        pconfig = {
            "id": "read_length_plot",
            "title": "Read length distribution",
            "xlab": "Read length (bp)",
            "ymin": 0,
            "xmin": 0,
            "data_labels": [
                {"name": "Bases", "ylab": "Number of bases", "tt_label": "<b>~%{x:,.0f} bp</b>: %{y:,.0f} bases"},
                {"name": "Reads", "ylab": "Number of reads", "tt_label": "<b>~%{x:,.0f} bp</b>: %{y:,.0f} reads"},
            ],
        }

        # Add plot to the report
        self.add_section(
            name="Read length distribution",
            anchor="read_length_plot",
            description="Distribution of the number of reads / bases per sample.",
            comment="A good sample should have one big peak around the expected full length of the plasmid.",
            helptext="""
            Each line is one sample. Reads are placed into evenly spaced and equal sized bins.
            Use the buttons above the plot to switch views.

            **Bases**: Every read adds the number of bases to its read length bin. The y-axis is the total number of bases that fall in each read length bin.

            **Reads**: The y-axis is the number of reads in each read length bin.
            """,
            plot=linegraph.plot([bases_by_sample, reads_by_sample], pconfig),
        )
