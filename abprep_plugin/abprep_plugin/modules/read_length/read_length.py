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
            info="Comparison of read lengths across samples.",
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

        # self.read_length_violin_plot(read_len_data)

        self.read_length_plots(read_len_data)

    def parse_read_length(self, f) -> List[int]:
        """Parse read length files"""

        parsed_data = [int(x) for x in f.split("\n") if x.strip()]

        return parsed_data

    def read_length_violin_plot(self, read_len_data):
        """Generate an interactive violin plot of read lengths per sample"""

        # log10-transform (drop zero-length reads: log10(0) is undefined)
        lengths = {
            s_name: np.log10(v[v > 0])
            for s_name, v in ((s_name, np.asarray(vals, dtype=np.int64)) for s_name, vals in read_len_data.items())
        }

        # Create plot
        fig = go.Figure()
        for s_name, vals in lengths.items():
            fig.add_trace(
                go.Violin(
                    x=[s_name] * len(vals),
                    y=vals,
                    name=s_name,
                    points=False,
                    spanmode="hard",
                )
            )

        # Ticks at powers of 10, labelled as real read lengths
        lo = int(np.floor(min(v.min() for v in lengths.values())))
        hi = int(np.ceil(max(v.max() for v in lengths.values())))
        ticks = list(range(lo, hi + 1))

        # Update plot layout
        fig.update_layout(
            showlegend=False,
            autosize=True,
            height=520,
            yaxis=dict(
                title="Log-transformed read length",
                tickvals=ticks,
                ticktext=[f"{10**t:,}" for t in ticks],
            ),
            margin=dict(l=60, r=20, t=30, b=60),
        )

        # Only the download button in the modebar
        config = {
            "displaylogo": False,
            "responsive": True,
            "toImageButtonOptions": {
                "format": "svg",  # or "png"
                "filename": "read_length_violin",
                "scale": 2,
            },
        }

        # Add plot section to the report
        self.add_section(
            name="Read length distribution by barcode",
            anchor="barcode_readlen_violin",
            description="Violin plot of log-transformed read lengths per barcode.",
            content=fig.to_html(
                full_html=False,
                include_plotlyjs=False,  # MultiQC report already loads plotly.js
                div_id="barcode_readlen_violin_plot",
                config=config,
            ),
        )

    def read_length_plots(self, read_len_data):
        """
        Read length distribution per sample, with a toggle between read counts and bases per length bin.
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
            "showlegend": False,
            "style": "lines+markers",
            "x_bands": [
                {"from": 4000, "to": 6000, "color": "#009500", "opacity": 0.13},
                {"from": 7000, "to": 9000, "color": "#a07300", "opacity": 0.13},
            ],
            "data_labels": [
                {"name": "Bases", "ylab": "Number of bases", "tt_label": "<b>~%{x:,.0f} bp</b>: %{y:,.0f} bases"},
                {"name": "Reads", "ylab": "Number of reads", "tt_label": "<b>~%{x:,.0f} bp</b>: %{y:,.0f} reads"},
            ],
        }

        # Add plot to the report
        self.add_section(
            name="Read length distribution",
            anchor="read_length_plot",
            description="Number of reads / bases per sample.",
            helptext="""
            Each line is one sample. Reads are placed into equal size bins.
            Use the buttons above the plot to switch views.

            **Bases**: Every read adds its own length to its bin instead of a count
            of 1, so the y-axis shows how many sequenced bases fall in each bin.
            This shows where most of the data is, which the read count view hides
            because short reads dominate by number.

            **Reads**: The y-axis is the number of reads in each length bin.
            """,
            plot=linegraph.plot([bases_by_sample, reads_by_sample], pconfig),
        )
