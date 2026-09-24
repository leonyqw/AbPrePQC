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

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        self.read_length_violin_plot(read_len_data)

        self.read_length_line_plot(read_len_data)

    def parse_read_length(self, f) -> List[int]:
        """Parse read length files"""

        parsed_data = [int(x) for x in f.split("\n") if x.strip()]

        return parsed_data

    def read_length_line_plot(self, read_len_data):
        """Generate line graph plot of read lengths per sample"""

        # {sample: {length: count}}
        data_by_sample = {s_name: dict(sorted(Counter(lengths).items())) for s_name, lengths in read_len_data.items()}

        ####### FASTQC line plot example

        # data_by_sample: Dict[str, Dict[int, int]] = dict()
        # for s_name, sd in self.fastqc_data.items():
        #     if sd.get("per_sequence_quality_scores") is None:
        #         continue
        #     data_by_sample[s_name] = {d["quality"]: d["count"] for d in sd["per_sequence_quality_scores"]}

        # Convert status dict format
        # status_dict: Dict[Literal["pass", "warn", "fail"], List[str]] = {"pass": [], "warn": [], "fail": []}
        # for s_name, status in section_statuses.items():
        #     if status in status_dict:
        #         status_dict[status].append(s_name)

        pconfig = {
            # "id": f"{self.anchor}_per_sequence_quality_scores_plot",
            "id": "read_length_line_plot",
            "title": "Per sample read length distribution",
            "ylab": "Number of reads",
            "xlab": "Read length (bp)",
            "ymin": 0,
            "xmin": 0,
            # "x_decimals": False,
            # "tt_label": "<b>Phred {point.x}</b>: {point.y} reads",
            "showlegend": False,
            # "colors": self.get_status_cols("per_sequence_quality_scores"),
            "x_bands": [
                {"from": 6000, "to": 10000, "color": "#009500", "opacity": 0.13},
                # {"from": 20, "to": 28, "color": "#a07300", "opacity": 0.13},
                # {"from": 0, "to": 20, "color": "#990101", "opacity": 0.13},
            ],
        }

        self.add_section(
            name="Per sample read length distribution",
            anchor="read_length_line_plot",
            description="Number of reads at each read length, per sample.",
            helptext="""
            Each line is one sample. The x-axis is read length and the y-axis 
            is the number of reads with that length.
            """,
            plot=linegraph.plot(data_by_sample, pconfig),
            # statuses=status_dict if section_statuses else None,
        )

    def read_length_violin_plot(self, read_len_data):
        """Generate an interactive violin plot of read lengths per sample"""

        # log10-transform (drop zero-length reads: log10(0) is undefined)
        lengths = {
            s: np.log10(v[v > 0])
            for s, v in ((s, np.asarray(vals, dtype=np.int64)) for s, vals in read_len_data.items())
        }

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
