import logging
from typing import Dict

# from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table
from multiqc.plots.table_object import ColumnDict, ValueT

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    MultiQC with parse outputs from matchbox
    """

    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name="matchbox",
            anchor="matchbox",
            href="https://github.com/jakob-schuster/matchbox",
            info="Read processor that matches and transforms reads.",
            doi="https://doi.org/10.1101/2025.11.09.685711",
        )

        # Find and load any matchbox reports
        matchbox_data: Dict[str, Dict[str, int]] = dict()
        for f in self.find_log_files("matchbox"):
            s_name = f["s_name"]
            matchbox_data[s_name] = self.parse_matchbox(f["f"])
            if s_name in matchbox_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f)
        # Report if no samples found
        if len(matchbox_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(matchbox_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Add matchbox summary to the general stats table
        self.matchbox_general_stats_table(matchbox_data)

        # Alignment Rate Plot
        # self.matchbox_alignment_plot()

    def parse_matchbox(self, f) -> Dict[str, int]:
        """Parse matchbox files"""

        parsed_data = {}

        for line in f.splitlines():
            s = line.strip().split(",")

            if s[0] != "value":
                parsed_data[s[0]] = int(s[1])

        return parsed_data

    def matchbox_general_stats_table(self, matchbox_data):
        """Take the parsed stats from the matchbox report and add it to the
        basic stats table at the top of the report"""

        headers = {
            "heavy": {
                "title": "Heavy chains",
                "description": "Number of heavy chains found",
                "min": 0,
                "scale": "OrRd",
            },
            "heavy + kappa": {
                "title": "Heavy + kappa light chains",
                "description": "Number of heavy and kappa light chains found",
                "min": 0,
                "scale": "Greens",
            },
            "heavy + lambda": {
                "title": "Heavy + lambda light chains",
                "description": "Number of heavy and lambda light chains found",
                "min": 0,
                "scale": "Greens",
            },
            "rotated": {
                "title": "Reads rotated",
                "description": "Number of reads that have been rotated",
                "min": 0,
                "scale": "BuPu",
            },
            "total reads": {
                "title": "Total reads",
                "description": "Total number of reads parsed",
                "min": 0,
                "scale": "Blues",
            },
        }

        self.general_stats_addcols(matchbox_data, headers)

    # def matchbox_func1(self):
    #     """Generate plot for the matchbox plot"""

    #     p_config = {"id": "mirtop_read_count_plot",
    #                 "title": "mirtop: IsomiR read counts",
    #                 "ylab": "Read counts"}

    #     self.add_section(
    #         name = "matchbox test section",
    #         anchor = "matchbox test",
    #         description = "Total counts of chains over all reads.",
    #         helptext = """
    #         Breakdown of total reads and heavy and light chains extracted.
    #         """,

    #         plot = bargraph.plot(self.filter_plot_data("sum"),
    #                              self.get_plot_cats("sum"),
    #                              p_config),
    #     )
