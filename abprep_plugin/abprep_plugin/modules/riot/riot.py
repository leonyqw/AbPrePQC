import csv
import logging
from typing import Dict

# from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table
from multiqc.plots.table_object import ColumnDict, ValueT

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    MultiQC with parse outputs from riot
    """

    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name="riot",
            anchor="riot",
            href="https://pypi.org/project/riot-na/",
            info="Rapid Immunoglobulin Overview Tool for antibody numbering.",
            doi="https://doi.org/10.1093/bib/bbae632",
        )

        # Find and load any riot data
        riot_data: Dict[str, Dict[str, int]] = dict()
        for f in self.find_log_files("riot", filehandles=True):
            s_name = f["s_name"]
            riot_data[s_name] = self.parse_riot(f["f"])
            if s_name in riot_data:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f)

        # Report if no samples found
        if len(riot_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(riot_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Add riot summary to the general stats table
        self.riot_general_stats_table(riot_data)

    def parse_riot(self, f) -> Dict[str, int]:
        """Parse riot files"""

        file = csv.reader(f)
        header = next(file)
        productive_idx = header.index("productive")

        productive_count = 0
        total = 0

        # Iterate and get productive counts from riot output
        for row in file:
            if row[productive_idx] == "True":
                productive_count += 1
            total += 1

        return {
            "productive": productive_count,
            "total": total,
            "productive_percent": ((productive_count / total) * 100) if total > 0 else 0,
        }

    def riot_general_stats_table(self, riot_data):
        """Take the parsed stats from the riot report and add it to the
        basic stats table at the top of the report"""

        # Create new dictionary for barcode and productive counts
        productive_data = {}

        for sample_chain in riot_data:
            barcode = sample_chain.split("_")[0]

            if barcode not in productive_data:
                productive_data[barcode] = {}

            if "heavy" in sample_chain:
                productive_data[barcode]["productive_heavy"] = riot_data[sample_chain]["productive_percent"]
            elif "light" in sample_chain:
                productive_data[barcode]["productive_light"] = riot_data[sample_chain]["productive_percent"]

        headers = {
            "productive_heavy": {
                "title": "Productive heavy chains",
                "description": "Percentage of productive heavy chains",
                "min": 0,
                "scale": "OrRd",
            },
            "productive_light": {
                "title": "Productive light chains",
                "description": "Percentage of productive light chains",
                "min": 0,
                "scale": "Greens",
            },
        }

        self.general_stats_addcols(productive_data, headers)


#     # def riot_func1(self):
#     #     """Generate plot for the riot plot"""

#     #     p_config = {"id": "mirtop_read_count_plot",
#     #                 "title": "mirtop: IsomiR read counts",
#     #                 "ylab": "Read counts"}

#     #     self.add_section(
#     #         name = "riot test section",
#     #         anchor = "riot test",
#     #         description = "Total counts of chains over all reads.",
#     #         helptext = """
#     #         Breakdown of total reads and heavy and light chains extracted.
#     #         """,

#     #         plot = bargraph.plot(self.filter_plot_data("sum"),
#     #                              self.get_plot_cats("sum"),
#     #                              p_config),
#     #     )
