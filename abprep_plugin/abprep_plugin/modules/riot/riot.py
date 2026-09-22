import csv
import logging
from typing import Dict, Union

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
            name="RIOT",
            anchor="riot",
            href="https://pypi.org/project/riot-na/",
            info="RIOT - Rapid Immunoglobulin Overview Tool for antibody numbering.",
            doi="https://doi.org/10.1093/bib/bbae632",
        )

        # Find and load any riot data
        riot_data: Dict[str, Dict[str, Union[int, float]]] = dict()

        for f in self.find_log_files("riot", filehandles=True):
            s_name = f["s_name"].split("_")[0]
            riot_data.setdefault(s_name, {}).update(self.parse_riot(f))

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

        # Add riot section to the report
        self.riot_plot(riot_data)

    def parse_riot(self, f) -> Dict[str, Union[int, float]]:
        """Parse riot files"""

        file = csv.reader(f["f"])
        header = next(file)
        productive_idx = header.index("productive")
        chain_type = f["s_name"].split("_")[2]

        productive_count = 0
        total = 0

        # Iterate and get productive counts from riot output
        for row in file:
            if row[productive_idx] == "True":
                productive_count += 1
            total += 1

        return {
            "total": total,
            (chain_type + "_productive"): productive_count,
            (chain_type + "_unproductive"): (total - productive_count),
            (chain_type + "_productive_percent"): ((productive_count / total) * 100) if total > 0 else 0,
            (chain_type + "_unproductive_percent"): 100 - (((productive_count / total) * 100) if total > 0 else 0),
        }

    def riot_general_stats_table(self, riot_data):
        """Take the parsed stats from the riot report and add it to the
        basic stats table at the top of the report"""

        # # Create new dictionary for barcode and productive counts
        # productive_data = {}

        # for sample_chain in riot_data:
        #     barcode = sample_chain.split("_")[0]

        #     if barcode not in productive_data:
        #         productive_data[barcode] = {}

        #     if "heavy" in sample_chain:
        #         productive_data[barcode]["productive_heavy"] = riot_data[sample_chain]["productive_percent"]
        #     elif "light" in sample_chain:
        #         productive_data[barcode]["productive_light"] = riot_data[sample_chain]["productive_percent"]

        headers = {
            "heavy_productive": {
                "title": "Productive heavy chains",
                "description": "Percentage of productive heavy chains",
                "min": 0,
                "suffix": "%",
                "scale": "OrRd",
            },
            "light_productive": {
                "title": "Productive light chains",
                "description": "Percentage of productive light chains",
                "min": 0,
                "suffix": "%",
                "scale": "Greens",
            },
        }

        self.general_stats_addcols(riot_data, headers)

    def riot_plot(self, riot_data):
        """Generate plot for the riot plot"""

        p_config = {"id": "riot_productivity_plot", "title": "Riot TITLE", "xlab": "Read counts"}

        headers = {
            "total": {
                "title": "Total reads",
                "description": "Total heavy and light chain pairs found",
                "min": 0,
                # "format": "{:,.0f}",  # No decimal places please
            },
            "heavy_productive": {
                "title": "Productive heavy chains",
                "description": "Total number of productive heavy chains found",
                "min": 0,
            },
            "heavy_unproductive": {
                "title": "Unproductive heavy chains",
                "description": "Total number of unproductive heavy chains found",
                "min": 0,
            },
            "heavy_productive_percent": {
                "title": "Productive heavy chains (%)",
                "description": "Percentage of productive heavy chains found",
                "min": 0,
                "suffix": "%",
            },
            "light_productive": {
                "title": "Productive light chains",
                "description": "Total number of productive light chains found",
                "min": 0,
            },
            "light_unproductive": {
                "title": "Unproductive light chains",
                "description": "Total number of unproductive light chains found",
                "min": 0,
            },
            "light_productive_percent": {
                "title": "Productive light chains (%)",
                "description": "Percentage of productive light chains found",
                "min": 0,
                "suffix": "%",
            },
        }

        self.add_section(
            name="RIOT: productivity",
            anchor="riot_productivity",
            description="Number and percentage of productive heavy and light chains.",
            helptext="""
            Number and percentage of productive heavy and light chains (no stop codons).
            """,
            plot=table.plot(riot_data, headers=headers, pconfig=p_config),
        )
