import logging
import pandas as pd
from typing import List
import hashlib

# from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, table, heatmap

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

        # Find and load riot data
        riot_data: List[pd.DataFrame] = []
        report_num = 0

        for f in self.find_log_files("riot", filehandles=True):
            s_name = f["s_name"].split("_")[0]
            df = self.parse_riot(f)
            if not df.empty:
                report_num += 1
            df["barcode"] = s_name
            riot_data.append(df)

            # if s_name in riot_data:
            #     log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            # self.add_data_source(f)

        # Concatenate data into one dataframe
        riot_data = pd.concat(riot_data, ignore_index=True)

        # Report if no samples found
        if report_num == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {report_num} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version("4.0.2")

        # Remove duplicate sequence headers
        riot_data = riot_data.drop_duplicates(subset=["sequence_header", "barcode", "chain"], keep="first")

        # Join data by sequence header to pair heavy and light chains together for each read
        riot_wide_data = riot_data.pivot(
            index=["sequence_header", "barcode"], columns=["chain"], values=["v_call", "productive"]
        )

        # Rename columns
        riot_wide_data.columns = pd.MultiIndex.from_tuples(
            [
                ("v_call", "v_heavy"),
                ("v_call", "v_light"),
                ("productive", "prod_heavy"),
                ("productive", "prod_light"),
            ]
        )
        riot_wide_data.columns = riot_wide_data.columns.droplevel(0)
        riot_wide_data = riot_wide_data.reset_index()

        # Create new column (True / False) if read has both productive heavy and light chains
        riot_wide_data["prod_vh_vl"] = riot_wide_data["prod_heavy"] & riot_wide_data["prod_light"]

        # Summarise riot data
        riot_data_summary = self.summarise_data(riot_wide_data)

        # Add riot summary to the general stats table
        self.riot_general_stats_table(riot_data_summary)

        # Add full riot summary table to the report
        self.riot_summary_table(riot_data_summary)

        # Add riot V gene heat map to the report
        self.riot_v_gene_plots(riot_wide_data)

    def parse_riot(self, f):
        """Parse riot files"""

        chain_type = f["s_name"].split("_")[2]
        df = pd.read_csv(f["f"], usecols=["sequence_header", "productive", "v_call"])
        df["chain"] = chain_type

        return df

    def summarise_data(self, riot_data):
        """Summarise riot data"""

        # Aggregate the data by barcode
        riot_agg_data = riot_data.groupby(["barcode"]).agg(
            total=("sequence_header", "count"),
            productive_vh_vl=("prod_vh_vl", lambda x: (x).sum()),
            productive_heavy=("prod_heavy", lambda x: (x).sum()),
            unique_vh_genes=("v_heavy", lambda x: (x.nunique())),
            productive_light=("prod_light", lambda x: (x).sum()),
            unique_vl_genes=("v_light", lambda x: (x.nunique())),
        )

        # Calculate unproductive and productive percent
        riot_agg_data["unproductive_heavy"] = riot_agg_data["total"] - riot_agg_data["productive_heavy"]
        riot_agg_data["unproductive_light"] = riot_agg_data["total"] - riot_agg_data["productive_light"]
        riot_agg_data["productive_perc_heavy"] = (riot_agg_data["productive_heavy"] / riot_agg_data["total"]) * 100
        riot_agg_data["productive_perc_light"] = (riot_agg_data["productive_light"] / riot_agg_data["total"]) * 100
        riot_agg_data["productive_perc_vh_vl"] = (riot_agg_data["productive_vh_vl"] / riot_agg_data["total"]) * 100

        # Summarise the riot data
        riot_data_summary = {}
        for barcode, row in riot_agg_data.iterrows():
            riot_data_summary[barcode] = row.to_dict()

        return riot_data_summary

    def riot_general_stats_table(self, riot_data_summary):
        """Take the parsed summarised riot data and add it to the
        basic stats table at the top of the report"""

        headers = {
            "total": {
                "title": "Total VH-VL pairs",
                "description": "Total number of heavy and light chain pairs found",
                "min": 0,
                "max": 100,
                "format": "{:,.0f}",
                "scale": "Greens",
            },
            "productive_perc_vh_vl": {
                "title": "% Productive VH-VL pairs",
                "description": "Percentage of productive heavy and light chain pairs",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Greens",
            },
        }

        self.general_stats_addcols(riot_data_summary, headers)

    def riot_summary_table(self, riot_data_summary):
        """Generate riot data table"""

        p_config = {"id": "riot_productivity_plot", "title": "RIOT: Productivity", "xlab": "Read counts"}

        headers = {
            "total": {
                "title": "Total VH-VL pairs",
                "description": "Total number of heavy and light chain pairs found",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Greens",
            },
            "productive_vh_vl": {
                "title": "Productive VH-VL pairs",
                "description": "Total number of productive heavy and light chain pairs found",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Greens",
            },
            "productive_perc_vh_vl": {
                "title": "% Productive VH-VL pairs",
                "description": "Percentage of productive heavy and light chain pairs found",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Greens",
            },
            "productive_heavy": {
                "title": "Productive VH chains",
                "description": "Total number of productive heavy chains found",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Oranges",
            },
            "unproductive_heavy": {
                "title": "Unproductive VH chains",
                "description": "Total number of unproductive heavy chains found",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Oranges",
            },
            "productive_perc_heavy": {
                "title": "% Productive VH chains",
                "description": "Percentage of productive heavy chains found",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Oranges",
            },
            "unique_vh_genes": {
                "title": "Unique VH genes",
                "description": "Number of unique variable heavy genes found",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Oranges",
            },
            "productive_light": {
                "title": "Productive VL chains",
                "description": "Total number of productive light chains found",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Blues",
            },
            "unproductive_light": {
                "title": "Unproductive VL chains",
                "description": "Total number of unproductive light chains found",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Blues",
            },
            "productive_perc_light": {
                "title": "% Productive VL chains",
                "description": "Percentage of productive light chains found",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Blues",
            },
            "unique_vl_genes": {
                "title": "Unique VL genes",
                "description": "Number of unique variable light genes found",
                "min": 0,
                "format": "{:,.0f}",
                "scale": "Blues",
            },
        }

        self.add_section(
            name="RIOT: Productivity",
            anchor="riot_productivity",
            description="Number and percentage of productive heavy and light chains, and unique V genes identified.",
            comment="Low productive % may be due to sequencing error.",
            helptext="""
            Number and percentage of productive (no stop codons) and unproductive heavy and light chains.
            """,
            plot=table.plot(riot_data_summary, headers=headers, pconfig=p_config),
        )

    def riot_v_gene_plots(self, riot_wide_data):
        """Generate heat map for v gene pairing"""

        # Keep only rows with both vh and vl genes, and specific columns
        riot_wide_data = riot_wide_data[["barcode", "v_heavy", "v_light"]].dropna(subset=["v_heavy", "v_light"])

        # Clean up v gene name. Remove allele information
        for col in ["v_heavy", "v_light"]:
            riot_wide_data.loc[:, col] = riot_wide_data[col].str.split("*").str[0]

        # Create dictionary of VH gene and VL gene pair counts
        v_genes_dict = {}
        for (v_heavy, v_light), count in riot_wide_data[["v_heavy", "v_light"]].value_counts().items():
            v_genes_dict.setdefault(v_heavy, {})[v_light] = count

        # Create dictionary of vh gene counts per barcode
        vh_counts = pd.concat(
            [
                riot_wide_data.groupby("barcode")["v_heavy"].value_counts(),
                riot_wide_data.groupby("barcode")["v_heavy"].value_counts(normalize=True),
            ],
            axis=1,
            keys=["count", "prop"],
        ).reset_index("v_heavy")

        vh_counts["v_heavy"] = vh_counts["v_heavy"].mask((vh_counts["prop"] < 0.05), "other")
        vh_counts = vh_counts.groupby(["barcode", "v_heavy"])["count"].sum().reset_index()

        vh_counts_dict = {
            barcode: dict(zip(group["v_heavy"], group["count"])) for barcode, group in vh_counts.groupby("barcode")
        }

        # Create dictionary of vl gene counts per barcode
        vl_counts = pd.concat(
            [
                riot_wide_data.groupby("barcode")["v_light"].value_counts(),
                riot_wide_data.groupby("barcode")["v_light"].value_counts(normalize=True),
            ],
            axis=1,
            keys=["count", "prop"],
        ).reset_index("v_light")

        vl_counts["v_light"] = vl_counts["v_light"].mask((vl_counts["prop"] < 0.05), "other")
        vl_counts = vl_counts.groupby(["barcode", "v_light"])["count"].sum().reset_index()

        vl_counts_dict = {
            barcode: dict(zip(group["v_light"], group["count"])) for barcode, group in vl_counts.groupby("barcode")
        }

        # Assign different colours to each gene for the plot
        genes = sorted(set(vh_counts["v_heavy"]).union(set(vl_counts["v_light"])))
        gene_cats = {}

        for gene in genes:
            gene_cats[gene] = {"name": gene, "color": self.color_from_name(gene)}

        # Config for bar plots
        bar_config = {
            "id": "v_gene_plot",
            "title": "RIOT: Germline V gene composition",
            "xlab": "Counts",
            "data_labels": [
                {"name": "VH genes"},
                {"name": "VL genes"},
            ],
            "tt_decimals": 0,
            "tt_suffix": "reads",
            "cpswitch_c_active": False,
        }

        # Add V gene composition bar plots to the report
        self.add_section(
            name="RIOT: Germline V gene composition",
            anchor="riot_genes_composition",
            description="Germline V genes identified in the heavy and light chains.",
            helptext="""
            Number of each germline V gene identified in the heavy and light chains.

            The other category includes all genes that make up less than 5% of the total V gene counts for that sample.
            """,
            plot=bargraph.plot([vh_counts_dict, vl_counts_dict], cats=gene_cats, pconfig=bar_config),
        )

        # Config for heatmap plot
        heatmap_config = {
            "id": "v_gene_pairing_heatmap_plot",
            "title": "RIOT: Germline VH and VL gene pairing",
            "xlab": "VL gene",
            "ylab": "VH gene",
            "zlab": "Counts",
            "xcats_samples": False,
            "ycats_samples": False,
            "colstops": [
                [0, "#feedde"],
                [0.25, "#fdbe85"],
                [0.5, "#fd8d3c"],
                [0.75, "#e6550d"],
                [1, "#a63603"],
            ],
        }

        # Add V gene pairing heatmap to the report
        self.add_section(
            name="RIOT: Germline V gene pairing",
            anchor="riot_genes_pairing",
            description="Germline variable (V) gene pairing between heavy and light chains. Darker colours indicate greater abundance of that V gene pairing.",
            helptext="""
            Number of each germline V gene pairing between the heavy and light chains.

            The x-axis is the light chain V genes.
            The y-axis is the heavy chain V genes.
            """,
            plot=heatmap.plot(v_genes_dict, pconfig=heatmap_config),
        )

    def color_from_name(self, name):
        """Deterministic colour from a string — same name always gives same colour."""

        h = hashlib.md5(name.encode()).hexdigest()
        r = int(h[0:2], 16)
        g = int(h[2:4], 16)
        b = int(h[4:6], 16)

        # Avoid very light (invisible) or very dark colours
        r = max(50, min(210, r))
        g = max(50, min(210, g))
        b = max(50, min(210, b))

        # TODO delete useless code
        # # Map 0-255 into a lighter range, e.g. 100-220
        # def lighten(col):
        #     return 50 + int(col / 127.5 * 100)

        # # Convert the hash to a color
        # hue = int(h[0:4], 16) / 0xFFFF
        # r, g, b = colorsys.hls_to_rgb(hue, 0.75, .55)
        # r, g, b = int(r*255), int(g*255), int(b*255)

        # return
        # return f'#{lighten(r):02x}{lighten(g):02x}{lighten(b):02x}'
        return f"#{(r):02x}{(g):02x}{(b):02x}"
