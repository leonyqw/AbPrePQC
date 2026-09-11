#!/usr/bin/env python
"""AbPreP plugin

We can add any custom Python functions here and call them
using the setuptools plugin hooks.

File to hold custom functions that can tie into the main MultiQC execution flow.
In this file, we define some new config defaults, including the search patterns used by the example module
"""

# from __future__ import print_function

from multiqc import config
import importlib_metadata
import logging

# Initialise the main MultiQC logger
log = logging.getLogger("multiqc")


# Add default config options for the things that are used in MultiQC_NGI
def abprep_plugin_execution_start():
    """Code to execute after the config files and command line flags have been parsed."""

    # Plugin's version number defined in pyproject.toml:
    version = importlib_metadata.version("multiqc_abprep_plugin")
    log.info(f"Running AbPreP MultiQC Plugin v{version}")

    # Add to the main MultiQC config object.

    # Add to the search patterns used by modules
    if "matchbox" not in config.sp:
        config.update_dict(config.sp, {"matchbox": {"fn": "*_matchbox_counts.csv", "num_lines": 10}})

    if "riot" not in config.sp:
        config.update_dict(config.sp, {"riot": {"fn": "*_annot_*.csv"}})

    if "read_length" not in config.sp:
        config.update_dict(config.sp, {"read_length": {"fn": "*read_lengths.tsv"}})

    # # Some additional filename cleaning
    # config.fn_clean_exts.extend([".my_tool_extension", ".removeMetoo"])
