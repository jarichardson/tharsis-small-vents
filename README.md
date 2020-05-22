# tharsis-small-vents
Code and data for Richardson et al. "Small volcanic vents of the Tharsis Volcanic Province, Mars"

## Setup
Before running tharsis_vent_analysis.py, expand MOLA.tar.bz2 and ensure MOLA_128ppd_Tharsis.tif is a file within the data subdirectory.

## tharsis_vent_analysis.py
This script takes in an ascii catalog of observed vents and some measurements (e.g., data/TharsisVentCatalog_200sample.tsv) and returns an enhanced database with derived fields.

## tharsis_intervent_analysis.py
This script takes an enhanced database from tharsis_vent_analysis.py and identifies intervent alignments between nearby vent features.

## cluster_analysis_tools.py
This script holds functions that are used by the previous two python scripts.
