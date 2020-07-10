# tharsis-small-vents
Code and data for Richardson et al. "Small volcanic vents of the Tharsis Volcanic Province, Mars"

## Data
The Enhanced Vent Catalog used for analysis in Richardson et al. is available in data/TharsisVentCatalog_20200606_enhanced.tsv. Each vent entry in this catalog is assigned geographic and morphologic information including:

 - FID - vent number in catalog.
 - Longitude - degrees (positive east).
 - Latitude - degrees (positive north).
 - Elevation - MOLA gridded product elevation at Longitude, Latitude location.
 - v1 lat, v1 long, v2 lat, v2 long - latitude and longitude coordinates for the end-points of vent major axis.
 - Width - measured minor axis distance of vent, meters.
 - peak east, peak north - location of vent-adjacent peak, meters easting/northing using the Mars 2000 equirectangular projection.
 - peak elev - MOLA gridded product elevation at vent-adjacent peak location.
 - noVentInt, noVent - If a vent is observed, these values are 0, \Vent". If no vent is observed, this entry is a likely vent and values are 1, \No Vent."
 - Orientation - orientation of the vent major axis, degrees from north.
 - Length - major axis length, meters.
 - AspectRatio - Length/Width.
 - EquantInt, Equant - If vent aspect ratio is 1.5, values are 1, \Equant." If vent aspect ratio is >1.5, values are 0, \NotEquant."
 - Prominence - prominence of vent-adjacent peak, meters .
 - clusterID - numeral ID for Tharsis region, from clustering algorithm (see text).
 - clusterName - name of Tharsis region (e.g. Olympus).

Note, all fields after and including 'Orientation' are produced in the following python scripts. Fields before this in the catalog are measured and are included in data/TharsisVentCatalog_20200606.tsv

## Code

#### Setup
Before running tharsis_vent_analysis.py, expand MOLA.tar.bz2 and ensure MOLA_128ppd_Tharsis.tif is a file within the data subdirectory.

#### tharsis_vent_analysis.py
This script takes in an ascii catalog of observed vents and some measurements (e.g., data/TharsisVentCatalog_200sample.tsv) and returns an enhanced database with derived fields.

#### tharsis_intervent_analysis.py
This script takes an enhanced database from tharsis_vent_analysis.py and identifies intervent alignments between nearby vent features.

#### cluster_analysis_tools.py
This script holds functions that are used by the previous two python scripts.
