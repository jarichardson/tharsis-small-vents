# -*- coding: utf-8 -*-
"""
Volcanic vent analysis for Richardson et al, "Small volcanic vents of the Tharsis Volcanic Province, Mars"
Code by J Richardson

Python 3.7.6
requires numpy, pandas, sys, and cluster_analysis_tools (should be in same directory)

Calculate Mars Vent Data:
  Vent Bearing
  Vent Aspect Ratio
  Prominence
  Cluster ID


  ********************

Inputs:
  a csv file (vents) that contains:
    Vent ID (unique number) (updated to remove missing/deleted VentID numbers?)
    Long
    Lat
    Easting
    Northing
    Elevation
    long axis end coordinate pairs (lat, long)
    peak location coordinate pait (lat long)
    vent/no vent toggle

Outputs:
  a csv file with statistics for each point and all input information:
    Vent ID (unique number) (updated to remove missing/deleted VentID numbers?)
    Lat
    Long
    Easting
    Northing
    VentLeft (-999 if no vent, Long if equant)
    VentRight (-999 if no vent, Long if equant)
    VentNorth (-999 if no vent, Lat if equant)
    VentSouth (-999 if no vent, Lat if equant)
    Orientation
    Cluster Number
    Cluster Name
    Topographic Prominence

Cataloged small volcanoes have an X,Y location. All small volcanoes will have a
mapped prominence, even one that is 0.

Cataloged small volcanoes might not have observed vents. Observed vent lineations
might be present for vents that do not exist anymore (they've been consolidated)
"No vent" vent orientation lines have already been removed (in ArcMap).
"Equant" vents are listed in a file, by giving ventID.

"""

import cluster_analysis_tools as capy
import numpy as np
import pandas as pd
import sys, time

startTime = time.time()

### Data Files
topofile      = 'data/MOLA_128ppd_Tharsis.tif'
ventCatalog   = 'data/TharsisVentCatalog_200sample.tsv'
outputCatalog = 'data/TharsisVentCatalog_200sample_enhanced.tsv'

### LOAD DATA
sys.stdout.write("\nLoading Vent Data\n")
points = pd.read_csv(ventCatalog,sep='\t')
sys.stdout.write('Loaded %d Vents\n\n' % len(points))


### ASSEMBLE VENT DATA with points and vent lines

floatColumn = [np.nan]*len(points)
intColumn = [0]*len(points)
strColumn = ['']*len(points)

suppData = pd.DataFrame({'Orientation' : floatColumn,
                         'Length' : floatColumn,
                         'AspectRatio' : floatColumn,
                         'EquantInt' : intColumn,
                         'Equant' : strColumn,
                         'Prominence' : intColumn,
                         'clusterID' : intColumn,
                         'clusterName' : strColumn
                                    })

ventData = points.join(suppData)


### Vent-wise analysis
for idx, row in ventData.iterrows():
  ### Remove vent line information for NoVent features (noVentInt)
  if row.noVentInt: #if this value is not zero
    ventData.at[idx, 'v1_lat'] = np.nan
    ventData.at[idx, 'v1_long'] = np.nan
    ventData.at[idx, 'v2_lat'] = np.nan
    ventData.at[idx, 'v2_long'] = np.nan
    ventData.at[idx, 'Width'] = np.nan
    continue

  ### Calculate Vent Length, Aspect Ratio, Orientation
  length = capy.haversineM(row.v1_long,row.v1_lat,row.v2_long,row.v2_lat)
  aspectRatio = length/row.Width
  orientation = capy.Bearing(row.v1_long,row.v1_lat,row.v2_long,row.v2_lat)
  #record orientation as +/-90 degrees from North
  if abs(orientation) > 90.0: orientation += np.sign(orientation)*(-180.0);

  ### Record these values in DataFrame
  ventData.at[idx, 'Length'] = length
  ventData.at[idx, 'AspectRatio'] = aspectRatio

  # Is Vent Equant? If so, record the orientation
  if aspectRatio < 1.5:
    ventData.at[idx, 'EquantInt'] = 1
    ventData.at[idx, 'Equant']    = "Equant"
  else:
    ventData.at[idx, 'EquantInt']   = 0
    ventData.at[idx, 'Equant']      = "NotEquant"
    ventData.at[idx, 'Orientation'] = orientation

### Write out some statistics!
sys.stdout.write('Found all vent lengths, orientations, and aspect ratios\n')
#equantCt = np.sum(ventData['EquantInt'].astype(int))

sys.stdout.write('Statistics Summary:\n')
sys.stdout.write('  Length [min, max]:       %7.2f m,  %7.2f m\n' %
                 (ventData['Length'].min(),ventData['Length'].max()))
sys.stdout.write('  Width [min, max]:        %7.2f m,   %7.2f m\n' %
                 (ventData['Width'].min(),ventData['Width'].max()))
sys.stdout.write('  Aspect Ratio [min, max]: %7.2f,     %7.2f\n' %
                 (ventData['AspectRatio'].min(),ventData['AspectRatio'].max()))
sys.stdout.write('  Orientation [min, max]:  %7.2f deg, %7.2f deg\n' %
                 (ventData['Orientation'].min(),ventData['Orientation'].max()))
sys.stdout.write('  All Vents:      %5d\n' % int(ventData['EquantInt'].count()))
sys.stdout.write('  Equant Vents:   %5d\n' % int(ventData['EquantInt'].sum()))
sys.stdout.write('  "Likely" Vents: %5d\n' % int(ventData['noVentInt'].sum()))

sys.stdout.write('\n\nTotal executed time: %0.2f seconds\n\n' % (time.time() - startTime))


### Find Clusters
sys.stdout.write('Finding Vent clusters\n' )
maxClusterSize = 600000 #600 km radius

ventLocs = ventData[['Longitude','Latitude']].to_numpy()

ventData['clusterID'] = capy.doCluster(ventLocs,thresh=maxClusterSize,color='gray')
ventData['clusterName'] = capy.doCluster(ventLocs,thresh=maxClusterSize,color='gray')
sys.stdout.write('  Identified %d Clusters.\n' % max(ventData['clusterID']+1))

sys.stdout.write('\n\nTotal executed time: %0.2f seconds\n\n' % (time.time() - startTime))


### Find Prominences
# Load Topo data for Prominences
sys.stdout.write("\nLoading Topographic Data\n")
(topo,geoTrans) = capy.gdalLoad(topofile)
sys.stdout.write('Loaded topography data\n\n')


sys.stdout.write('Finding Vent Prominences\n' )
# Convert peak locations to raster locations
peaklocs = ventData[['peak_east','peak_north']].to_numpy()
# Easting to Column
peaklocs[:,0] = (peaklocs[:,0]-geoTrans[0])/geoTrans[1]
# Northing to Row
peaklocs[:,1] = (peaklocs[:,1]-geoTrans[3])/geoTrans[5]
#Reduce to integers
peaklocs = peaklocs.astype('int')

#Calculate Prominences
ventData['Prominence'] = capy.getAllProms(peaklocs,topo)


### Write out everything!
ventData.to_csv(outputCatalog, sep="\t")
sys.stdout.write('\nWrote out all data to enhanced catalog at: %s\n\n' % outputCatalog)

sys.stdout.write('\n\nTotal executed time: %0.2f seconds\n\n' % (time.time() - startTime))
