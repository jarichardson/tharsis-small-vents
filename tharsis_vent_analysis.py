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
topofile       = 'data/MOLA_128ppd_Tharsis.tif'
ventCatalog    = 'data/TharsisVentCatalog_20200606.tsv'
outputCatalog  = 'data/TharsisVentCatalog_20200606_enhanced.tsv'
azimuthCatalog = 'data/Azimuths_20200606.tsv'

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
  if aspectRatio < 1.0:  
      '''If length is shorter than width, major/minor axes are mismapped. 
      Aspect ratio is always defined as major axis / minor axis.'''
      aspectRatio = row.Width/length
  orientation = capy.Bearing(row.v1_long,row.v1_lat,row.v2_long,row.v2_lat)
  #record orientation as +/-90 degrees from North
  if abs(orientation) > 90.0: orientation += np.sign(orientation)*(-180.0);

  ### Record these values in DataFrame
  ventData.at[idx, 'Length'] = length
  ventData.at[idx, 'AspectRatio'] = aspectRatio

  # Is Vent Equant? If not, record the orientation
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


clusterNames = ["Olympus", "NE Ascreaus", "Ceraunius", "East Pavonis", "Arsia", \
                "Tempe-Mareotis", "Syria Planum", "Ulysses", "Daedalia Planum", \
                "Fortuna", "Labeatis Mons"]
maxClusterSize = 600000 #600 km radius

#Find clusters and assign vents ClusterID Numbers
ventLocs = ventData[['Longitude','Latitude']].to_numpy()
ventData['clusterID'] = capy.doCluster(ventLocs,thresh=maxClusterSize,color='color')
# Assign cluster names
for idx, row in ventData.iterrows():
    ventData.at[idx, 'clusterName'] = clusterNames[int(ventData.at[idx,'clusterID'])]
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


sys.stdout.write('\n All Prominences Calculated.\n')
sys.stdout.write('\nTotal executed time: %0.2f seconds\n\n' % (time.time() - startTime))


sys.stdout.write('Identifying Intervent alignments.\n\n')

### Intervent Analysis ###

### Cluster Specific Analysis



# Set up a Two-Point Azimuth Data Frame        
TPAdf = pd.DataFrame({ 'clusterID' : 1000,
                                'Bearing' : 0.0,
                                'Lat1' : 0.0,
                                'Lon1' : 0.0,
                                'Lat2' : 0.0,
                                'Lon2' : 0.0
                            }, index=[0])

ventData['clusterID']
for C in np.unique(ventData['clusterID']).astype(int):
  cVents = ventData[ventData['clusterID'] == C]
  ptCt = len(cVents)
  print("%s Cluster has %d vents" % (clusterNames[C],len(cVents)))
  
  ### Two-Point Azimuth Analysis
  ptMin = 4
  #distance matrix
  dCt = int(0.5 * ptCt * (ptCt-1))
  dMatrix = np.zeros(dCt)
  azims = np.zeros((dCt,5))
  
  if ptCt > ptMin:
    d = 0
    ptLocs = np.asarray([cVents['Longitude'],cVents['Latitude']]).T
    for p in range(ptCt):
      p1 = ptLocs[p]
      for r in np.arange(p+1,ptCt):
        p2 = ptLocs[r]
        
        #Calculate distance using spheroid
        dMatrix[d] = capy.haversineM(p1[0],p1[1],p2[0],p2[1])
        azims[d] = [p1[0],p1[1],p2[0],p2[1],
             capy.Bearing(p1[0],p1[1],p2[0],p2[1])]
        d+=1
    
    CebriaLength = (np.mean(dMatrix) - np.std(dMatrix))/3
    LocalDists = dMatrix[np.where(dMatrix < CebriaLength)]
    LocalAzims = azims[np.where(dMatrix < CebriaLength)]
    
    ### Adjust azimuths to keep in North semicircle of rose plots
    for v in LocalAzims:
      if v[4] >= 90.0:
        v[4] += 180.0 #add additional 0.01 to keep the azimuths within the North semicircle of windrose plot
      elif v[4] <= -90.0:
        v[4] += 179.99
      if v[4] < 0:
        v[4] += 360
  
    TPA = pd.DataFrame({ 'clusterID' : C,
                                'Bearing' : LocalAzims[:,4],
                                'Lon1' : LocalAzims[:,0],
                                'Lat1' : LocalAzims[:,1],
                                'Lon2' : LocalAzims[:,2],
                                'Lat2' : LocalAzims[:,3]
                            })
    
    print("  %d Azimuths are less than Cebria length: %0.3f m" % (len(TPA), CebriaLength))
    
    # Concatenate TPA dataframes
    dfs = [TPAdf, TPA]
    TPAdf = pd.concat(dfs)
    
    histazims = np.histogram(LocalAzims[:,4],bins=np.arange(0,360,20))
    print('  Azimuths     (0-360)      :  ',histazims[0])
    
  orients = np.asarray([cVents['Orientation']])
  orients = orients[np.where(np.isfinite(orients))] #remove nans
  historients = np.histogram(orients,bins=np.arange(-180,180,20))
  print('  Orientations (-180 to 180): ', historients[0])
    
#Get rid of first (placeholder) entry in TPA results.
TPAdf = TPAdf[TPAdf['clusterID']<500]

### Write Out all enhanced results!

# Write out enhanced vent information
ventData.to_csv(outputCatalog, sep="\t", float_format='%g')
sys.stdout.write('\nWrote out all data to enhanced catalog at: %s\n' % outputCatalog)


# Write Two Point Azimuth Results to csv
TPAdf.to_csv(azimuthCatalog, sep="\t")
print("Two-point azimuth results written to file: %s" % azimuthCatalog)

sys.stdout.write('\n\nTotal executed time: %0.2f seconds\n\n' % (time.time() - startTime))

