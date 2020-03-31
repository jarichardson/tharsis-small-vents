# -*- coding: utf-8 -*-
"""
Intervent alignment analysis for Richardson et al, "Small volcanic vents of the Tharsis Volcanic Province, Mars"
Code by J Richardson

Python 3.7.6
requires numpy, pandas, sys, and cluster_analysis_tools (should be in same directory)

Input: Enhanced Vent database, output from "tharsis_vent_analysis.py"
Output: A list of intervent alignments following Cebria et al 2009 methodology.
Also provides some histogrammable statistics.

"""

import cluster_analysis_tools as capy
import numpy as np
import pandas as pd
import sys

vent_db = 'C:/Users/jaricha4/Documents/research/tharsis_catalog/vent_catalog/TharsisVentCatalog_20190412_enhanced.tsv'
output_azim_db = 'C:/Users/jaricha4/Documents/research/tharsis_catalog/vent_catalog/Azimuths_20190412.tsv'
ventData = pd.read_csv(vent_db, sep='\t')

### Cluster Specific Analysis

clusterNames = ["Olympus", "Ceraunius", "East Pavonis", "NE Ascreaus", "Arsia", "Tempe-Mareotis",
        "North Uranius", "Syria Planum", "Ulysses", "Daedalia Planum", "Fortuna",
        "Labeatis Mons", "East Tharsis"]
        
TPAdf = pd.DataFrame({ 'clusterID' : 0,
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
  print "%s Cluster has %d vents" % (clusterNames[C],len(cVents))
  
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
    
    print("  %d Azimuths are less than Cebria length: %0.3f km" % (len(TPA), CebriaLength))
    
    # Concatenate TPA dataframes
    dfs = [TPAdf, TPA]
    TPAdf = pd.concat(dfs)
    
    histazims = np.histogram(LocalAzims[:,4],bins=np.arange(0,360,20))
    print '  Azimuths:     ',histazims[0]

    orients = np.asarray([cVents['Prominence']])
    orients = orients[np.where(np.isfinite(orients))] #remove nans
    historients = np.histogram(orients,bins=np.arange(0,360,20))
    print '  Orientations: ', historients[0]

    logprom = np.log10(np.asarray([cVents['Prominence']]))
    logprom = logprom[np.where(np.isfinite(logprom))] #remove -infs
    histproms = np.histogram(logprom,bins=np.arange(0,4,0.25))
    print '  Prominences: ', histproms[0]
    
    
TPAdf.drop(TPAdf.index[0]) #Get rid of first (placeholder) entry in TPA results.
# Write Two Point Azimuth Results to csv
TPAdf.to_csv(output_azim_db, sep='\t')

print("Two-point azimuth results written to file")