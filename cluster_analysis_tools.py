# -*- coding: utf-8 -*-
"""
Cluster Analysis Tools for Richardson et al, "Small volcanic vents of the Tharsis Volcanic Province, Mars"
Code by J Richardson

Python 3.7.6
requires numpy, gdal, matplotlib, scipy, cv2

Cluster Analysis Tools provides tools to describe the morphology and 
arrangement of distributed volcanic vents in a large volcanic province on Mars,
the Tharsis Volcanic Province. Tools might need to be modified for other planets
or point data sets.
"""

def Bearing(lon1, lat1, lon2, lat2):
  '''
  Bearing gives the initial bearing from pt 1 to pt 2, along a great circle of 
  a sphere. Initial Bearing returned in Clockwise Degrees from North. 
  -999 is returned if pts 1 and 2 are the same.
  '''

  import numpy as np
  if lon1 == lon2 and lat1 == lat2:
    return -999
  
  lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
  y = np.sin(lon2-lon1) * np.cos(lat2)
  x = np.cos(lat1)*np.sin(lat2) -\
        np.sin(lat1)*np.cos(lat2)*np.cos(lon2-lon1)
  brng = np.degrees(np.arctan2(y, x))
  return brng

def testBearing():
  '''Runs some tests to make sure Bearing is working right.'''
  print(Bearing(0,0,0,0), "-999")
  print(Bearing(0,0,0,1), "0")
  print(Bearing(0,0,1,0), "90")
  print(Bearing(0,0,0,-1), "180")
  print(Bearing(0,0,-1,0), "-90")
  print(Bearing(0,0,1,1), "45ish")
  print(Bearing(0,0,1,-1), "135ish")
  print(Bearing(0,0,-1,1), "-45ish")
  print(Bearing(0,0,-1,-1), "-135ish")
 

def haversineM(lon1, lat1, lon2, lat2):
    '''
    Calculate the great circle distance between two points 
    on Mars (specified in decimal degrees)
    '''
    import numpy as np
    
    radius = 3390000 #mars radius, m
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = np.sin(dlat/2)**2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon/2)**2
    c = 2 * np.arcsin(a**0.5) 
    m = radius * c
    return m

    
def doCluster(pts,thresh=600,color='color',treefig='',mapfig=''):
    '''
    Clusters a group of lon-lat points [pts] using the UPGMA algorithm. Clusters 
    smaller than [thresh] km are identified. Dendrogram is saved as figure if 
    [treefig] is defined. Map is saved if [mapfig] is defined. Returns 1D array 
    of cluster numbers corresponding to each input point of [pts]. Max # of clusters
    is # of steps in color chosen (use color=gray for up to 20 clusters).
    
    Color options: gray, 20 gray steps separated by 5%. 
                   by, 8 blue to yellow steps.
                   color, 11 simple colors.
    '''
    
    from scipy.cluster import hierarchy
    import matplotlib.pyplot as plt
    import numpy as np

    #Assemble Distance Matrix
    ptCt = len(pts)
    dCt = int(0.5 * ptCt * (ptCt-1))
    dMatrix = np.zeros(dCt)
    d = 0
    for p in range(ptCt):
        p1 = pts[p]
        for r in np.arange(p+1,ptCt):
            p2 = pts[r]
            #Calculate distance in some way (euclidean here, do spheroid next)
            #dMatrix[d] = ((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)**0.5
            dMatrix[d] = haversineM(p1[0],p1[1],p2[0],p2[1])
            d+=1
    
    #Hierarchy linkage with UPGMA algorithm
    Z = hierarchy.linkage(dMatrix, 'average')
    
    plt.figure()
    
    #Define Cluster Colors
    if color=='gray':
        #20 group, 0.05 Gray Color Palette
        hierarchy.set_link_color_palette(['0.0','0.05','0.1','0.15','0.2','0.25',
                                          '0.3','0.35','0.4','0.45','0.5','0.55',
                                          '0.6','0.65','0.7','0.75','0.8','0.85',
                                          '0.9','0.95'])
    elif color=='lgray':
        #20 group, 0.05 Gray Color Palette
        hierarchy.set_link_color_palette(['0.0','0.1','0.2',
                                          '0.3','0.4','0.5',  '0.6','0.7','0.8',
                                          '0.9','0.95','0.15','0.35','0.55','0.75'])
    elif color=='by':
        #Yellow Blue Color Palette
        hierarchy.set_link_color_palette(['#0743CA','#2357C1','#3f6bb9','#5b7fb0',
                                          '#7794a8','#93a89f','#afbc97','#cbd18f'])
    elif color=='color':
        #Colorful Color Palette
        hierarchy.set_link_color_palette(['k','orange','b','g','m','r','c','y',
                                          'peru','indigo','lime'])
    else:
        print("Color not found (options: gray, yb, color). Using 'color'.")
        #Colorful Color Palette
        hierarchy.set_link_color_palette(['k','orange','b','g','m','r','c','y',
                                          'peru','indigo','lime'])
    
    #Draw Dendrogram    
    dn = hierarchy.dendrogram(Z, distance_sort='ascending', orientation='top',
                              above_threshold_color='red',color_threshold=thresh,
                              no_labels=True)
    #Save Dendrogram Figure
    if treefig != '':
        plt.savefig(treefig)
    plt.show()
    
    #find point color corresponding to dendrogram
    leafColor = [None] * ptCt
    
    xlink = np.asarray(dn['icoord'])
    ylink = np.asarray(dn['dcoord'])
    for i in range(len(xlink)):
        #find if this is a line to a leaf
        if (xlink[i,0]-5.0)%10.0==0 and ylink[i,0]==0:
            leafColor[dn['leaves'][int((xlink[i,0]-5)/10)]] = dn['color_list'][i]
        if (xlink[i,3]-5.0)%10.0==0 and ylink[i,3]==0:
            leafColor[dn['leaves'][int((xlink[i,3]-5)/10)]] = dn['color_list'][i]
    
    #Map Points with Cluster Colors
    plt.clf()
    plt.scatter(pts[:,0], pts[:,1], s=15, c=leafColor, lw=0)
    if mapfig != '':
        plt.savefig(mapfig)
    plt.show()
    
    #create list of groups and assign group number to point
    colorInd = []
    groupList = np.zeros(len(pts))
    for i,c in enumerate(leafColor):
        try:
            groupList[i] = colorInd.index(c)
        except ValueError:
            colorInd.append(c)
            groupList[i] = colorInd.index(c)
    
    return groupList
    
def gdalLoad(rasterfile):
  '''
  Loads a Raster file into a 2D data array, using the GDAL API
  '''
  
  import gdal
  from gdalconst import GA_ReadOnly#, GDT_Float32
  #import struct
  import numpy as np
  
  dataset = gdal.Open( rasterfile, GA_ReadOnly )
  if dataset is None:
    print("No File Found")
  else:
    print('Driver: ', dataset.GetDriver().ShortName,'/', \
          dataset.GetDriver().LongName)
    print('Size is ',dataset.RasterXSize,'x',dataset.RasterYSize)
  
  #print 'Projection is ',dataset.GetProjection()
  geotransform = dataset.GetGeoTransform()
  if not geotransform is None:
      print('Origin = (',geotransform[0], ',',geotransform[3],')')
      print('Pixel Size = (',geotransform[1], ',',geotransform[5],')')
  

  band = dataset.GetRasterBand(1)
  '''
  min = band.GetMinimum()
  max = band.GetMaximum()
  if min is None or max is None:
      (min,max) = band.ComputeRasterMinMax(1)
  (minR,maxR) = band.ComputeRasterMinMax(1)
  print 'Min=%.3f, Max=%.3f' % (minR,maxR)
  '''

  #LOAD DATA INTO AN ARRAY
  
  data = band.ReadAsArray(0, 0, dataset.RasterXSize, dataset.RasterYSize).astype(np.float)
  '''
  data = np.zeros((dataset.RasterXSize,dataset.RasterYSize))

  for i in range(dataset.RasterYSize):
    scanline = band.ReadRaster( 0, 0, band.XSize, (i+1), \
                               band.XSize, 1, GDT_Float32 )
    line = np.asarray(struct.unpack('f' * band.XSize, scanline))
    data[:,i] = line
  '''
  metadata = np.array([geotransform[0],geotransform[1], dataset.RasterXSize,
                      geotransform[3], dataset.RasterYSize ,geotransform[5]])

  return data, metadata

def peakLoc(pt,topo,geoTrans):
    '''
    Finds a peak nearby a given point on a Topo map. This is created to find a
    summit for a point, when the point might be 1 to n grid cells off.
    '''
    
    import numpy as np
    ptCell = [int((pt[0]-geoTrans[0])/geoTrans[1]),int((pt[1]-geoTrans[3])/geoTrans[5])]
    #topo[ptCell[1],ptCell[0]]
    n = 4 # This is for 2 km radius on ~500 m MOLA topography on Mars
    n = 10 # This is for 0.1 km radius on 10 m SRTM topography on Earth...
    peakNeigh = np.meshgrid(np.arange(ptCell[1]-n,ptCell[1]+n+1),np.arange(ptCell[0]-n,ptCell[0]+n+1))
    peakLoc = np.where(topo[peakNeigh] == np.max(topo[peakNeigh]))
    if len(peakLoc[0])>1:
        #find peak closest to vent, if there is more than one pixel with peak height
        dist = np.zeros(len(peakLoc[0]))
        for i in range(len(peakLoc[0])):
            dist[i] = ((peakLoc[0][i]-n)**2+(peakLoc[1][i]-n)**2)**0.5
        peakLoc = [peakLoc[0][np.where(dist==np.min(dist))],peakLoc[1][np.where(dist==np.min(dist))]]
    
    return [peakNeigh[0][peakLoc[0][0],peakLoc[1][0]], peakNeigh[1][peakLoc[0][0],peakLoc[1][0]]],np.max(topo[peakNeigh])

def findProminence(topoWindow,cElev):
    '''
    Within a topographic map, find the prominence of the center location. 
    This assumes there is a higher point in the map
    '''
    import cv2
    import numpy as np
    c = int(topoWindow.shape[0]/2)
    hip = 0  
    prom = 0
    while hip == 0:
        #if prom<50:
        prom += 1
        #else:
        #    prom += 10
        conLvl = cElev - prom
        frame = np.copy(topoWindow)
        frame[np.where(frame<conLvl)] = 0
        frame[np.where(frame>=conLvl)] = 255
        frame = frame.astype(np.uint8)      
        
        higher = np.where(topoWindow > cElev)
        
        contours, hierarchy = cv2.findContours(frame,cv2.RETR_TREE,cv2.CHAIN_APPROX_SIMPLE)
        
        for ctr in contours:
            if cv2.pointPolygonTest(ctr,(c,c),0)>=0: #if center location is in the contour
                for i in range(len(higher[0])): #Check to see if a higher location is also in the contour
                    if cv2.pointPolygonTest(ctr,(higher[0][i],higher[1][i]),0)==0:
                        hip = 1
                        break
        if prom > 30000:
            break
    
    return prom


def findPeaks(pts, topoFile, peakNeighborhood, peakMinHeight, peakMaxRadius):
  '''
  Use a topographic grid and point locations to identify peak locations 
  corresponding to points.
  
  Inputs:    
    pts, topoFile, neighborhood, relief threshold,  search radius
    
  pts is Nx2 [easting, northing] array
  
  topoFile is an elevation raster file
  
  Output:
    list of peak locations and prominence (long, lat, prominence), corresponding
    to input point list. 0, 0, 0 means no peak was found
  '''
  
  import sys
  import numpy as np
  from scipy.spatial import distance_matrix
  ### Load topofile
  (topo, geoTrans) = gdalLoad(topoFile)
  
  
  ### Convert point locations to raster locations
  #Easting to Column
  pts[:,0] = (pts[:,0]-geoTrans[0])/geoTrans[1]
  #Northing to Row
  pts[:,1] = (pts[:,1]-geoTrans[3])/geoTrans[5]
  
  outputPts = np.zeros((len(pts),3))
  
  ### Find all maxima in topo map
  Maxes = getMaxima(topo, peakNeighborhood, peakMinHeight)
  
  if len(Maxes) > 0:
    Peaks = np.zeros((len(Maxes),5))
    Peaks[:,0], Peaks[:,1] = Maxes[:,0], Maxes[:,1]
  else:
    sys.stderr.write('No Peaks Found\n\n')
    return False
  
  
  ### Find peaks nearby points, keep track of closest point
  # Calculate the distance matrix. elements list distances from each peak to all pts
  sys.stdout.write('Creating the distance matrix...\n')
  distMatrix = distance_matrix(Peaks[:,0:2],pts[:,0:2])
  
  # Cut all peaks farther than the radius away from any point
  Peaks = Peaks[np.where(np.min(distMatrix,axis=1) <= peakMaxRadius)]
  distMatrix = distMatrix[np.where(np.min(distMatrix,axis=1) <= peakMaxRadius)]
  if len(Peaks)<1:
    sys.stderr.write('No peaks near vents!! Check that input DEM and peak locs'+
                            ' are coregistered!\n\n')
    return None
  
  # Assign closest pt index to Peaks
  Peaks[:,4] = np.argmin(distMatrix,axis=1)
  sys.stdout.write('  %d peaks are within %0.2f units of a vent\n' % (len(Peaks),
                            peakMaxRadius))
  
  ### Find candidate peaks for some points that don't have their own cand. peak
  nopeaks = filter(lambda x: x not in set(Peaks[:,4]), range(len(pts)))
  
  #find closest peak to each vent
  ventMinDist = np.min(distMatrix[:,nopeaks],axis=0)
  for i, emptyvent in enumerate(nopeaks):
    # If there is at least one close peak
    if ventMinDist[i] < peakMaxRadius:
      newPeak = np.argmin(distMatrix[:,emptyvent])
      #If the closest peaks' corresponding vent has more than 1 cand. peak
      #Then it's ok to "steal" it:
      if len(np.where(Peaks[:,4]==Peaks[newPeak,4])[0]) > 1:
        Peaks[newPeak,4] = emptyvent

  nopeaks = filter(lambda x: x not in set(Peaks[:,4]), range(len(pts)))
  
  sys.stdout.write('  %d vents have no candidate peaks\n' % len(nopeaks))
  
  ### Get Prominence for ALL candidate peaks
  Peaks = getAllProms(Peaks,topo)
  
  ### Find peak with highest prominence for each vent
  for i, pt in enumerate(outputPts):
    if pt[2] < 1:
      candPeaks = Peaks[np.where(Peaks[:,4] == i)]
      if len(candPeaks) >= 1:
        thisProm = np.max(candPeaks[:,2])
        realPeak = candPeaks[np.where(candPeaks[:,2] == thisProm)[0][0]]
        # Get peak position in lat long from peak position column, row
        try:
          pt[0] = realPeak[:,0]*geoTrans[1] + geoTrans[0]
        except IndexError:
          pt[0] = realPeak[0]*geoTrans[1] + geoTrans[0]
        try:
          pt[1] = realPeak[:,1]*geoTrans[5] + geoTrans[3]
        except IndexError:
          pt[1] = realPeak[1]*geoTrans[5] + geoTrans[3]
        pt[2] = thisProm
  
  return outputPts


def getAllProms(Peaks,topo):
  '''
  Calculate Prominence of many points using findProminence
  '''
  import numpy as np
  
  #Make a list for all prominences
  proms = np.zeros(len(Peaks))
  
  #Find prominence of ALL Candidate Peaks
  for i, pk in enumerate(Peaks):
    #find grid cell of peak
    peakCell = [int(pk[1]),int(pk[0])]
    peakElev = topo[peakCell[0],peakCell[1]]
    
    topoShape = topo.shape
    #find a neighborhood where there's at least one higher point
    n=50
    minR = max(peakCell[0]-n, 0)
    maxR = min(peakCell[0]+n+1,topoShape[0]-1)
    minC = max(peakCell[1]-n, 0)
    maxC = min(peakCell[1]+n+1,topoShape[1]-1)
    
    while np.max(topo[minR:maxR,minC:maxC])==peakElev:
      n += 100
      minR = max(peakCell[0]-n, 0)
      maxR = min(peakCell[0]+n+1,topoShape[0]-1)
      minC = max(peakCell[1]-n, 0)
      maxC = min(peakCell[1]+n+1,topoShape[1]-1)
      
      
    #add 150 cells to neighborhood for better chance of finding real key col  
    n += 150
    minR = max(peakCell[0]-n, 0)
    maxR = min(peakCell[0]+n+1,topoShape[0]-1)
    minC = max(peakCell[1]-n, 0)
    maxC = min(peakCell[1]+n+1,topoShape[1]-1)
    
    #calculate topographic prominence in the vent's neighborhood
    proms[i] = findProminence(topo[minR:maxR,minC:maxC],peakElev)
    
    if (i+1)%10==0:
      print ("   Found Prominence of %-6d/%d Peaks" % ((i+1),len(Peaks)))
  return proms
  


def getMaxima(topo, neighborhood_size, vertical_threshold):
  '''
  Find maxima of a raster
  '''
  import numpy as np
  import scipy.ndimage as ndimage
  import scipy.ndimage.filters as filters
  
  #find potential maxima from a maximum filter
  data_max = filters.maximum_filter(topo, neighborhood_size)
  maxima = (topo == data_max)
  print ("  Found %d Potential Peaks" % len(np.where(maxima)[0]))
  print ("Removing low lying peaks")
  #remove maxima that are too small
  data_min = filters.minimum_filter(topo, neighborhood_size)
  diff = ((data_max - data_min) > vertical_threshold)
  maxima[diff == 0] = 0
  
  labeled, num_maxima = ndimage.label(maxima)
  print ("  Confirmed %d Real Peaks" % num_maxima)
  
  x, y = [], []
  if num_maxima > 0:
      slices = ndimage.find_objects(labeled)
      for dy,dx in slices:
          x_center = (dx.start + dx.stop - 1)/2
          x.append(x_center)
          y_center = (dy.start + dy.stop - 1)/2
          y.append(y_center)
  
      maxes = np.zeros((num_maxima,2)) - 1
      maxes[:,0], maxes[:,1] = x, y
  return maxes









