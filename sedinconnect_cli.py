#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  8 13:45:31 2022

@author: François Clapuyt
"""

import sys, os
from time import time, sleep
import time
import osgeo
from osgeo import gdal, ogr
import numpy
import math
import struct
import shutil
from osgeo.gdal import *
import threading
import gdal_rasterize as poly2ras
import base64
import scipy.signal
from scipy.special import _ufuncs_cxx #this for dependencies, otherwise pyinstaller does not find _ufuncs_cxx trying to create scipy dll
import matplotlib.pyplot

"""
--------------------------------------
###########################################################################       

COMPUTE SINKS

###########################################################################
--------------------------------------
"""       
# Compute sinks
def sinks(dtm_f, sink_shp):
    #
    filename = dtm_f.replace('\\','/')
    filename = str(filename)
    tif = osgeo.gdal.Open(filename)#opening the file (1 band, bit not specified, but usually 32)
    #
    #folder path
    dir_pat=os.path.dirname(os.path.realpath(filename))#path of the selected input file
    dir_path=dir_pat.encode('ascii','ignore')
    dir_path=dir_path.replace('\\','/')
    #
    # defining a global function to use for passing the filenames created inside the functions
    #without making the user specify each filepath/filename
    global path_sink_dtm 
    path_sink_dtm = dir_path
    #
    # Opening Messages
    if tif is None:
      print(time.strftime("%d/%m/%Y %H:%M:%"    ),"couldn't open input DTM dataset")
      sys.exit(1)
    else:
      print(time.strftime("%d/%m/%Y %H:%M:%S    "),"opening DTM for sinks was successful!")
    #
    #Colonne, righe, bande
    cols = tif.RasterXSize
    rows = tif.RasterYSize
    bands =  tif.RasterCount
    #
    #dtm as an array                #################################       ATTENZIONE     ##########################
    tif_ar = tif.ReadAsArray() #con GDAL 1.6 e 1.7 non funzia 90%, con gdal bindings per windows 1.9.2 sembra funzionare provare anche a metter su ultime versioni di numpy
    tif_ar = tif_ar.astype(float)
    band = tif.GetRasterBand(1)#bands retrieving
    geoinf = tif.GetGeoTransform()#upper left coordinates, cellsize
    proj = tif.GetProjection()
    tipo = tif_ar.dtype
    del(tif)
    #
    file_sk = sink_shp.replace('\\','/')
    file_sk =str(file_sk)
    
    source_sk = ogr.Open(file_sk, 1)
    layer_sk = source_sk.GetLayer()
    layer_defn_sk= layer_sk.GetLayerDefn()
    field_names_sk = [layer_defn_sk.GetFieldDefn(i).GetName() for i in range(layer_defn_sk.GetFieldCount())]

    try:
        in_sk_field = field_names_sk.index('sink_id')
        msg_sk = ('sink_id field found in field n° ' +  str(in_sk_field + 1) + '!')
        print (time.strftime("%d/%m/%Y %H:%M:%S    "), msg_sk)
    except:
        msg_sk = ('sink_id field not found, ignoring input')
        print (time.strftime("%d/%m/%Y %H:%M:%S    "), msg_sk)

    #
    sk=poly2ras.rasterize(file_sk, dir_path+"/sinks.tif", prototype=filename, options=["ATTRIBUTE=sink_id"])#assign value to variable or arbitrary one
    sink_ar=sk.GetRasterBand(1).ReadAsArray()
    del (sk)
    #
    #
    #D8 flow directions & outputs to derive watersheds
    os.system((("mpiexec -n 8 D8Flowdir -p ").lower())+filename[0:-4]+"p.tif"+ " -sd8 "+filename[0:-4]+"sd8.tif -fel "+filename)#unix case
    #
    #reclassifying D8 flow directions
    tif_fdir8=osgeo.gdal.Open(filename[0:-4]+"p.tif")#opening the flowdir file (1 band, bit not specified, but usually 32)
    if tif_fdir8 is None:
      print(time.strftime("%d/%m/%Y %H:%M:%S    "),"couldn't open input dir8 dataset")
      sys.exit(1)
    else:
      print(time.strftime("%d/%m/%Y %H:%M:%S    "),"Flowdir opening was successful!")
    #Array 
    tif_fdir8_ar=tif_fdir8.ReadAsArray()
    tif_fdir8_ar=tif_fdir8_ar.astype(numpy.float)
    ndv=numpy.min(tif_fdir8_ar)
    tif_fdir8_ar[tif_fdir8_ar==ndv]=-1#to Basin Code
    del(tif_fdir8)
    os.remove(filename[0:-4]+"p.tif")
    os.remove(dir_path+"/sinks.tif")
    #
    tif_sd8=osgeo.gdal.Open(filename[0:-4]+"sd8.tif")#opening the flowdir file (1 band, bit not specified, but usually 32)
    if tif_sd8 is None:
      print(time.strftime("%d/%m/%Y %H:%M:%S    "),"couldn't open input sd8 dataset")
      sys.exit(1)
    else:
      print(time.strftime("%d/%m/%Y %H:%M:%S    "),"Sd8 opening was successful!")
    #Array 
    tif_sd8_ar=tif_sd8.ReadAsArray()
    tif_sd8_ar=tif_sd8_ar.astype(numpy.float)
    ndv=numpy.min(tif_sd8_ar)
    tif_sd8_ar[tif_sd8_ar==ndv]=-1#to Basin Code
    del(tif_sd8)
    #
    #removing files
    os.remove(filename[0:-4]+"sd8.tif")#we need to sink this file too and save as global, it will then be used in the connectivity with sinks
    #
    #Dinf_flow directions I will then use in the connectivity with sinks to avoid double runs of flowdir  that delete cells and create boundary watersheds
    os.system((("mpiexec -n 8 DinfFlowdir -ang ").lower())+filename[0:-4]+"ang.tif"+ " -slp "+filename[0:-4]+"slp.tif -fel "+filename)#unix case
    #removing unecessary Dinf slope tif file
    os.remove(filename[0:-4]+"slp.tif")
    #
    tif_dirinf=osgeo.gdal.Open(filename[0:-4]+"ang.tif")#opening the flowdir file (1 band, bit not specified, but usually 32)
    if tif_dirinf is None:
      print(time.strftime("%d/%m/%Y %H:%M:%S    "),"couldn't open input dirinf dataset")
      sys.exit(1)
    else:
      print(time.strftime("%d/%m/%Y %H:%M:%S    "),"Dirinf opening was successful!")
    #Array 
    tif_dirinf_ar=tif_dirinf.ReadAsArray()
    tif_dirinf_ar=tif_dirinf_ar.astype(numpy.float)
    ndv=numpy.min(tif_dirinf_ar)
    tif_dirinf_ar[tif_dirinf_ar==ndv]=-1#to Basin Code
    del(tif_dirinf)
    #
    #
    #bigger matrix to move with indexes, will store the propagated basin code
    big_wat_sink = numpy.zeros(shape=((tif_ar.shape[0])+2,(tif_ar.shape[1])+2),dtype=numpy.float32)#add rows and columns of zeros to elaborate better matrix indexes
    big_wat_sink[1:-1,1:-1] = sink_ar
    #
    #zero matrix bigger than F_dir8, to avoid border indexing problems
    Fd8 = numpy.zeros(shape=((tif_fdir8_ar.shape[0])+2,(tif_fdir8_ar.shape[1])+2),dtype=numpy.float32)#add rows and columns of zeros to elaborate better matrix indexes
    Fd8[1:(Fd8.shape[0]-1),1:(Fd8.shape[1]-1)] = Fd8[1:(Fd8.shape[0]-1),1:(Fd8.shape[1]-1)]+tif_fdir8_ar
    Fdir8 = Fd8
    #
    #Creating a bigger matrix as large as weight(and all the matrices) to store the weighted flow length values
    SK_m = numpy.zeros(shape=((tif_ar.shape[0])+2,(tif_ar.shape[1])+2),dtype=numpy.float32)
    #
    #
    #DOING THE PROCEDURE
    #Let's go for the search and algo-rhytm for the weighted-Flow-Length
    #
    start = time.time.time()#take your time
    #
    print (msg_sk)
    #
    print (time.strftime("%d/%m/%Y %H:%M:%S    "),'Computing sinks watershed extraction ... ...')
    #I retrieve the basin code theat I will propagate
    SV = numpy.where(big_wat_sink > 0)#fast coordinates all the Basin Code values, starting from them to go forward and compute flow length
    #
    Y = SV[0]#rows, code indexes
    X = SV[1]#columns, code indexes pay attention not to invert values !!!!!!!!!!!!!!
    #
    #initializing lists for outlet and moving cell coordinates, in function of their position
    YC1=[];YC2=[];YC3=[];YC4=[];YC5=[];YC6=[];YC7=[];YC8=[]
    XC1=[];XC2=[];XC3=[];XC4=[];XC5=[];XC6=[];XC7=[];XC8=[]
    #
    #   Flow Directions Taudem
    #   4   3   2
    #   5   -   1
    #   6   7   8
    #
    #   Draining in Direction Matrix
    #   8   7   6
    #   1   -   5
    #   2   3   4
    #
    i1=Fdir8[Y,X-1]#Searching for Basin Code with cells draining into them, 8 directions
    D1=numpy.where(i1==1)#
    YC1.extend(Y[D1])#coordinates satisfacting the conditions
    XC1.extend(X[D1])
    SK_m[YC1,XC1]=big_wat_sink[YC1,XC1]#initialize basin code at cells draining to Basin Code
    #
    i2=Fdir8[Y+1,X-1]#Searching for Basin Code with cells draining into them, 8 directions
    D2=numpy.where(i2==2)#
    YC2.extend(Y[D2])#coordinates satisfacting the conditions
    XC2.extend(X[D2])
    SK_m[YC2,XC2]=big_wat_sink[YC2,XC2]#initialize basin code at cells draining to Basin Code
    #
    i3=Fdir8[Y+1,X]#Searching for Basin Code with cells draining into them, 8 directions
    D3=numpy.where(i3==3)#
    YC3.extend(Y[D3])#coordinates satisfacting the conditions
    XC3.extend(X[D3])
    SK_m[YC3,XC3]=big_wat_sink[YC3,XC3]#initialize basin code at cells draining to Basin Code
    #
    i4=Fdir8[Y+1,X+1]#Searching for Basin Code with cells draining into them, 8 directions
    D4=numpy.where(i4==4)#
    YC4.extend(Y[D4])#coordinates satisfacting the conditions
    XC4.extend(X[D4])
    SK_m[YC4,XC4]=big_wat_sink[YC4,XC4]#initialize basin code at cells draining to Basin Code
    #
    i5=Fdir8[Y,X+1]#Searching for Basin Code with cells draining into them, 8 directions
    D5=numpy.where(i5==5)#
    YC5.extend(Y[D5])#coordinates satisfacting the conditions
    XC5.extend(X[D5])
    SK_m[YC5,XC5]=big_wat_sink[YC5,XC5]#initialize basin code at cells draining to Basin Code
    #
    i6=Fdir8[Y-1,X+1]#Searching for Basin Code with cells draining into them, 8 directions
    D6=numpy.where(i6==6)#
    YC6.extend(Y[D6])#coordinates satisfacting the conditions
    XC6.extend(X[D6])
    SK_m[YC6,XC6]=big_wat_sink[YC6,XC6]#initialize basin code at cells draining to Basin Code
    #
    i7=Fdir8[Y-1,X]#Searching for Basin Code with cells draining into them, 8 directions
    D7=numpy.where(i7==7)#
    YC7.extend(Y[D7])#coordinates satisfacting the conditions
    XC7.extend(X[D7])
    SK_m[YC7,XC7]=big_wat_sink[YC7,XC7]#initialize basin code at cells draining to Basin Code
    #
    i8=Fdir8[Y-1,X-1]#Searching for Basin Code with cells draining into them, 8 directions
    D8=numpy.where(i8==8)#
    YC8.extend(Y[D8])#coordinates satisfacting the conditions
    XC8.extend(X[D8])
    SK_m[YC8,XC8]=big_wat_sink[YC8,XC8]#initialize basin code at cells draining to Basin Code
    #
    #start =time.time()#da cancellare poi.....!!!!!! Solo per check
    count=1# "0" passage already done during the previous step
    while len (YC1) or len(YC2) or len(YC3) or len(YC4) or len(YC5) or len(YC6) or len(YC7) or len(YC8) >0:
        #Converting into array to be able to do operations
        YYC1=numpy.asarray(YC1);XXC1=numpy.asarray(XC1)
        YYC2=numpy.asarray(YC2);XXC2=numpy.asarray(XC2)
        YYC3=numpy.asarray(YC3);XXC3=numpy.asarray(XC3)
        YYC4=numpy.asarray(YC4);XXC4=numpy.asarray(XC4)
        YYC5=numpy.asarray(YC5);XXC5=numpy.asarray(XC5)
        YYC6=numpy.asarray(YC6);XXC6=numpy.asarray(XC6)
        YYC7=numpy.asarray(YC7);XXC7=numpy.asarray(XC7)
        YYC8=numpy.asarray(YC8);XXC8=numpy.asarray(XC8)
        #
        #Now I will propagate always the same basin code
        #
        YYC1=(YYC1);XXC1=(XXC1-1)#
        YYC2=(YYC2+1);XXC2=(XXC2-1)#
        YYC3=(YYC3+1);XXC3=(XXC3)#l
        YYC4=(YYC4+1);XXC4=(XXC4+1)#
        YYC5=(YYC5);XXC5=(XXC5+1)#
        YYC6=(YYC6-1);XXC6=(XXC6+1)#
        YYC7=(YYC7-1);XXC7=(XXC7)#
        YYC8=(YYC8-1);XXC8=(XXC8-1)#
        #  
        if len(YYC1) > 0:
            SK_m[YYC1,XXC1]=SK_m[YC1,XC1]
        else:
            pass
        if len(YYC2) > 0:
            SK_m[YYC2,XXC2]=SK_m[YC2,XC2]
        else:
            pass
        if len(YYC3) > 0:
            SK_m[YYC3,XXC3]=SK_m[YC3,XC3]    
        else:
            pass
        if len(YYC4) > 0:
            SK_m[YYC4,XXC4]=SK_m[YC4,XC4]
        else:
            pass
        if len(YYC5) > 0:
            SK_m[YYC5,XXC5]=SK_m[YC5,XC5]    
        else:
            pass
        if len(YYC6) > 0:
            SK_m[YYC6,XXC6]=SK_m[YC6,XC6]   
        else:
            pass
        if len(YYC7) > 0:
            SK_m[YYC7,XXC7]=SK_m[YC7,XC7]   
        else:
            pass
        if len(YYC8) > 0:
            SK_m[YYC8,XXC8]=SK_m[YC8,XC8]
        else:
            pass
        #
        #Reconstructing all X and Y of this step and moving on upwards (Downstream if you think in GIS, right?)
        YY=[];XX=[]
        YY.extend(YYC1);XX.extend(XXC1)
        YY.extend(YYC2);XX.extend(XXC2)
        YY.extend(YYC3);XX.extend(XXC3)
        YY.extend(YYC4);XX.extend(XXC4)
        YY.extend(YYC5);XX.extend(XXC5)
        YY.extend(YYC6);XX.extend(XXC6)
        YY.extend(YYC7);XX.extend(XXC7)
        YY.extend(YYC8);XX.extend(XXC8)
        #
        YY=numpy.asarray(YY)
        XX=numpy.asarray(XX)
        #
        i1=Fdir8[YY,XX-1]#Searching for cells draining into them, 8 directions
        D1=numpy.where(i1==1)#
        YC1=YY[D1]#coordinates satisfacting the conditions
        XC1=XX[D1]
        #
        i2=Fdir8[YY+1,XX-1]#Searching for cells draining into them, 8 directions
        D2=numpy.where(i2==2)#
        YC2=YY[D2]#coordinates satisfacting the conditions
        XC2=XX[D2]
        #
        i3=Fdir8[YY+1,XX]#Searching for cells draining into them, 8 directions
        D3=numpy.where(i3==3)#
        YC3=YY[D3]#coordinates satisfacting the conditions
        XC3=XX[D3]
        #
        i4=Fdir8[YY+1,XX+1]#Searching for cells draining into them, 8 directions
        D4=numpy.where(i4==4)#
        YC4=YY[D4]#coordinates satisfacting the conditions
        XC4=XX[D4]
        #
        i5=Fdir8[YY,XX+1]#Searching for cells draining into them, 8 directions
        D5=numpy.where(i5==5)#
        YC5=YY[D5]#coordinates satisfacting the conditions
        XC5=XX[D5]
        #
        i6=Fdir8[YY-1,XX+1]#Searching for cells draining into them, 8 directions
        D6=numpy.where(i6==6)#
        YC6=YY[D6]#coordinates satisfacting the conditions
        XC6=XX[D6]
        #
        i7=Fdir8[YY-1,XX]#Searching for cells draining into them, 8 directions
        D7=numpy.where(i7==7)#
        YC7=YY[D7]#coordinates satisfacting the conditions
        XC7=XX[D7]
        #
        i8=Fdir8[YY-1,XX-1]#Searching for cells draining into them, 8 directions
        D8=numpy.where(i8==8)#
        YC8=YY[D8]#coordinates satisfacting the conditions
        XC8=XX[D8]
        count=count+1
    #
    #
    elapsed=(time.time()-start)#computational time
    print(time.strftime("%d/%m/%Y %H:%M:%S    "),"Process concluded succesfully \n","%.2f" % elapsed,'seconds for sink watersheds calculation with ',int(count),' iterations') #truncating the precision
    #os.system("pause")
    #
    del Fdir8
    #
    SK_m= SK_m[1:SK_m.shape[0]-1,1:SK_m.shape[1]-1] #reshaping weigthed flow length, we need this step to homogenize matrices dimensions!!!!!!!!!!
    #I put NoData in DTM where sink watershed >0
    min_tif_ar = tif_ar.min()
    tif_ar[numpy.where(SK_m>0)] = min_tif_ar
    tif_ar[numpy.where(tif_ar == min_tif_ar)] = -9999
    # I need to "cut" also the flow direction array and to write it to a file, it will be used then into the Connectivity
    #otherwise another run of flow direction would cut a round of pixels creating strange interbasins
    min_fdir8_array =tif_fdir8_ar.min()#valore minimo era ex no-data
    tif_fdir8_ar[numpy.where(SK_m>0)] = min_fdir8_array
    tif_fdir8_ar[numpy.where(tif_fdir8_ar == min_fdir8_array)] = -9999
    #
    #same procedure for slope D8 and later on for Dinf directions
    min_sd8_array =tif_sd8_ar.min()#valore minimo era ex no-data
    tif_sd8_ar[numpy.where(SK_m>0)] = min_sd8_array
    tif_sd8_ar[numpy.where(tif_sd8_ar == min_sd8_array)] = -9999
    #
    #same procedure for Dinf directions
    min_dirinf_array =tif_dirinf_ar.min()#valore minimo era ex no-data
    tif_dirinf_ar[numpy.where(SK_m>0)] = min_dirinf_array
    tif_dirinf_ar[numpy.where(tif_dirinf_ar == min_dirinf_array)] = -9999
    #
    #
    #write sinked DTM
    SK_m_ds = gdal.GetDriverByName('GTiff').Create((dir_path+"/sinked_dtm.tif"), tif_ar.shape[1], tif_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    SK_m_ds.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    SK_m_ds.SetProjection(proj)
    SK_m_ds.GetRasterBand(1).SetNoDataValue(-9999) 
    SK_m_ds.GetRasterBand(1).WriteArray(tif_ar,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
    SK_m_ds=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del SK_m_ds
    #
    #write sinked Flow Dir
    SK_fdir8 = gdal.GetDriverByName('GTiff').Create((dir_path+"/sinked_fdir8.tif"), tif_ar.shape[1], tif_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    SK_fdir8.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    SK_fdir8.SetProjection(proj)
    SK_fdir8.GetRasterBand(1).SetNoDataValue(-9999) 
    SK_fdir8.GetRasterBand(1).WriteArray(tif_fdir8_ar,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
    SK_fdir8=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del SK_fdir8 
    #
    #write sd8 array
    #write sinked Flow Dir
    SK_sd8 = gdal.GetDriverByName('GTiff').Create((dir_path+"/sinked_sd8.tif"), tif_ar.shape[1], tif_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    SK_sd8.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    SK_sd8.SetProjection(proj)
    SK_sd8.GetRasterBand(1).SetNoDataValue(-9999) 
    SK_sd8.GetRasterBand(1).WriteArray(tif_sd8_ar,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
    SK_sd8=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del SK_sd8
    #
    #write sinked Flow Dirinf
    SK_dirinf= gdal.GetDriverByName('GTiff').Create((dir_path+"/sinked_dirinf.tif"), tif_ar.shape[1], tif_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    SK_dirinf.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    SK_dirinf.SetProjection(proj)
    SK_dirinf.GetRasterBand(1).SetNoDataValue(-9999) 
    SK_dirinf.GetRasterBand(1).WriteArray(tif_dirinf_ar,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
    SK_dirinf=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del SK_dirinf 
    
    os.remove(filename[0:-4]+"ang.tif")
    
    print(time.strftime("%d/%m/%Y %H:%M:%"),'Sinks computation concluded... great!')
    
    
"""
--------------------------------------
###########################################################################        

CAVALLI ROUGHNESS FUNCTION

###########################################################################
--------------------------------------
"""
# Roughness index calculation
def rw_cavalli(dtm_r, f_s, ri_out, w_out, sk_flag = 0):
    filename = dtm_r.replace('\\','/')
    filename = str(filename)
    tif = osgeo.gdal.Open(filename)#opening the file (1 band, bit not specified, but usually 32)
    #
    #folder path
    dir_pat = os.path.dirname(os.path.realpath(filename))#path of the selected input file
    dir_path = dir_pat.encode('ascii','ignore')
    dir_path = dir_path.replace('\\','/')
    #
    # Opening Messages
    if tif is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't open input DTM dataset, for Weighting factor computatin")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "opening DTM for weighting factor computation was successful!")
    #
    #Colonne, righe, bande
    cols = tif.RasterXSize
    rows = tif.RasterYSize
    bands =  tif.RasterCount
    #
    #dtm as an array                #################################       ATTENZIONE     ##########################
    I_ar = tif.ReadAsArray() #con GDAL 1.6 e 1.7 non funzia 90%, con gdal bindings per windows 1.9.2 sembra funzionare provare anche a metter su ultime versioni di numpy
    I_ar = I_ar.astype(float)
    band = tif.GetRasterBand(1)#bands retrieving
    geoinf = tif.GetGeoTransform()#upper left coordinates, cellsize
    proj = tif.GetProjection()
    tipo = I_ar.dtype
    del(tif)
    #
    #mask_nodata
    ndv = numpy.min(I_ar)
    I_ar[I_ar == ndv] = numpy.nan
    #
    #kernel size
    size_filter = numpy.float32(int(f_s))#importante che resti float32, per omogeneita'
    #
    start =time.time()#take the time
    #
    I_ar_p = numpy.ones((rows,cols), dtype=tipo) #matrice di numeri uno e zero sui nodata da usare per far la media poi
    I_ar_p[numpy.isnan(I_ar)==1] = 0
    I_ar_d = numpy.ones((rows,cols), dtype=tipo) #build up array for summing correctly elevation data later
    I_ar_d[numpy.isnan(I_ar)==1] = 0
    I_ar_d[numpy.isnan(I_ar)==0] = I_ar[numpy.isnan(I_ar)==0]
    #
    #
    ker = numpy.ones((size_filter,size_filter))#creating weights for sum and then an averaging filter
    #C = scipy.signal.convolve2d(I_ar,ker, 'same')#borders kept the same but not realistic values, it divides by n^2 even if it finds only 2 elements scipy version
    D_num_el = scipy.signal.convolve2d(I_ar_p, ker, 'same')#i get the number of elements per window
    D_num_el[numpy.isnan(I_ar)==1] = numpy.nan
    
    E_sum = scipy.signal.convolve2d(I_ar_d,ker, 'same') #sum of values in the DTM
    E_sum[numpy.isnan(I_ar)==1] = numpy.nan
    
    M = E_sum/D_num_el #average =sum/number of elements
    DTM_R  = numpy.float64(I_ar) - numpy.float64(M) #residual DTM, now I need to compute the Std dev of this averaged DTM (so I need again mean values of this and then squared values to fasten the computation)
    #
    """Use float64 to have more precision in decimals"""
    #
    DTM_R[numpy.isnan(I_ar)==1] = 0
    # mean of residual topography
    E_R = scipy.signal.convolve2d(DTM_R, ker, 'same') #sum of values in the DTM
    E_R[numpy.isnan(I_ar)==1] = numpy.nan
    M_R = E_R/D_num_el #average =sum/number of elements
    #Now playing for standard deviation following: http://stackoverflow.com/questions/18419871/improving-code-efficiency-standard-deviation-on-sliding-windows
    #
    #DTM_SQ = numpy.float64(I_ar_d) * numpy.float64(I_ar_d) #build up array for summing correctly elevation data later, this way I squared the values otherwise I loose the decimal precision if array is starting with 4 digits+ decimal + sign
    #
    DTM_R_SQ = numpy.float64(DTM_R) *  numpy.float64(DTM_R) #residuals squared
    #
    M_R_SQ = numpy.power(M_R, 2)
    E_R_SQ = scipy.signal.convolve2d(DTM_R_SQ, ker, 'same') #sum of values in the DTM
    E_R_SQ[numpy.isnan(I_ar)==1] = numpy.nan
    #
    M_R_SQ = E_R_SQ/D_num_el# E[X^2] average of squared residual DTM values
    MRSQ = numpy.power(M_R, 2) # square of average residual DTM values
    Ri = numpy.sqrt(M_R_SQ - MRSQ) # E[X^2] - E[X]^2       WC = Weight Cavalli
    #
    """Pay attention to the number we're treating, if we need to compute std dev via convolve and squared trick, for a DTM we need to use numpy.float64(DTM), otherwise the digits are not enough"""
            
    del (I_ar_p, I_ar_d, D_num_el, ker, E_sum, M, DTM_R, E_R, M_R, DTM_R_SQ, M_R_SQ, E_R_SQ, MRSQ)
    # 
    #
    #Check flag for use of sinks and sink the wighting factor accordingly
    if sk_flag == 0:
        pass
    else:
        sk_tmp = path_sink_dtm+"/sinked_dtm.tif"
        sk_tmp = str(sk_tmp)
        sk_dtm = osgeo.gdal.Open(sk_tmp)#opening the flowdir file (1 band, bit not specified, but usually 32)
        if sk_dtm is None:
          print (time.strftime("%d/%m/%Y %H:%M:%S    "),"couldn't open input sinked_dtm dataset")
          sys.exit(1)
        else:
          print (time.strftime("%d/%m/%Y %H:%M:%S    "),"sinked_dtm opening was successful!")
        #Array 
        sk_dtm_ar=sk_dtm.ReadAsArray()
        sk_dtm_ar=sk_dtm_ar.astype(numpy.float)
        ndv=numpy.min(sk_dtm_ar)
        sk_dtm_ar[sk_dtm_ar==ndv] = numpy.NaN #to Basin Code
        del(sk_dtm)
        Ri[sk_dtm_ar == numpy.NaN] = numpy.NaN
    #
    #
    sur_rough = gdal.GetDriverByName('GTiff').Create(str(ri_out), I_ar.shape[1], I_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    sur_rough.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    sur_rough.SetProjection(proj)
    sur_rough.GetRasterBand(1).WriteArray(Ri,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
    sur_rough.GetRasterBand(1).GetStatistics(0,1) #calculate statistics for visualization
    sur_rough = None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del sur_rough
    #mean_t[y,x]=numpy.average(isnan([I_ar_t[y,x], I_ar_t[y-1,x], I_ar_t[y+1,x], I_ar_t[y,x-1], I_ar_t[y,x+1], I_ar_t[y-1,x-1], I_ar_t[y+1,x+1], I_ar_t[y+1,x-1], I_ar_t[y-1,x+1]])==False)
    #
    ##compute the weighting factor as 1-(R/(R+0.001))
    print (time.strftime("%d/%m/%Y %H:%M:%S    "), 'computing weighting factor ... ...')
    #
    maxr = float(Ri[numpy.isnan(Ri) == False].max())#max element in surface roughness
    #Mr = maxr + 0.001#increase a bit, in this release, old style of standardization
    if normseba.isChecked():
        #
        # minima in Roughness set to 0.001
        Ri[numpy.where(Ri < 0.001)] = 0.001
        # detect Roughness minimum or pick a 0.001 minimum
        minr = float(Ri[numpy.isnan(Ri) == False].min())# min element in surface roughness
        #
        weig_fac = 1.0 - ((numpy.log(Ri) - numpy.log(minr)) / (numpy.log(maxr) - numpy.log(minr)))
        weig_fac[numpy.where(weig_fac < 0.001)] = 0.001 #all avalues in the range 0-0-001 are shifted to 0.001, the rest remains NoData
    else:
        weig_fac = 1.0 - (Ri/maxr)#1-... to avoid 1-1=0 at denominator in downslope component
        weig_fac[numpy.where(weig_fac < 0.001)] = 0.001 #all avalues in the range 0-0-001 are shifted to 0.001, the rest remains NoData
    #
    print (time.strftime("%d/%m/%Y %H:%M:%S    "), 'Weighting factor calculated!')
    #
    wf = gdal.GetDriverByName('GTiff').Create(str(w_out), I_ar.shape[1], I_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    wf.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    wf.SetProjection(proj)
    wf.GetRasterBand(1).WriteArray(weig_fac,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
    wf.GetRasterBand(1).GetStatistics(0,1) #calculate statistics for visualization
    wf = None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del wf
    del I_ar
    del Ri
    del weig_fac
    print (time.strftime("%d/%m/%Y %H:%M:%S    "),'RI and W files saved!')

"""
--------------------------------------
###########################################################################       

CONNECTIVITY TO THE OUTLET

###########################################################################
--------------------------------------
"""        
def CavalliConnectivityout (dtm_f, c_s, w_f, out_ic, sink_flag = 0):#I will flag sink_are to 1 in case I need to retrieve the fdir_8 from the sink analysis    
    
    filename=dtm_f.replace('\\','/')
    filename=str(filename)
    tif=osgeo.gdal.Open(filename)
    #folder path
    dir_path=os.path.dirname(os.path.realpath(filename))#path of the selected input file
    dir_path=dir_path.replace('\\','/')
    #
    # Opening Messages
    if tif is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't open input dataset (connectivity main DTM)")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "opening DTM for Connectivity was successful!")
    #
    #Colonne, righe, bande
    cols = tif.RasterXSize
    rows = tif.RasterYSize
    bands =  tif.RasterCount
    #
    #dtm as an array                #################################       ATTENZIONE     ##########################
    tif_ar=tif.ReadAsArray() #con GDAL 1.6 e 1.7 non funzia 90%, con gdal bindings per windows 1.9.2 sembra funzionare provare anche a metter su ultime versioni di numpy
    tif_ar = tif_ar.astype(float)
    band=tif.GetRasterBand(1)#bands retrieving
    geoinf = tif.GetGeoTransform()#upper left coordinates, cellsize
    proj=tif.GetProjection()
    tipo=tif_ar.dtype
    del(tif)
    #
    #
    ndv=numpy.min(tif_ar)
    tif_ar[tif_ar==ndv]=numpy.NaN
    #cell_size
    cell_s=numpy.float32(str(c_s))#importante che resti float32, per omogeneita'
    #
    #create constant array to trasform into raster
    const_ar=tif_ar*0+cell_s#array cellsize
    #
    if sink_flag == 0:
        #D8 flow directions & outputs
        #os.system("mpiexec -n 8 D8Flowdir -p "+filename[0:-4]+"p.tif"+ " -sd8 "+filename[0:-4]+"sd8.tif -fel "+filename)#windows case
        os.system((("mpiexec -n 8 D8Flowdir -p ").lower())+filename[0:-4]+"p.tif"+ " -sd8 "+filename[0:-4]+"sd8.tif -fel "+filename)#unix case
        """
        ---------------------------------------------
        WORKING ON D_DOWN COMPONENT
        ---------------------------------------------
        """    
        tif_fdir8=osgeo.gdal.Open(filename[0:-4]+"p.tif")#opening the file (1 band, bit not specified, but usually 32)
        if tif_fdir8 is None:
          print (time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't open input fd8 dataset")
          sys.exit(1)
        else:
          print (time.strftime("%d/%m/%Y %H:%M:%S    "), "opening fd8 was successful!")
        #Array 
        tif_fdir8_ar=tif_fdir8.ReadAsArray()
        tif_fdir8_ar=tif_fdir8_ar.astype(numpy.float)#otherwise overflow in future operations
        ndv=numpy.min(tif_fdir8_ar)
        tif_fdir8_ar[tif_fdir8_ar==ndv]=numpy.NaN
        del(tif_fdir8)
        #
        #
        #Dinf flow directions & outputs
        #os.system("mpiexec -n 8 DinfFlowdir -ang "+filename[0:-4]+"ang.tif"+ " -slp "+filename[0:-4]+"slp.tif -fel "+filename)#windows case
        os.system((("mpiexec -n 8 DinfFlowdir -ang ").lower())+filename[0:-4]+"ang.tif"+ " -slp "+filename[0:-4]+"slp.tif -fel "+filename)#unix case
        #removing unecessary Dinf slope tif file
        os.remove(filename[0:-4]+"slp.tif")
        #
        #Slope D8 file modifications
        tif_sd8=osgeo.gdal.Open(filename[0:-4]+"sd8.tif")#opening the file (1 band, bit not specified, but usually 32)
        if tif_sd8 is None:
          print (time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't open input sd8 dataset")
          sys.exit(1)
        else:
          print (time.strftime("%d/%m/%Y %H:%M:%S    "),"opening sd8 was successful!")
        #Array 
        tif_sd8_ar=tif_sd8.ReadAsArray()
        del(tif_sd8)
    else:
        sink_dir8 = path_sink_dtm + "/sinked_fdir8.tif"
        sink_dir8 = str(sink_dir8)
        tif_fdir8=osgeo.gdal.Open(sink_dir8)#opening the file (1 band, bit not specified, but usually 32)
        if tif_fdir8 is None:
          print (time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't open input sk_fd8 dataset")
          sys.exit(1)
        else:
          print (time.strftime("%d/%m/%Y %H:%M:%S    "), "opening sk_fd8 was successful!")
        #Array 
        tif_fdir8_ar=tif_fdir8.ReadAsArray()
        tif_fdir8_ar=tif_fdir8_ar.astype(numpy.float)#otherwise overflow in future operations
        ndv=numpy.min(tif_fdir8_ar)
        tif_fdir8_ar[tif_fdir8_ar==ndv]=numpy.NaN
        del(tif_fdir8)
        #os.remove(path_sink_dtm + "/sinked_fdir8.tif")
        #
        #Slope D8 sinked
        sink_sd8 = path_sink_dtm + "/sinked_sd8.tif"
        sink_sd8 = str(sink_sd8)
        tif_sd8=osgeo.gdal.Open(sink_sd8)
        if tif_sd8 is None:
          print (time.strftime("%d/%m/%Y %H:%M:%S    "),"couldn't open input sk_sd8 dataset")
          sys.exit(1)
        else:
          print (time.strftime("%d/%m/%Y %H:%M:%S    "),"opening sk_sd8 was successful!")
        #Array 
        tif_sd8_ar=tif_sd8.ReadAsArray()
        del(tif_sd8)
        os.remove(path_sink_dtm + "/sinked_sd8.tif")
        #
        #copying and renaming the d8 and dirinf file so that from here on it is as it was the usual ang Taudem file
        shutil.copy2(path_sink_dtm + "/sinked_dirinf.tif", filename[0:-4]+"ang.tif")
        shutil.copy2(path_sink_dtm + "/sinked_fdir8.tif", filename[0:-4]+"p.tif")
    
    #
    #imposing upper and lower limits to slope, no data here are -1
    tif_sd8_ar[(tif_sd8_ar>=0) & (tif_sd8_ar<0.005)]=0.005
    tif_sd8_ar[(tif_sd8_ar>1)]=1
    ndv=numpy.min(tif_sd8_ar)
    tif_sd8_ar[tif_sd8_ar==ndv]=numpy.NaN
    #
    #Create a reclassified slope D8 tiff array from numpy using gdal(better after removing it if it exists)
    dst_d8sl_ds = gdal.GetDriverByName('GTiff').Create((filename[0:-4]+'s.tif'), const_ar.shape[1], const_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    dst_d8sl_ds.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    dst_d8sl_ds.SetProjection(proj)
    dst_d8sl_ds.GetRasterBand(1).WriteArray(tif_sd8_ar,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
    dst_d8sl_ds=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del dst_d8sl_ds
    #
    #Selecting weighting factor
    filewgt=w_f.replace('\\','/')
    filewgt=str(filewgt)
    tif_wgt=osgeo.gdal.Open(filewgt)#opening the file (1 band, bit not specified, but usually 32)
    #
    # Opening Messages
    if tif_wgt is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't open input weight dataset")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "opening Weight was successful!")
    
    #Wgt factor matrix
    tif_wgt_ar=tif_wgt.ReadAsArray()
    del(tif_wgt)
    ndv=numpy.min(tif_wgt_ar)
    tif_wgt_ar[tif_wgt_ar==ndv]=numpy.NaN
    #
    #Computing 1/(W*S)
    Ws_1=1/(tif_wgt_ar*tif_sd8_ar)
    #
    #
    #Calculating flow acumulation in order to retrieve after the coordinates of the outlet
    #os.system("mpiexec -n 8 AreaD8 -p "+filename[0:-4]+"p.tif"+ " -ad8 "+filename[0:-4]+"ad8.tif -nc")#windows case
    os.system((("mpiexec -n 8 AreaD8 -p ").lower())+filename[0:-4]+"p.tif"+ " -ad8 "+filename[0:-4]+"ad8.tif -nc")#unix case
    #
    #
    #Trying to do a weighted flow length, USING FLOW DIR MATRIX, that removes the border cells
    #addind top row and left column with zeros
    #
    #
    #zero matrix bigger than F_dir8, to avoid border indexing problems
    Fd8 = numpy.zeros(shape=((tif_fdir8_ar.shape[0])+1,(tif_fdir8_ar.shape[1])+1),dtype=numpy.float32)
    Fd8[1:Fd8.shape[0],1:Fd8.shape[1]] = Fd8[1:Fd8.shape[0],1:Fd8.shape[1]]+tif_fdir8_ar
    #adding bottom row and right y axis with zeros
    Fdir8=numpy.zeros(shape=((Fd8.shape[0])+1,(Fd8.shape[1])+1),dtype=numpy.float32)
    Fdir8[:Fdir8.shape[0]-1,:Fdir8.shape[1]-1] = Fd8
    #
    #zero matrix bigger than weight, to avoid border indexing problems, and have same indexing as Fdir8
    Wg = numpy.zeros(shape=((tif_wgt_ar.shape[0])+1,(tif_wgt_ar.shape[1])+1),dtype=numpy.float32)
    Wg[1:Wg.shape[0],1:Wg.shape[1]] = Wg[1:Fd8.shape[0],1:Wg.shape[1]]+Ws_1#the weigth to weigth tha flow length
    #adding bottom row and right y axis with zeros
    Wgt=numpy.zeros(shape=((Wg.shape[0])+1,(Wg.shape[1])+1),dtype=numpy.float32)
    Wgt[:Wgt.shape[0]-1,:Wgt.shape[1]-1] = Wg
    #
    start =time.time()#for computational time
    #Creating a bigger matrix as large as weight(and all the matrices) to store the weighted flow length values
    W_Fl=numpy.zeros(shape=((Wgt.shape[0]),(Wgt.shape[1])),dtype=numpy.float32)
    W_Fl=W_Fl-1#to give -1 to NoData after the while loop calculation
    #
    #Let's go for the search and algo-rhytm for the weighted-Flow-Length
    ND=numpy.where(numpy.isnan(Fdir8)==True)#fast coordinates all the NoData values, starting from them to go forward and compute flow length
    #
    Y=ND[0]#rows, NoData indexes
    X=ND[1]#columns, NoData indexes pay attention not to invert values !!!!!!!!!!!!!!
    #
    #initializing lists for outlet and moving cell coordinates, in function of their position
    YC1=[];YC2=[];YC3=[];YC4=[];YC5=[];YC6=[];YC7=[];YC8=[]
    XC1=[];XC2=[];XC3=[];XC4=[];XC5=[];XC6=[];XC7=[];XC8=[]
    #
    #   Flow Directions Taudem
    #   4   3   2
    #   5   -   1
    #   6   7   8
    #
    #   Draining in Direction Matrix
    #   8   7   6
    #   1   -   5
    #   2   3   4
    #
    i1=Fdir8[Y,X-1]#Searching for NoData with cells draining into them, 8 directions
    D1=numpy.where(i1==1)#l
    YC1.extend(Y[D1])#coordinates satisfacting the conditions
    XC1.extend(X[D1])
    W_Fl[YC1,XC1]=0#initialize flow length at cells draining to NoData
    #
    i2=Fdir8[Y+1,X-1]#Searching for NoData with cells draining into them, 8 directions
    D2=numpy.where(i2==2)#lrad2
    YC2.extend(Y[D2])#coordinates satisfacting the conditions
    XC2.extend(X[D2])
    W_Fl[YC2,XC2]=0#initialize flow length at cells draining to NoData
    #
    i3=Fdir8[Y+1,X]#Searching for NoData with cells draining into them, 8 directions
    D3=numpy.where(i3==3)#l
    YC3.extend(Y[D3])#coordinates satisfacting the conditions
    XC3.extend(X[D3])
    W_Fl[YC3,XC3]=0#initialize flow length at cells draining to NoData
    #
    i4=Fdir8[Y+1,X+1]#Searching for NoData with cells draining into them, 8 directions
    D4=numpy.where(i4==4)#lrad2
    YC4.extend(Y[D4])#coordinates satisfacting the conditions
    XC4.extend(X[D4])
    W_Fl[YC4,XC4]=0#initialize flow length at cells draining to NoData
    #
    i5=Fdir8[Y,X+1]#Searching for NoData with cells draining into them, 8 directions
    D5=numpy.where(i5==5)#l
    YC5.extend(Y[D5])#coordinates satisfacting the conditions
    XC5.extend(X[D5])
    W_Fl[YC5,XC5]=0#initialize flow length at cells draining to NoData
    #
    i6=Fdir8[Y-1,X+1]#Searching for NoData with cells draining into them, 8 directions
    D6=numpy.where(i6==6)#lrad2
    YC6.extend(Y[D6])#coordinates satisfacting the conditions
    XC6.extend(X[D6])
    W_Fl[YC6,XC6]=0#initialize flow length at cells draining to NoData
    #
    i7=Fdir8[Y-1,X]#Searching for NoData with cells draining into them, 8 directions
    D7=numpy.where(i7==7)#l
    YC7.extend(Y[D7])#coordinates satisfacting the conditions
    XC7.extend(X[D7])
    W_Fl[YC7,XC7]=0#initialize flow length at cells draining to NoData
    #
    i8=Fdir8[Y-1,X-1]#Searching for NoData with cells draining into them, 8 directions
    D8=numpy.where(i8==8)#lrad2
    YC8.extend(Y[D8])#coordinates satisfacting the conditions
    XC8.extend(X[D8])
    W_Fl[YC8,XC8]=0#initialize flow length at cells draining to NoData
    #
    #start =time.time()#da cancellare poi.....!!!!!! Solo per check
    count=1# "0" passage already done during the previous step
    while len (YC1) or len(YC2) or len(YC3) or len(YC4) or len(YC5) or len(YC6) or len(YC7) or len(YC8) >0:
        #Converting into array to be able to do operations
        YYC1=numpy.asarray(YC1);XXC1=numpy.asarray(XC1)
        YYC2=numpy.asarray(YC2);XXC2=numpy.asarray(XC2)
        YYC3=numpy.asarray(YC3);XXC3=numpy.asarray(XC3)
        YYC4=numpy.asarray(YC4);XXC4=numpy.asarray(XC4)
        YYC5=numpy.asarray(YC5);XXC5=numpy.asarray(XC5)
        YYC6=numpy.asarray(YC6);XXC6=numpy.asarray(XC6)
        YYC7=numpy.asarray(YC7);XXC7=numpy.asarray(XC7)
        YYC8=numpy.asarray(YC8);XXC8=numpy.asarray(XC8)
        #
        #Now I can do operations and moving towards the right cell!!!!!!!! Weigthing flow length, weights are half sum of pixels weight * travelled length
        #I'm chosing the directions accordingly to Flow_dir step by step going from outlet-nodata to the ridges,
        #each time account for distance (l or l*rad2) multiplied by the half of the weigths of the 2 travelled cells.
        #Then, with variables substitution I'm moving a step further, and adding the prevous pixel value to the new calculated.
        #
        YYC1=(YYC1);XXC1=(XXC1-1)#l
        YYC2=(YYC2+1);XXC2=(XXC2-1)#lrad2
        YYC3=(YYC3+1);XXC3=(XXC3)#l
        YYC4=(YYC4+1);XXC4=(XXC4+1)#lrad2
        YYC5=(YYC5);XXC5=(XXC5+1)#l
        YYC6=(YYC6-1);XXC6=(XXC6+1)#lrad2
        YYC7=(YYC7-1);XXC7=(XXC7)#l
        YYC8=(YYC8-1);XXC8=(XXC8-1)#lrad2
        #
        if count==1:#first run zero, like TauDEM, need to check if there is a Nodata pixel receiving flow for all the 8 directions
            if len(YYC1)>0:
                W_Fl[YYC1,XXC1]=0
            else:
                pass
            if len(YYC2)>0:
                W_Fl[YYC2,XXC2]=0
            else:
                pass
            if len(YYC3)>0:
                W_Fl[YYC3,XXC3]=0
            else:
                pass
            if len(YYC4)>0:
                W_Fl[YYC4,XXC4]=0
            else:
                pass
            if len(YYC5)>0:
                W_Fl[YYC5,XXC5]=0
            else:
                pass
            if len(YYC6)>0:
                W_Fl[YYC6,XXC6]=0
            else:
                pass
            if len(YYC7)>0:
                W_Fl[YYC7,XXC7]=0
            else:
                pass
            if len(YYC8)>0:
                W_Fl[YYC8,XXC8]=0
            else:
                pass
        else:               
            W_Fl[YYC1,XXC1]=W_Fl[YC1,XC1]+(cell_s * ((Wgt[YC1,XC1] + Wgt[YYC1,XXC1])/2))
            W_Fl[YYC2,XXC2]=W_Fl[YC2,XC2]+(cell_s * math.sqrt(2) * ((Wgt[YC2,XC2] + Wgt[YYC2,XXC2])/2))   
            W_Fl[YYC3,XXC3]=W_Fl[YC3,XC3]+(cell_s * ((Wgt[YC3,XC3] + Wgt[YYC3,XXC3])/2))    
            W_Fl[YYC4,XXC4]=W_Fl[YC4,XC4]+(cell_s * math.sqrt(2) * ((Wgt[YC4,XC4] + Wgt[YYC4,XXC4])/2))    
            W_Fl[YYC5,XXC5]=W_Fl[YC5,XC5]+(cell_s * ((Wgt[YC5,XC5] + Wgt[YYC5,XXC5])/2))    
            W_Fl[YYC6,XXC6]=W_Fl[YC6,XC6]+(cell_s * math.sqrt(2) * ((Wgt[YC6,XC6] + Wgt[YYC6,XXC6])/2))    
            W_Fl[YYC7,XXC7]=W_Fl[YC7,XC7]+(cell_s * ((Wgt[YC7,XC7] + Wgt[YYC7,XXC7])/2))    
            W_Fl[YYC8,XXC8]=W_Fl[YC8,XC8]+(cell_s * math.sqrt(2) * ((Wgt[YC8,XC8] + Wgt[YYC8,XXC8])/2))
            #
        #
        #Reconstructing all X and Y of this step and moving on upwards (Downstream if you think in GIS, right?)
        YY=[];XX=[]
        YY.extend(YYC1);XX.extend(XXC1)
        YY.extend(YYC2);XX.extend(XXC2)
        YY.extend(YYC3);XX.extend(XXC3)
        YY.extend(YYC4);XX.extend(XXC4)
        YY.extend(YYC5);XX.extend(XXC5)
        YY.extend(YYC6);XX.extend(XXC6)
        YY.extend(YYC7);XX.extend(XXC7)
        YY.extend(YYC8);XX.extend(XXC8)
        #
        YY=numpy.asarray(YY)
        XX=numpy.asarray(XX)
        #
        i1=Fdir8[YY,XX-1]#Searching for cells draining into them, 8 directions
        D1=numpy.where(i1==1)#l
        YC1=YY[D1]#coordinates satisfacting the conditions, HERE i NEED TO ADD ACTUAL LENGTH VALUE + PREVIOUS ONE
        XC1=XX[D1]
        #
        i2=Fdir8[YY+1,XX-1]#Searching for cells draining into them, 8 directions
        D2=numpy.where(i2==2)#lrad2
        YC2=YY[D2]#coordinates satisfacting the conditions
        XC2=XX[D2]
        #
        i3=Fdir8[YY+1,XX]#Searching for cells draining into them, 8 directions
        D3=numpy.where(i3==3)#l
        YC3=YY[D3]#coordinates satisfacting the conditions
        XC3=XX[D3]
        #
        i4=Fdir8[YY+1,XX+1]#Searching for cells draining into them, 8 directions
        D4=numpy.where(i4==4)#lrad2
        YC4=YY[D4]#coordinates satisfacting the conditions
        XC4=XX[D4]
        #
        i5=Fdir8[YY,XX+1]#Searching for cells draining into them, 8 directions
        D5=numpy.where(i5==5)#l
        YC5=YY[D5]#coordinates satisfacting the conditions
        XC5=XX[D5]
        #
        i6=Fdir8[YY-1,XX+1]#Searching for cells draining into them, 8 directions
        D6=numpy.where(i6==6)#lrad2
        YC6=YY[D6]#coordinates satisfacting the conditions
        XC6=XX[D6]
        #
        i7=Fdir8[YY-1,XX]#Searching for cells draining into them, 8 directions
        D7=numpy.where(i7==7)#l
        YC7=YY[D7]#coordinates satisfacting the conditions
        XC7=XX[D7]
        #
        i8=Fdir8[YY-1,XX-1]#Searching for cells draining into them, 8 directions
        D8=numpy.where(i8==8)#lrad2
        YC8=YY[D8]#coordinates satisfacting the conditions
        XC8=XX[D8]
        count=count+1
    #
    #
    elapsed=(time.time()-start)#computational time
    print (time.strftime("%d/%m/%Y %H:%M:%S    "), "Process concluded succesfully \n","%.2f" % elapsed,'seconds for Weighted-Flow Length calculation with ',int(count),' iterations') #truncating the precision
    #os.system("pause")
    #
    #
    W_fl= W_Fl[1:W_Fl.shape[0]-1,1:W_Fl.shape[1]-1] #reshaping weigthed flow length, we need this step to homogenize matrices dimensions!!!!!!!!!!
    del W_Fl
    del Fdir8
    #
    w_fl_ds = gdal.GetDriverByName('GTiff').Create((dir_path+"/w_flow_length.tif"), const_ar.shape[1], const_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    w_fl_ds.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    w_fl_ds.SetProjection(proj)
    w_fl_ds.GetRasterBand(1).WriteArray(W_fl,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
    w_fl_ds=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del w_fl_ds
    #
    #
    #Downslope component, imposing 1 values where weigthed flow length==0, thus avoiding division by 0
    D_down_ar=W_fl
    del W_fl
    D_down_ar[D_down_ar==0]=1
    #
    #
    #
    """
    --------------------------------------
    WORKING ON D_UP COMPONENT
    --------------------------------------
    """
    #
    #Calculating D_inf flow acumulation to obtain dtmfilloksca
    #os.system("mpiexec -n 8 AreaDinf -ang "+filename[0:-4]+"ang.tif"+ " -sca "+filename[0:-4]+"sca.tif -nc")#windows case
    os.system((("mpiexec -n 8 AreaDinf -ang ").lower())+filename[0:-4]+"ang.tif"+ " -sca "+filename[0:-4]+"sca.tif -nc")#unix case
    #
    #dtmfilloksca file modifications--> conversion to array for manipulations
    tif_dtmsca=osgeo.gdal.Open(filename[0:-4]+"sca.tif")#opening the file (1 band, bit not specified, but usually 32)
    if tif_dtmsca is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't open input dtmsca dataset")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "opening dtmsca was successful!")
    #Array 
    tif_sca_ar=tif_dtmsca.ReadAsArray()
    del(tif_dtmsca)
    ndv=numpy.min(tif_sca_ar)
    tif_sca_ar[tif_sca_ar==ndv]=numpy.NaN
    #
    acc_final_ar=tif_sca_ar/const_ar
    #
    #Calculating accW-->Dinf contributing area wrigthed with weigth raster
    #os.system("mpiexec -n 8 AreaDinf -ang "+filename[0:-4]+"ang.tif"+ " -sca " +dir_path+"/accW.tif -wg "+filewgt+" -nc")#windows case
    os.system((("mpiexec -n 8 AreaDinf -ang ").lower())+filename[0:-4]+"ang.tif"+ " -sca " +dir_path+"/accW.tif -wg "+filewgt+" -nc")#unix case
    #
    acc_W=osgeo.gdal.Open(dir_path+"/accW.tif")#opening the file (1 band, bit not specified, but usually 32)
    if acc_W is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't open input acc_W dataset")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "opening acc_W was successful!")
    #Array 
    acc_W_ar=acc_W.ReadAsArray()
    del(acc_W)
    #
    #Calculating accS-->Dinf contributing area wrigthed with D8 slope value
    os.system("mpiexec -n 8 AreaDinf -ang "+filename[0:-4]+"ang.tif"+ " -sca " +dir_path+"/accS.tif -wg "+filename[0:-4]+'s.tif'+" -nc")#windows case
    os.system((("mpiexec -n 8 AreaDinf -ang ").lower())+filename[0:-4]+"ang.tif"+ " -sca " +dir_path+"/accS.tif -wg "+filename[0:-4]+'s.tif'+" -nc")#unix case
    #
    acc_S=osgeo.gdal.Open(dir_path+"/accS.tif")#opening the file (1 band, bit not specified, but usually 32)
    if acc_S is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't open input acc_S dataset")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "opening acc_S was successful!")
    #Array 
    acc_S_ar=acc_S.ReadAsArray()
    del(acc_S)
    #
    #
    #
    #
    #
    #
    #
    #Some cleanings removing taudem created files
    os.remove(dir_path+"/accS.tif")
    os.remove(filename[0:-4]+"ang.tif")
    os.remove(dir_path+"/w_flow_length.tif")
    os.remove(dir_path+"/accW.tif")
    os.remove(filename[0:-4]+"sca.tif")
    os.remove(filename[0:-4]+"ad8.tif")
    os.remove(filename[0:-4]+"s.tif")
    os.remove(filename[0:-4]+"p.tif")
    #
    if sink_flag ==1:        
        os.remove(dir_path+"/sinked_dtm.tif")
        os.remove(dir_path+"/sinked_fdir8.tif")
        #os.remove(dir_path+"/sinked_sd8.tif")
        os.remove(dir_path+"/sinked_dirinf.tif")
    else:
        os.remove(filename[0:-4]+"sd8.tif")
        pass
    #
    #
    #
    #
    #
    #
    #
    #Computing C_mean as (accW+weigth)/acc_final
    C_mean_ar=(acc_W_ar+tif_wgt_ar)/acc_final_ar
    del(acc_W_ar)#free memory
    #
    #Computing S mean (accS+s)/acc_final
    S_mean_ar=(acc_S_ar+tif_sd8_ar)/acc_final_ar
    del(acc_S_ar,tif_sd8_ar)#free memory
    #
    #Computing D_up as "%cmean.tif%" * "%smean.tif%" * SquareRoot("%ACCfinal.tif%" * "%resolution.tif%" * "%resolution.tif%")
    cell_area=(const_ar)**2#change of variables, to be sure
    D_up_ar=C_mean_ar*S_mean_ar*numpy.sqrt(acc_final_ar*cell_area)#to transform from unit values to square units
    #
    #Computing Connectivity index
    ic_ar=numpy.log10(D_up_ar/D_down_ar)
    #
    #saving upslope and downslope components if reqiured, filename is fixed
    if updown.isChecked():
        upsl_comp = gdal.GetDriverByName('GTiff').Create((dir_path+"/D_up.tif"), const_ar.shape[1], const_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
        upsl_comp.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
        upsl_comp.SetProjection(proj)
        upsl_comp.GetRasterBand(1).WriteArray(D_up_ar,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
        upsl_comp=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
        del upsl_comp
        del D_up_ar
        #
        dwsl_comp = gdal.GetDriverByName('GTiff').Create((dir_path+"/D_down.tif"), const_ar.shape[1], const_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
        dwsl_comp.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
        dwsl_comp.SetProjection(proj)
        dwsl_comp.GetRasterBand(1).WriteArray(D_down_ar,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
        dwsl_comp=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
        del dwsl_comp
        del D_down_ar
    else:
        del D_up_ar
        del D_down_ar
    #
    #Create the tiff raster of ic
    o_ic=out_ic.replace('\\','/')
    o_ic=str(o_ic)
    dst_ic_ds = gdal.GetDriverByName('GTiff').Create(o_ic, const_ar.shape[1], const_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    dst_ic_ds.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    dst_ic_ds.SetProjection(proj)
    dst_ic_ds.GetRasterBand(1).WriteArray(ic_ar,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
    dst_ic_ds.GetRasterBand(1).GetStatistics(0,1) #calculate statistics for visualization
    dst_ic_ds=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del dst_ic_ds
    #
    print (time.strftime("%d/%m/%Y %H:%M:%S    "), "Calculation finished!, Compliments!")
    #END?....YUPPY!!!!!
    #
    # Figure
    imgplot = matplotlib.pyplot.imshow(ic_ar)
    imgplot.set_cmap('Paired')
    matplotlib.pyplot.colorbar()
    matplotlib.pyplot.xlabel('x', fontsize=18)
    matplotlib.pyplot.ylabel('y', fontsize=18)
    matplotlib.pyplot.suptitle('Connectivity output map', fontsize=18)
    matplotlib.pyplot.show()
    
    # easygui.msgbox("Calculation finished!", title="Compliments")
    return   

"""
--------------------------------------
###########################################################################       

CONNECTIVITY TO TARGETS

###########################################################################
--------------------------------------
"""
def CavalliConnectivitytg (dtm_f, c_s, tg_f, w_f, out_ic_tg, sink_flag_tg = 0):
     #
     #file selection
     #
    filename=dtm_f.replace('\\','/')
    filename=str(filename)
    tif=osgeo.gdal.Open(filename)#opening the file (1 band, bit not specified, but usually 32)
    #
    #folder path
    dir_path=os.path.dirname(os.path.realpath(filename))#path of the selected input file
    dir_path=dir_path.replace('\\','/')
    #
    # Opening Messages
    if tif is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't open input DTM dataset (Connectivity to targets)")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "opening DTM (Connectivity to targets) was successful!")
    #
    #Colonne, righe, bande
    cols = tif.RasterXSize
    rows = tif.RasterYSize
    bands =  tif.RasterCount
    #
    #dtm as an array                #################################       ATTENZIONE     ##########################
    tif_ar=tif.ReadAsArray() #con GDAL 1.6 e 1.7 non funzia 90%, con gdal bindings per windows 1.9.2 sembra funzionare provare anche a metter su ultime versioni di numpy
    tif_ar = tif_ar.astype(float)
    band=tif.GetRasterBand(1)#bands retrieving
    geoinf = tif.GetGeoTransform()#upper left coordinates, cellsize
    proj=tif.GetProjection()
    tipo=tif_ar.dtype
    del(tif)
    #
    #cell_size
    cell_s=numpy.float32(str(c_s))#importante che resti float32, per omogeneita'
    #
    #create constant array to trasform into raster
    const_ar=(tif_ar>=0)*cell_s#in base all'array true/false, costruisco il raster costante
    c1=numpy.where(const_ar==0)
    const_ar[c1]=-1#sostituisco il valore -1 per i NoData
    #
    #target file selection
    #
    NoData_value=-9999
    file_tar = tg_f.replace('\\','/')
    file_tar = str(file_tar)
    temp_tar = gdal.GetDriverByName('GTiff').Create((dir_path+"/targets.tif"), const_ar.shape[1], const_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    temp_tar.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    temp_tar.SetProjection(proj)
    band = temp_tar.GetRasterBand(1)
    band.SetNoDataValue(NoData_value)
    
    targ_ds = ogr.Open(file_tar)
    targ_layer = targ_ds.GetLayer()
    
    gdal.RasterizeLayer(temp_tar, [1], targ_layer, None, None, burn_values=[10])#,options= ['ALL_TOUCHED=TRUE'])
    temp_tar=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del temp_tar
    
    #tg=poly2ras.rasterize(file_tar, dir_path+"/targets.tif", prototype=filename, burn_values=[10], options=["MERGE_ALG=ADD"])#assign value to variable or arbitrary one
    tg = osgeo.gdal.Open(dir_path+"/targets.tif")
    target_ar=tg.GetRasterBand(1).ReadAsArray()
    target_ar=target_ar.astype(float)
    #target_ar[numpy.where(tif_ar == tif_ar.min())] = numpy.NaN #where input has nodata, also the target must have, useful in case of sinks that clip the target too, no need for connectivity to the outlet
    del (tg)
    #
    #
    #create constant array to trasform into Fl length raster
    Fl_ar=(tif_ar>=0)*1#in base all'array true/false, costruisco il raster costante
    c1=numpy.where(Fl_ar==0)#1 data, -1 NoData
    Fl_ar[c1]=-1
    #
    if sink_flag_tg == 0:
        #D8 flow directions & outputs
        #os.system("mpiexec -n 8 D8Flowdir -p "+filename[0:-4]+"p.tif"+ " -sd8 "+filename[0:-4]+"sd8.tif -fel "+filename)#windows case
        os.system((("mpiexec -n 8 D8Flowdir -p ").lower())+filename[0:-4]+"p.tif"+ " -sd8 "+filename[0:-4]+"sd8.tif -fel "+filename)#unix case
        #Dinf flow directions & outputs
        #os.system("mpiexec -n 8 DinfFlowdir -ang "+filename[0:-4]+"angt.tif"+ " -slp "+filename[0:-4]+"slp.tif -fel "+filename)#windows case
        os.system((("mpiexec -n 8 DinfFlowdir -ang ").lower())+filename[0:-4]+"angt.tif"+ " -slp "+filename[0:-4]+"slp.tif -fel "+filename)#unix case
        #removing unecessary Dinf slope tif file
        os.remove(filename[0:-4]+"slp.tif")
        #
    else:
        shutil.copy2(path_sink_dtm + "/sinked_dirinf.tif", filename[0:-4]+"angt.tif")
        shutil.copy2(path_sink_dtm + "/sinked_fdir8.tif", filename[0:-4]+"p.tif")
        shutil.copy2(path_sink_dtm + "/sinked_sd8.tif", filename[0:-4]+"sd8.tif")
        #
    """
    ---------------------------------------------
    WORKING ON D_DOWN COMPONENT
    ---------------------------------------------
    """
    #opening D8 flow directions
    tif_fdir8=osgeo.gdal.Open(filename[0:-4]+"p.tif")#opening the file (1 band, bit not specified, but usually 32)
    if tif_fdir8 is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "couldn't open input fdir_8 dataset")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "opening fdir_8 was successful!")
    #Array 
    tif_fdir8_ar=tif_fdir8.ReadAsArray()
    tif_fdir8_ar=tif_fdir8_ar.astype(numpy.float)
    ndv=numpy.min(tif_fdir8_ar)
    tif_fdir8_ar[tif_fdir8_ar==ndv]=numpy.NaN
    del (tif_fdir8)
    #Now_setting_NoData_where I have targets==10
    tif_fdir8_ar[(target_ar==10)]=-1000#I will then consider all negative values as NoData
    #
    #opening Dinf flow directions
    tif_fdirinf=osgeo.gdal.Open(filename[0:-4]+"angt.tif")#opening the file (1 band, bit not specified, but usually 32)
    if tif_fdirinf is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "),"couldn't open input dirinf_tg dataset")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "),"opening dirinf_tg was successful!")
    #Array 
    tif_fdirinf_ar=tif_fdirinf.ReadAsArray()
    ndv=numpy.min(tif_fdirinf_ar)
    tif_fdirinf_ar[tif_fdirinf_ar==ndv]=-9999
    del (tif_fdirinf) #otherwise it is in use and will not be deleted
    #Now_setting_NoData_where I have targets==10
    tif_fdirinf_ar[(target_ar==10)]=-1000#I will then consider all negative values as NoData
    #removing  Dinf  tif file, and rewrite it, the one with NoData
    #os.remove(filename[0:-4]+"angt.tif")
    #creating the new Dinf
    FD_inf_ds = gdal.GetDriverByName('GTiff').Create((filename[0:-4]+'ang.tif'), const_ar.shape[1], const_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    FD_inf_ds.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    FD_inf_ds.SetProjection(proj)
    FD_inf_ds.GetRasterBand(1).SetNoDataValue(-9999) 
    FD_inf_ds.GetRasterBand(1).WriteArray(tif_fdirinf_ar,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
    FD_inf_ds=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del FD_inf_ds
    #
    #
    #
    #Slope D8 file modifications
    tif_sd8=osgeo.gdal.Open(filename[0:-4]+"sd8.tif")#opening the file (1 band, bit not specified, but usually 32)
    if tif_sd8 is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "),"couldn't open input sd8 dataset")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "),"opening sd8 was successful!")
    #Array 
    tif_sd8_ar=tif_sd8.ReadAsArray()
    del (tif_sd8)
    #
    #imposing upper and lower limits to slope, no data here are -1
    tif_sd8_ar[(tif_sd8_ar>=0) & (tif_sd8_ar<0.005)]=0.005
    tif_sd8_ar[(tif_sd8_ar>1)]=1
    tif_sd8_ar[(tif_sd8_ar<0)]=-1
    #
    #Create a reclassified slope D8 tiff array from numpy using gdal(better after removing it if it exists)
    dst_d8sl_ds = gdal.GetDriverByName('GTiff').Create((filename[0:-4]+'s.tif'), const_ar.shape[1], const_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    dst_d8sl_ds.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    dst_d8sl_ds.SetProjection(proj)
    dst_d8sl_ds.GetRasterBand(1).WriteArray(tif_sd8_ar,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
    dst_d8sl_ds=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del dst_d8sl_ds
    #
    #Selecting weighting factor
    filewgt=w_f.replace('\\','/')
    filewgt=str(filewgt)
    tif_wgt=osgeo.gdal.Open(filewgt)#opening the file (1 band, bit not specified, but usually 32)
    #
    # Opening Messages
    if tif_wgt is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "),"couldn't open input Weight dataset")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "),"opening Weight was successful!")
    
    #Wgt factor matrix
    tif_wgt_ar=tif_wgt.ReadAsArray()
    del (tif_wgt)
    #
    #Computing 1/(W*S)
    Ws_1=1/(tif_wgt_ar*tif_sd8_ar)
    #
    #
    #Calculating flow acumulation in order to retrieve after the coordinates of the outlet
    #os.system("mpiexec -n 8 AreaD8 -p "+filename[0:-4]+"p.tif"+ " -ad8 "+filename[0:-4]+"ad8.tif -nc")#windows case
    os.system((("mpiexec -n 8 AreaD8 -p ").lower())+filename[0:-4]+"p.tif"+ " -ad8 "+filename[0:-4]+"ad8.tif -nc")#unix case
    #
    #
    #Trying to do a weighted flow length, USING FLOW DIR MATRIX, that removes the border cells
    #addind top row and left column with zeros
    #
    start =time.time()#for computational time
    #
    #zero matrix bigger than F_dir8, to avoid border indexing problems
    Fd8 = numpy.zeros(shape=((tif_fdir8_ar.shape[0])+1,(tif_fdir8_ar.shape[1])+1),dtype=numpy.float32)
    Fd8[1:Fd8.shape[0],1:Fd8.shape[1]] = Fd8[1:Fd8.shape[0],1:Fd8.shape[1]]+tif_fdir8_ar
    #adding bottom row and right y axis with zeros
    Fdir8=numpy.zeros(shape=((Fd8.shape[0])+1,(Fd8.shape[1])+1),dtype=numpy.int16)
    Fdir8[:Fdir8.shape[0]-1,:Fdir8.shape[1]-1] = Fd8
    #
    #zero matrix bigger than weight, to avoid border indexing problems, and have same indexing as Fdir8
    Wg = numpy.zeros(shape=((tif_wgt_ar.shape[0])+1,(tif_wgt_ar.shape[1])+1),dtype=numpy.float32)
    Wg[1:Wg.shape[0],1:Wg.shape[1]] = Wg[1:Fd8.shape[0],1:Wg.shape[1]]+Ws_1#the weigth to weigth tha flow length
    #adding bottom row and right y axis with zeros
    Wgt=numpy.zeros(shape=((Wg.shape[0])+1,(Wg.shape[1])+1),dtype=numpy.float32)
    Wgt[:Wgt.shape[0]-1,:Wgt.shape[1]-1] = Wg
    #
    #Creating a bigger matrix as large as weight(and all the matrices) to store the weighted flow length values
    W_Fl=numpy.zeros(shape=((Wgt.shape[0]),(Wgt.shape[1])),dtype=numpy.float32)
    W_Fl=W_Fl-1#to give -1 to NoData after the while loop calculation
    #
    #Let's go for the search and algo-rhytm for the weighted-Flow-Length
    ND=numpy.where(Fdir8<0)#fast coordinates all the NoData values, starting from them to go forward and compute flow length
    #
    Y=ND[0]#rows, NoData indexes
    X=ND[1]#columns, NoData indexes pay attention not to invert values !!!!!!!!!!!!!!
    #
    #initializing lists for outlet and moving cell coordinates, in function of their position
    YC1=[];YC2=[];YC3=[];YC4=[];YC5=[];YC6=[];YC7=[];YC8=[]
    XC1=[];XC2=[];XC3=[];XC4=[];XC5=[];XC6=[];XC7=[];XC8=[]
    #
    #   Flow Directions Taudem
    #   4   3   2
    #   5   -   1
    #   6   7   8
    #
    #   Draining in Direction Matrix
    #   8   7   6
    #   1   -   5
    #   2   3   4
    #
    #
    i1=Fdir8[Y,X-1]#Searching for NoData with cells draining into them, 8 directions
    D1=numpy.where(i1==1)#l
    YC1.extend(Y[D1])#coordinates satisfacting the conditions
    XC1.extend(X[D1])
    W_Fl[YC1,XC1]=0#initialize flow length at cells draining to NoData
    #
    i2=Fdir8[Y+1,X-1]#Searching for NoData with cells draining into them, 8 directions
    D2=numpy.where(i2==2)#lrad2
    YC2.extend(Y[D2])#coordinates satisfacting the conditions
    XC2.extend(X[D2])
    W_Fl[YC2,XC2]=0#initialize flow length at cells draining to NoData
    #
    i3=Fdir8[Y+1,X]#Searching for NoData with cells draining into them, 8 directions
    D3=numpy.where(i3==3)#l
    YC3.extend(Y[D3])#coordinates satisfacting the conditions
    XC3.extend(X[D3])
    W_Fl[YC3,XC3]=0#initialize flow length at cells draining to NoData
    #
    i4=Fdir8[Y+1,X+1]#Searching for NoData with cells draining into them, 8 directions
    D4=numpy.where(i4==4)#lrad2
    YC4.extend(Y[D4])#coordinates satisfacting the conditions
    XC4.extend(X[D4])
    W_Fl[YC4,XC4]=0#initialize flow length at cells draining to NoData
    #
    i5=Fdir8[Y,X+1]#Searching for NoData with cells draining into them, 8 directions
    D5=numpy.where(i5==5)#l
    YC5.extend(Y[D5])#coordinates satisfacting the conditions
    XC5.extend(X[D5])
    W_Fl[YC5,XC5]=0#initialize flow length at cells draining to NoData
    #
    i6=Fdir8[Y-1,X+1]#Searching for NoData with cells draining into them, 8 directions
    D6=numpy.where(i6==6)#lrad2
    YC6.extend(Y[D6])#coordinates satisfacting the conditions
    XC6.extend(X[D6])
    W_Fl[YC6,XC6]=0#initialize flow length at cells draining to NoData
    #
    i7=Fdir8[Y-1,X]#Searching for NoData with cells draining into them, 8 directions
    D7=numpy.where(i7==7)#l
    YC7.extend(Y[D7])#coordinates satisfacting the conditions
    XC7.extend(X[D7])
    W_Fl[YC7,XC7]=0#initialize flow length at cells draining to NoData
    #
    i8=Fdir8[Y-1,X-1]#Searching for NoData with cells draining into them, 8 directions
    D8=numpy.where(i8==8)#lrad2
    YC8.extend(Y[D8])#coordinates satisfacting the conditions
    XC8.extend(X[D8])
    W_Fl[YC8,XC8]=0#initialize flow length at cells draining to NoData
    #
    #start =time.time()#da cancellare poi.....!!!!!! Solo per check
    count=1# "0" passage already done during the previous step
    while len (YC1) or len(YC2) or len(YC3) or len(YC4) or len(YC5) or len(YC6) or len(YC7) or len(YC8) >0:
        #Converting into array to be able to do operations
        YYC1=numpy.asarray(YC1);XXC1=numpy.asarray(XC1)
        YYC2=numpy.asarray(YC2);XXC2=numpy.asarray(XC2)
        YYC3=numpy.asarray(YC3);XXC3=numpy.asarray(XC3)
        YYC4=numpy.asarray(YC4);XXC4=numpy.asarray(XC4)
        YYC5=numpy.asarray(YC5);XXC5=numpy.asarray(XC5)
        YYC6=numpy.asarray(YC6);XXC6=numpy.asarray(XC6)
        YYC7=numpy.asarray(YC7);XXC7=numpy.asarray(XC7)
        YYC8=numpy.asarray(YC8);XXC8=numpy.asarray(XC8)
        #
        #Now I can do operations and moving towards the right cell!!!!!!!! Weigthing flow length, weights are half sum of pixels weight * travelled length
        #I'm chosing the directions accordingly to Flow_dir step by step going from outlet-nodata to the ridges,
        #each time account for distance (l or l*rad2) multiplied by the half of the weigths of the 2 travelled cells.
        #Then, with variables substitution I'm moving a step further, and adding the prevous pixel value to the new calculated.
        #
        YYC1=(YYC1);XXC1=(XXC1-1)#l
        YYC2=(YYC2+1);XXC2=(XXC2-1)#lrad2
        YYC3=(YYC3+1);XXC3=(XXC3)#l
        YYC4=(YYC4+1);XXC4=(XXC4+1)#lrad2
        YYC5=(YYC5);XXC5=(XXC5+1)#l
        YYC6=(YYC6-1);XXC6=(XXC6+1)#lrad2
        YYC7=(YYC7-1);XXC7=(XXC7)#l
        YYC8=(YYC8-1);XXC8=(XXC8-1)#lrad2
        #
        if count==1:#first run zero, like TauDEM
            if len(YYC1)>0:
                W_Fl[YYC1,XXC1]=0
            else:
                pass
            if len(YYC2)>0:
                W_Fl[YYC2,XXC2]=0
            else:
                pass
            if len(YYC3)>0:
                W_Fl[YYC3,XXC3]=0
            else:
                pass
            if len(YYC4)>0:
                W_Fl[YYC4,XXC4]=0
            else:
                pass
            if len(YYC5)>0:
                W_Fl[YYC5,XXC5]=0
            else:
                pass
            if len(YYC6)>0:
                W_Fl[YYC6,XXC6]=0
            else:
                pass
            if len(YYC7)>0:
                W_Fl[YYC7,XXC7]=0
            else:
                pass
            if len(YYC8)>0:
                W_Fl[YYC8,XXC8]=0
            else:
                pass
        else:               
            W_Fl[YYC1,XXC1]=W_Fl[YC1,XC1]+(cell_s * ((Wgt[YC1,XC1] + Wgt[YYC1,XXC1])/2))
            W_Fl[YYC2,XXC2]=W_Fl[YC2,XC2]+(cell_s * math.sqrt(2) * ((Wgt[YC2,XC2] + Wgt[YYC2,XXC2])/2))   
            W_Fl[YYC3,XXC3]=W_Fl[YC3,XC3]+(cell_s * ((Wgt[YC3,XC3] + Wgt[YYC3,XXC3])/2))    
            W_Fl[YYC4,XXC4]=W_Fl[YC4,XC4]+(cell_s * math.sqrt(2) * ((Wgt[YC4,XC4] + Wgt[YYC4,XXC4])/2))    
            W_Fl[YYC5,XXC5]=W_Fl[YC5,XC5]+(cell_s * ((Wgt[YC5,XC5] + Wgt[YYC5,XXC5])/2))    
            W_Fl[YYC6,XXC6]=W_Fl[YC6,XC6]+(cell_s * math.sqrt(2) * ((Wgt[YC6,XC6] + Wgt[YYC6,XXC6])/2))    
            W_Fl[YYC7,XXC7]=W_Fl[YC7,XC7]+(cell_s * ((Wgt[YC7,XC7] + Wgt[YYC7,XXC7])/2))    
            W_Fl[YYC8,XXC8]=W_Fl[YC8,XC8]+(cell_s * math.sqrt(2) * ((Wgt[YC8,XC8] + Wgt[YYC8,XXC8])/2))
            #
        #Reconstructing all X and Y of this step and moving on upwards (Downstream if you think in GIS, right?)
        YY=[];XX=[]
        YY.extend(YYC1);XX.extend(XXC1)
        YY.extend(YYC2);XX.extend(XXC2)
        YY.extend(YYC3);XX.extend(XXC3)
        YY.extend(YYC4);XX.extend(XXC4)
        YY.extend(YYC5);XX.extend(XXC5)
        YY.extend(YYC6);XX.extend(XXC6)
        YY.extend(YYC7);XX.extend(XXC7)
        YY.extend(YYC8);XX.extend(XXC8)
        #
        YY=numpy.asarray(YY)
        XX=numpy.asarray(XX)
        #
        i1=Fdir8[YY,XX-1]#Searching for cells draining into them, 8 directions
        D1=numpy.where(i1==1)#l
        YC1=YY[D1]#coordinates satisfacting the conditions, HERE i NEED TO ADD ACTUAL LENGTH VALUE + PREVIOUS ONE
        XC1=XX[D1]
        #
        i2=Fdir8[YY+1,XX-1]#Searching for cells draining into them, 8 directions
        D2=numpy.where(i2==2)#lrad2
        YC2=YY[D2]#coordinates satisfacting the conditions
        XC2=XX[D2]
        #
        i3=Fdir8[YY+1,XX]#Searching for cells draining into them, 8 directions
        D3=numpy.where(i3==3)#l
        YC3=YY[D3]#coordinates satisfacting the conditions
        XC3=XX[D3]
        #
        i4=Fdir8[YY+1,XX+1]#Searching for cells draining into them, 8 directions
        D4=numpy.where(i4==4)#lrad2
        YC4=YY[D4]#coordinates satisfacting the conditions
        XC4=XX[D4]
        #
        i5=Fdir8[YY,XX+1]#Searching for cells draining into them, 8 directions
        D5=numpy.where(i5==5)#l
        YC5=YY[D5]#coordinates satisfacting the conditions
        XC5=XX[D5]
        #
        i6=Fdir8[YY-1,XX+1]#Searching for cells draining into them, 8 directions
        D6=numpy.where(i6==6)#lrad2
        YC6=YY[D6]#coordinates satisfacting the conditions
        XC6=XX[D6]
        #
        i7=Fdir8[YY-1,XX]#Searching for cells draining into them, 8 directions
        D7=numpy.where(i7==7)#l
        YC7=YY[D7]#coordinates satisfacting the conditions
        XC7=XX[D7]
        #
        i8=Fdir8[YY-1,XX-1]#Searching for cells draining into them, 8 directions
        D8=numpy.where(i8==8)#lrad2
        YC8=YY[D8]#coordinates satisfacting the conditions
        XC8=XX[D8]
        count=count+1
    #
    #
    elapsed=(time.time()-start)#computational time
    print (time.strftime("%d/%m/%Y %H:%M:%S    "), "Process concluded succesfully \n","%.2f" % elapsed,'seconds for Weighted-Flow Length calculation with ',int(count),' iterations') #truncating the precision
    #os.system("pause")
    #
    #
    W_fl= W_Fl[1:W_Fl.shape[0]-1,1:W_Fl.shape[1]-1] #reshaping weigthed flow length, we need this step to homogenize matrices dimensions!!!!!!!!!!
    del W_Fl
    del Fdir8
    #
    w_fl_ds = gdal.GetDriverByName('GTiff').Create((dir_path+"/w_flow_length.tif"), const_ar.shape[1], const_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    w_fl_ds.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    w_fl_ds.SetProjection(proj)
    w_fl_ds.GetRasterBand(1).WriteArray(W_fl,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
    w_fl_ds=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del w_fl_ds
    #
    #
    #Downslope component, imposing 1 values where weigthed flow length==0, thus avoiding division by 0
    D_down_ar=W_fl
    del W_fl
    D_down_ar[D_down_ar==0]=1
    D_down_ar[numpy.where(target_ar==10)]=numpy.NaN#repeating NOData setting where I have the target, otherwise it will misplace values in the output raster... Not really sure actually
    #
    #
    #
    """
    --------------------------------------
    WORKING ON D_UP COMPONENT
    --------------------------------------
    """
    #
    #Calculating D_inf flow acumulation to obtain dtmfilloksca
    #os.system("mpiexec -n 8 AreaDinf -ang "+filename[0:-4]+"ang.tif"+ " -sca "+filename[0:-4]+"sca.tif -nc")#windows case
    os.system((("mpiexec -n 8 AreaDinf -ang ").lower())+filename[0:-4]+"ang.tif"+ " -sca "+filename[0:-4]+"sca.tif -nc")#unix case
    #
    #dtmfilloksca file modifications--> conversion to array for manipulations
    tif_dtmsca=osgeo.gdal.Open(filename[0:-4]+"sca.tif")#opening the file (1 band, bit not specified, but usually 32)
    if tif_dtmsca is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "),"couldn't open input tif_dtmsca dataset")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "), "opening tif_dtmsca was successful!")
    #Array 
    tif_sca_ar=tif_dtmsca.ReadAsArray()
    #
    acc_final_ar=tif_sca_ar/const_ar
    #
    #Calculating accW-->Dinf contributing area wrigthed with weigth raster
    #os.system("mpiexec -n 8 AreaDinf -ang "+filename[0:-4]+"ang.tif"+ " -sca " +dir_path+"/accW.tif -wg "+filewgt+" -nc")#windows case
    os.system((("mpiexec -n 8 AreaDinf -ang ").lower())+filename[0:-4]+"ang.tif"+ " -sca " +dir_path+"/accW.tif -wg "+filewgt+" -nc")#unix case
    #
    acc_W=osgeo.gdal.Open(dir_path+"/accW.tif")#opening the file (1 band, bit not specified, but usually 32)
    if tif_dtmsca is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "),"couldn't open input acc_W dataset")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "),"opening acc_W was successful!")
    #Array 
    acc_W_ar=acc_W.ReadAsArray()
    del(acc_W)
    #
    #Calculating accS-->Dinf contributing area wrigthed with D8 slope value
    #os.system("mpiexec -n 8 AreaDinf -ang "+filename[0:-4]+"ang.tif"+ " -sca " +dir_path+"/accS.tif -wg "+filename[0:-4]+'s.tif'+" -nc")#windows case
    os.system((("mpiexec -n 8 AreaDinf -ang ").lower())+filename[0:-4]+"ang.tif"+ " -sca " +dir_path+"/accS.tif -wg "+filename[0:-4]+'s.tif'+" -nc")#unix case
    #
    acc_S=osgeo.gdal.Open(dir_path+"/accS.tif")#opening the file (1 band, bit not specified, but usually 32)
    if tif_dtmsca is None:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "),"couldn't open input acc_S dataset")
      sys.exit(1)
    else:
      print (time.strftime("%d/%m/%Y %H:%M:%S    "),"opening acc_S was successful!")
    #Array 
    acc_S_ar=acc_S.ReadAsArray()
    del (acc_S)
    del (tif_dtmsca)
    #
    #Computing C_mean as (accW+weigth)/acc_final
    C_mean_ar=(acc_W_ar+tif_wgt_ar)/acc_final_ar
    #
    #Computing S mean (accS+s)/acc_final
    S_mean_ar=(acc_S_ar+tif_sd8_ar)/acc_final_ar
    #
    #Computing D_up as "%cmean.tif%" * "%smean.tif%" * SquareRoot("%ACCfinal.tif%" * "%resolution.tif%" * "%resolution.tif%")
    #i need a trick to elevate power of 2 anly the positive values of the constant resolution array, otherwise all the -1 will become 1
    in_pos_cos=numpy.where(const_ar>=0)
    const_ar[in_pos_cos]=(const_ar[in_pos_cos])**2
    cell_area=const_ar#change of variables, to be sure
    D_up_ar=C_mean_ar*S_mean_ar*numpy.sqrt(acc_final_ar*cell_area)#to transform from unit values to square units
    D_up_ar[numpy.where(target_ar==10)]=numpy.NaN#repeating NoData setting where I have the target, otherwise it will misplace values in the output raster... Not really sure actually
    #
    #Computing Connectivity index
    ic_ar=numpy.log10(D_up_ar/D_down_ar)
    #ic_ar[ic_ar==0]=numpy.nan#setting NoData values, maybe incorrect for occurrences within the matrix
    #
    #saving upslope and downslope components if reqiured, filename is fixed
    if updown.isChecked():
        upsl_comp = gdal.GetDriverByName('GTiff').Create((dir_path+"/D_up.tif"), const_ar.shape[1], const_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
        upsl_comp.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
        upsl_comp.SetProjection(proj)
        upsl_comp.GetRasterBand(1).WriteArray(D_up_ar,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
        upsl_comp=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
        del upsl_comp
        del D_up_ar
        #
        dwsl_comp = gdal.GetDriverByName('GTiff').Create((dir_path+"/D_down.tif"), const_ar.shape[1], const_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
        dwsl_comp.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
        dwsl_comp.SetProjection(proj)
        dwsl_comp.GetRasterBand(1).WriteArray(D_down_ar,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
        dwsl_comp=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
        del dwsl_comp
        del D_down_ar
    else:
        del D_up_ar
        del D_down_ar
    #
    #
    o_ic_tg=out_ic_tg.replace('\\','/')
    o_ic_tg=str(o_ic_tg)
    #
    #
    #Create the tiff raster of ic
    dst_ic_ds = gdal.GetDriverByName('GTiff').Create(o_ic_tg, const_ar.shape[1], const_ar.shape[0], 1, GDT_Float32)#shape sono righe e colonne, 1 è il numero di bande
    dst_ic_ds.GetRasterBand(1).SetNoDataValue(-9999) 
    dst_ic_ds.SetGeoTransform(geoinf)#conferisco coordinate e proiezione del raster in input
    dst_ic_ds.SetProjection(proj)
    dst_ic_ds.GetRasterBand(1).WriteArray(ic_ar,0,0)#scrivo effettivamente il raster, prima avevo solo allocato la memoria
    dst_ic_ds.GetRasterBand(1).GetStatistics(0,1) #calculate statistics for visualization
    dst_ic_ds=None#Sometimes I add this to ensure the file is fully deallocated, and to prevent running into some gotchas:
    del dst_ic_ds
    #
    # Figure
    imgplot = matplotlib.pyplot.figure(1)
    matplotlib.pyplot.imshow(ic_ar)
    matplotlib.pyplot.set_cmap('Paired')
    matplotlib.pyplot.colorbar()
    matplotlib.pyplot.xlabel('x', fontsize=18)
    matplotlib.pyplot.ylabel('y', fontsize=18)
    matplotlib.pyplot.suptitle('Connectivity output map', fontsize=18)
    matplotlib.pyplot.show()
    #
    # Figure_2
    dist_data = ic_ar[numpy.isnan(ic_ar)==False]
    #
    imgplot = matplotlib.pyplot.figure(2)
    matplotlib.pyplot.hist(dist_data, bins=150, normed=True)
    matplotlib.pyplot.grid()
    matplotlib.pyplot.xlabel('IC values', fontsize=18)
    matplotlib.pyplot.ylabel('F', fontsize=18)
    matplotlib.pyplot.suptitle('Distribution and values', fontsize=18)
    mu = dist_data.mean()
    median = numpy.median(dist_data)
    sigma = dist_data.std()
    min_data = dist_data.min()
    max_data = dist_data.max()
    textstr = '$min=%.4f$\n$max=%.4f$\n$\mu=%.4f$\n$\mathrm{median}=%.4f$\n$\sigma=%.4f$'%(min_data,max_data,mu, median, sigma)
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    matplotlib.pyplot.text(0.45, 0.45, textstr,  fontsize=14,
            verticalalignment='top', bbox=props)
    matplotlib.pyplot.show()
    #
    #Some cleanings removing taudem created files
    os.remove(dir_path+"/accS.tif")
    os.remove(filename[0:-4]+"ang.tif")
    os.remove(filename[0:-4]+"angt.tif")
    os.remove(dir_path+"/w_flow_length.tif")
    os.remove(dir_path+"/accW.tif")
    os.remove(filename[0:-4]+"sca.tif")
    os.remove(filename[0:-4]+"ad8.tif")
    os.remove(filename[0:-4]+"s.tif")
    os.remove(filename[0:-4]+"p.tif")
    os.remove(filename[0:-4]+"sd8.tif")
    os.remove(dir_path+"/targets.tif")
    #
    if sink_flag_tg ==1:
        os.remove(dir_path+"/sinked_dtm.tif")
        os.remove(dir_path+"/sinked_fdir8.tif")
        os.remove(dir_path+"/sinked_sd8.tif")
        os.remove(dir_path+"/sinked_dirinf.tif")

    else:
        pass
    #
    print (time.strftime("%d/%m/%Y %H:%M:%S    "), "Calculation finished!, Compliments!")
    # easygui.msgbox("Calculation finished!", title="Compliments")
    return