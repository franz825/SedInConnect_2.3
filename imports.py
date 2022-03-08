import sys, os
from time import time, sleep
import time
# import tkinter as Tk
# import tkinter.messagebox
# import easygui
# from easygui import *
import osgeo
from osgeo import gdal, ogr
import numpy
import math
import struct
import shutil
from osgeo.gdal import*
# from PyQt4 import QtGui, QtCore
# from PyQt4.QtGui import QLineEdit
#from PIL import Image
import threading
#lib_path=os.path.abspath('C:/Users/stefano/Dropbox/qt_python_test/')
#lib_path2=os.path.abspath('C:/Python27_64/Lib/site-packages/scipy/signal/')
#sys.path.append(lib_path)
#sys.path.append(lib_path2)
# import guimages_2_3
import gdal_rasterize as poly2ras
#import progressbar
import base64
import scipy.signal
from scipy.special import _ufuncs_cxx #this for dependencies, otherwise pyinstaller does not find _ufuncs_cxx trying to create scipy dll
import matplotlib.pyplot
import FileDialog

print("hello world")
