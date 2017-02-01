""" This code attempts to read a *.xml file and load the annotation from *.xml 
and extract fixed size tiles from *.tiff file. (Andy_Atypia Project) """

import openslide
from openslide import open_slide, ImageSlide
from openslide.deepzoom import DeepZoomGenerator
from optparse import OptionParser
import os
import glob
import re
import shutil
import math
import sys
from unicodedata import normalize
import matplotlib.pyplot as plt
import pickle
from PIL import ImageDraw
import scipy.cluster.hierarchy as hcluster
import xml.etree.ElementTree as ET

Tile_W = 2000
Tile_H = 2000
TilesPath = 'R:/Beck Lab/Atypia_Andy/LowAgreement/Tiles/'
os.chdir('R:/Beck Lab/Atypia_Andy/LowAgreement')
for WSI in glob.glob("*.tif"):
    WSIName = os.path.splitext(WSI)[0]
    image = openslide.open_slide(WSI)
    [WSI_W,WSI_H] = image.dimensions
    tree = ET.parse(WSIName+'.xml')
    root = tree.getroot()
    ROI_No = 1
    for Annotation in root.findall('Annotation'):
        Regions = Annotation.find('Regions')
        Region = Regions.find('Region')
        Text = Region.find('Text')
        ROIName = Text.get('Value')
        Vertices = Region.find('Vertices')
        
        X1 = int(Vertices[0].get('X'))
        ROI_X = int(math.ceil(X1/1000))*1000
        
        X2 = int(Vertices[1].get('X'))
        ROI_W = abs(X2 - ROI_X)
        #ROI_W = int(round(float(X2 - ROI_X)/1000)*1000)
    
        Y1 = int(Vertices[0].get('Y'))
        ROI_Y = int(math.ceil(Y1/1000))*1000
        
        Y2 = int(Vertices[1].get('Y'))
        #ROI_H = int(round(float(Y2 - ROI_Y)/1000)*1000)
        ROI_H = abs(Y2 - ROI_Y)
        
        ROI_W = int(round(ROI_W,-3))
        ROI_H = int(round(ROI_H,-3))
        
        ROI_X2 = ROI_X + ROI_W
        ROI_Y2 = ROI_Y + ROI_H
        
        if ROI_X <= 0:
            ROI_X = 1
        if ROI_Y <= 0:
            ROI_Y = 1
        if ROI_X2 > WSI_W:
            ROI_X2 = WSI_W
        if ROI_Y2 > WSI_H:
            ROI_Y2 = WSI_H
        print(WSIName, ROIName, WSI_W, WSI_H, ROI_X, ROI_X2, ROI_Y, ROI_Y2)
        x = ROI_X
        Tile_No = 0
        while x + Tile_W <= ROI_X2:
            y = ROI_Y
            while y + Tile_H <= ROI_Y2:
                Tile_No = Tile_No + 1
                Tile = image.read_region((x,y), int(0), (Tile_W,Tile_H))
                TileName = TilesPath + WSIName + '_' + ROIName + '_' + str(ROI_No) + '_' + str(x) + '_' + str(y) + '_' + str(Tile_No) + '.png'
                Tile.save(TileName)
                y = y + Tile_H
            x = x + Tile_W
        ROI_No = ROI_No + 1
 
    
