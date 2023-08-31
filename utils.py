"""
Utility functions for projection, file handling, data manipulation, etc.
"""

import pandas as pd
import numpy as np
from ast import literal_eval
import os
import skimage.io
from pyproj import CRS,  Proj, transform
import urllib.request
from osgeo import gdal
import shutil
from tqdm import tqdm
import yaml

def read_cfg(cfg_path):
    ''' Read INI file and return ConfigParser object with user-set options '''
    with open(cfg_path) as f:
        cfg = yaml.safe_load(f)
    return cfg

def getWKT_PRJ(epsg_code):
    # access projection information
    with urllib.request.urlopen("http://spatialreference.org/ref/epsg/{0}/prettywkt/".format(epsg_code)) as wkt:
        content = wkt.read().decode('utf-8')
        # remove spaces between charachters
        remove_spaces = content.replace(" ","")
        # place all the text on one line
        output = remove_spaces.replace("\n", "")
    return output

def extract_exif(img_path):
    # extract exif info to a temp file
    call = f"exiftool -a {img_path} > exif_tmp.txt"
    os.system(call)
    exif_dict = {}
    with open(f"exif_tmp.txt", 'r') as f:
        lines = f.readlines()
        for line in lines:
            line = line.split(':')
            key = line[0].strip()
            val = line[1].strip()
            exif_dict[key] = val

    return exif_dict

def convert_gps(gps):
    parts = gps.split(' ')
    res = 0.0
    res += float(parts[0])
    res += float(parts[2][:-1]) / 60
    res += float(parts[3][:-1]) / 3600

    if parts[4] in ['W','S']:
        res *= -1
    return res

def rot_matrix(K, O, P):
    # convert to radians
    K = np.deg2rad(K)
    O = np.deg2rad(O)
    P = np.deg2rad(P)

    # compute 
    m = np.array([
        [np.cos(P)*np.cos(K),                                -np.cos(P)*np.sin(K),                                  np.sin(P)],
        [np.cos(O)*np.sin(K) + np.sin(O)*np.sin(P)*np.cos(K), np.cos(O)*np.cos(K) - np.sin(O)*np.sin(P)*np.sin(K), -np.sin(O)*np.cos(P)],
        [np.sin(O)*np.sin(K) - np.cos(O)*np.sin(P)*np.cos(K), np.sin(O)*np.cos(K) + np.cos(O)*np.sin(P)*np.sin(K),  np.cos(O)*np.cos(P)]
    ])

    return m


