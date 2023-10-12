# Image-GIS-Projection
Python CLI tool to project polygons from an image (e.g. object bounding boxes) to GIS (shapefile), and projection in the other direction, GIS to image coordinates. Batch functionality also implemented.

[Link to orthomosaic for provided examples](https://drive.google.com/file/d/1bFHmdtsY0fvG47YHJ8e_5pslcQkR7r6U/view?usp=drive_link)


## Installation Instructions (for Windows)

1. Create conda environment (named `img2gis` here for example) with desired python version. **Note**: 3.8.13 is strongly recommended as that's what I programmed and tested with. 
```PowerShell
conda create --name img2gis python=3.8.13
conda activate img2gis
```

2. Download GDAL by going to the following repository: [Christoph Gohlke's Unofficial Windows Binaries for Python Extension Packages](https://github.com/cgohlke/geospatial-wheels/releases/tag/v2023.1.10.1) and downloading the version matching your Windows system (e.g. `GDAL-3.6.2-cp38-cp38-win_amd64.whl`
for 64 bit AMD computers, already included in this repo). Make sure it says cp38 for python version 3.8.

3. Place the downloaded `.whl` file in this repository (already done for example above). 

4. Install the file with pip
```PowerShell
pip install GDAL-3.6.2-cp38-cp38-win32.whl
 OR
pip install GDAL-3.6.2-cp38-cp38-win_amd64.whl
```

5. Install most of the other package requirements using pip, as listed in `requirements.txt`
```PowerShell
pip install -r requirements.txt
```

6. Finally, install shapely (through Conda, not pip to avoid conflicts)
```PowerShell
conda install -y shapely==2.0.1
```

And now you're done! :)

**Note**: for other OS, the only change would be from steps 2-5 onwards, i.e., need to download and install GDAL correctly.


## Usage Instructions

```PowerShell
usage: main.py [-h] (--img2gis | --gis2img) --config_file CONFIG_FILE [--batch] [--bare_ground]
```

```python
Image <--> GIS Projection Tool

optional arguments:
  -h, --help            show this help message and exit
  --img2gis             Desired direction of projection. Either one of these flags is required
  --gis2img             Forward: img2gis, Backward: gis2img
  --config_file CONFIG_FILE
                        You must provide a path to the config file with the required settings
  --batch               Toggle for batch vs single image processing
  --bare_ground         Toggle for using canopy-level points (using DSM) or ground level (using DTM)
```

##  Example Usage

For the provided data in `/examples` and the default settings in `config.yml`
Example for single-image projection of image coordinates to GIS.
```PowerShell
python main.py --config_file config.yml --img2gis --batch
```


### Config File

Settings can be changed in config.yml file, default values are provided for the provided examples. 
Note that provided DSM should have altitude in the same reference system (geodetic/ellipsoidal) as in the exif headers of the images OR the `cameras.csv`. 

`cameras.csv` can be obtained for a flight from Agisoft and helps provide more precise projection, but can be omitted if not available (leave blank in config file).
