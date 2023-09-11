# Image-GIS-Projection
Python CLI tool to project polygons from an image (e.g. object bounding boxes) to GIS (shapefile), and projection in the other direction, GIS to image coordinates. Batch functionality also implemented.

[Link to orthomosaic for provided examples](https://drive.google.com/file/d/1bFHmdtsY0fvG47YHJ8e_5pslcQkR7r6U/view?usp=drive_link)

## Usage
Required non-standard python packages: rasterio, geopandas, pyproj, numpy, yaml, shapely, pandas, pprint, scikit-image

`
python main.py --config_file config.yml --img2gis --batch
`
* `CFG_PATH`: path to config file
* `--img2gis`: project image coordinates to shapefile
* `--gis2img`: project opposite direction
* `--batch: flag to toggle batch processing

Example for single-image projection of image coordinates to GIS.
`
python main.py --config_file config.yml --img2gis --batch
`



### Config File
Settings can be changed in config.yml file, default values are provided for the provided examples. 
Note that provided DSM should have altitude in the same reference system (geodetic/ellipsoidal) as in the exif headers of the images OR the `cameras.csv`. 

`cameras.csv` can be obtained for a flight from Agisoft and helps provide more precise projection, but can be omitted if not available (leave blank in config file).
