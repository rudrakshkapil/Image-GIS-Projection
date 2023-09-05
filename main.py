""" 
Main script implementing the CLI
    Reads arguments, config file, determines which functions to call, 
    ... reads inputs, calls projection functions as required, creates output file(s). 

Actual implementation of the projection is found in /projection.py
"""

import argparse
from ast import literal_eval

from matplotlib import pyplot as plt
from projection import *
from pprint import pprint
import shapefile
import pandas as pd
import geopandas as gpd
import shapely
shapely.speedups.disable()
from shapely.geometry import mapping
import rasterio

# command line arguments
parser = argparse.ArgumentParser(description='Image <--> GIS Projection Tool')
group_mode = parser.add_mutually_exclusive_group(required=True)
group_mode.add_argument('--img2gis', action='store_true', help='Desired direction of projection. Either one of these flags is required')
group_mode.add_argument('--gis2img', action='store_true', help='Forward: img2gis, Backward: gis2img')
parser.add_argument('--config_file', action='store', type=str, required=True, help='You must provide a path to the config file with the required settings')
parser.add_argument('--batch', action='store_true', help='Toggle for batch vs single image processing')
parser.add_argument('--bare_ground', action='store_true', help='Toggle for using canopy-level points (using DSM) or ground level (using DTM)')

# ------------ SINGLE IMAGE IMG2GIS ------------ #
def single_img2gis(cfg):
    # get image path, csv_path, dsm_path
    img_path = cfg['paths']['single_cfg']['img_path']
    dsm_path = cfg['paths']['dsm_path']
    points_csv = cfg['paths']['single_cfg']['img2gis']['points_csv']
    points_df = pd.read_csv(points_csv)

    # load dsm
    dsm = rasterio.open(dsm_path)
    dsm_arr = dsm.read(1)
    try: dsm_epsg = int(str(dsm.crs)[5:])
    except: exit("Error: DSM does not seem to have CRS information. Chcek and try again")
    

    # get output file path and make enclosing dir if not exists
    output_file = f"{cfg['paths']['single_cfg']['img2gis']['out_shp']}"
    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))

    # create output file(s)
    w = shapefile.Writer(output_file, shapefile.POLYGON)
    w.field('id','C','40')
    if cfg['centroids']:
        w_centroid = shapefile.Writer(output_file+'_centroids', shapefile.POINT)
        w_centroid.field('id','C','40')


    # loop over rows (boxes in same image)
    for idx, row in tqdm(points_df.iterrows(), total=len(points_df)):
        # convert points to numpy array and prepend centroid
        points = np.array(literal_eval(row['geometry']), dtype=np.uint16)
        xlim = (np.amin(points[:,0]), np.amax(points[:,0]))
        ylim = (np.amin(points[:,1]), np.amax(points[:,1]))
        center = np.array([(xlim[0]+xlim[1])/2, (ylim[0]+ylim[1])/2]).reshape(-1,2)
        points = np.vstack([center, points])
        
        # project to GIS
        curr_poly, _ = img2gis(img_path, points, dsm, dsm_arr, cfg) # NOTE: remove depth

        # add to shapefile
        w.poly([curr_poly.tolist()[1:]])
        w.record(row['id'], 'Polygon')

        if cfg['centroids']:
            w_centroid.point(*curr_poly[0])
            w_centroid.record(row['id'], 'Point')
                
    # close shapefile, create .prj file
    w.close()
    wkt = getWKT_PRJ(dsm_epsg)
    prj = open(f"{output_file}.prj", "w")
    prj.write(wkt) # TODO: config file for EPSG (only do this if None)
    prj.close()

    # close centroid shapefile, create .prj file
    if cfg['centroids']:
        w.close()
        prj = open(f"{output_file}_centroids.prj", "w")
        prj.write(wkt)
        prj.close()


def batch_img2gis(cfg):
    # get image path, csv_path, dsm_path
    img_dir = cfg['paths']['batch_cfg']['img_dir']
    dsm_path = cfg['paths']['dsm_path']
    points_csv = cfg['paths']['batch_cfg']['img2gis']['points_csv']
    points_df = pd.read_csv(points_csv)

    dsm = rasterio.open(dsm_path)
    dsm_arr = dsm.read(1)
    try: dsm_epsg = int(str(dsm.crs)[5:])
    except: exit("Error: DSM does not seem to have CRS information. Chcek and try again")
    

    # get output file path and make enclosing dir if not exists
    output_file = f"{cfg['paths']['batch_cfg']['img2gis']['out_shp']}"
    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))

    # create output file(s) with required cols
    w = shapefile.Writer(output_file, shapefile.POLYGON)
    w.field('id','C','20')
    w.field('img_path', 'C', '40')
    w.field('box_id', 'C', '20')
    w.field('distance_to_img_center', 'F', decimal=3)
    if cfg['centroids']:
        w_centroid = shapefile.Writer(output_file+'_centroids', shapefile.POINT)
        w_centroid.field('id','C','20')
        w_centroid.field('img_path', 'C', '40')
        w_centroid.field('box_id', 'C', '20')
        w_centroid.field('distance_to_img_center', 'F', '20')


    # loop over rows (boxes in same image)
    for idx, row in tqdm(points_df.iterrows(), total=len(points_df)):
        # convert points to numpy array and prepend centroid
        points = np.array(literal_eval(row['geometry']), dtype=np.uint16)
        xlim = (np.amin(points[:,0]), np.amax(points[:,0]))
        ylim = (np.amin(points[:,1]), np.amax(points[:,1]))
        center = np.array([(xlim[0]+xlim[1])/2, (ylim[0]+ylim[1])/2]).reshape(-1,2)
        points = np.vstack([center, points])

        # get current img_path
        img_path = f"{img_dir}/{row['img_path']}"
        box_id = row['box_id']
        
        # project to GIS
        curr_poly, dist = img2gis(img_path, points, dsm, dsm_arr, cfg) 

        # add to shapefile
        w.poly([curr_poly.tolist()[1:]])
        w.record(row['id'], row['img_path'], box_id, dist, 'Polygon') # TODO: add elevation of box, add height (val from dsm) as Z, add distance to image it was, + elevation as a field

        if cfg['centroids']:
            w_centroid.point(*curr_poly[0])
            w_centroid.record(row['id'], row['img_path'], box_id, dist, 'Point')
                
    # close shapefile, create .prj file
    w.close()
    prj = open(f"{output_file}.prj", "w")
    wkt = getWKT_PRJ(dsm_epsg)
    prj.write(wkt)
    prj.close()

    # close centroid shapefile, create .prj file
    if cfg['centroids']:
        w.close()
        prj = open(f"{output_file}_centroids.prj", "w")
        prj.write(wkt)
        prj.close()



def single_gis2img(cfg):
    # get image path, csv_path, dsm_path
    img_path = cfg['paths']['single_cfg']['img_path']
    dsm_path = cfg['paths']['dsm_path']
    points_shp = cfg['paths']['single_cfg']['gis2img']['points_shp'] + '.shp'

    # read shapefile
    shapefile = gpd.read_file(points_shp)
    geoms = shapefile.geometry.values # list of shapely geometries
    geoms = [[mapping(g)] for g in geoms]
    shapefile_dict = shapefile.to_dict(orient='records')

    dsm = rasterio.open(dsm_path)
    dsm_arr = dsm.read(1)
    try: dsm_epsg = int(str(dsm.crs)[5:])
    except: exit("Error: DSM does not seem to have CRS information. Chcek and try again")
    
    # row_ids
    row_ids = shapefile['id'].tolist()
    

    # get output file path and make enclosing dir if not exists
    output_file = f"{cfg['paths']['single_cfg']['gis2img']['out_csv']}"
    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))


    # create dataframe for csv of image points (attributes copied from df)
    shp_cols = list(shapefile.columns)
    csv_cols = ['id', 'image_points'].extend(shp_cols)
    df = pd.DataFrame(columns=csv_cols)
    all_rows = []


    # loop over boxes
    for idx, geom in enumerate(geoms):
        # get coordinates from row (-1 to get rid of last (== first))
        poly = np.array(geom[0]['coordinates'][0][:-1])

        # append centroid
        xlim = (np.amin(poly[:,0]), np.amax(poly[:,0]))
        ylim = (np.amin(poly[:,1]), np.amax(poly[:,1]))
        center = np.array([(xlim[0]+xlim[1])/2, (ylim[0]+ylim[1])/2]).reshape(-1,2)
        poly = np.vstack([center, poly])

        
        # get image coordinates of polygon vertices
        pts = gis2img(img_path, poly, dsm, dsm_arr, cfg)[0].tolist()
        pts = [tuple(pt) for pt in pts[1:]] # skip centroid

        # get current row with 3 new attrs and all existing ones (remove 'geometry' -- UTM coords)
        row = {'id':row_ids[idx], 'geometry':str(pts)}
        # row.update(shapefile.iloc[idx])
        # del row['geometry']
        all_rows.append(row)
    df = pd.DataFrame(all_rows)
    df.to_csv(output_file)

    # plot
    img = skimage.io.imread(img_path)
    plt.imshow(img)
    for idx, row in df.iterrows():
        id = row['id']
        pts = literal_eval(row['geometry'])
        pts.append(pts[0])     
        pts = np.array(pts)
        # pts += [100,-100]
        plt.plot(pts[:,0],pts[:,1],c='r')
    plt.show()

    

def batch_gis2img(cfg):
    # get image path, csv_path, dsm_path
    img_dir = cfg['paths']['batch_cfg']['img_dir']
    dsm_path = cfg['paths']['dsm_path']
    points_shp = cfg['paths']['batch_cfg']['gis2img']['points_shp'] + '.shp'
    
    dsm = rasterio.open(dsm_path)
    dsm_arr = dsm.read(1)
    try: dsm_epsg = int(str(dsm.crs)[5:])
    except: exit("Error: DSM does not seem to have CRS information. Chcek and try again")

    # load image paths
    img_names = [fname for fname in os.listdir(img_dir) if fname[-4:].lower() in ['.tif','.jpg','tiff','jpeg']]
    img_paths = [f'{img_dir}/{fname}' for fname in img_names]
    # exit()

    # get all image extents
    extents = get_extent_of_all_images(img_paths, dsm, dsm_arr, cfg)
    extents = np.array(extents)

    # read shapefile
    shapefile = gpd.read_file(points_shp)
    geoms = shapefile.geometry.values # list of shapely geometries
    geoms = [[mapping(g)] for g in geoms]
    shapefile_dict = shapefile.to_dict(orient='records')
    
    # row_ids
    row_ids = shapefile['id'].tolist()
    

    # get output file path and make enclosing dir if not exists
    output_file = f"{cfg['paths']['batch_cfg']['gis2img']['out_csv']}"
    if not os.path.exists(os.path.dirname(output_file)):
        os.makedirs(os.path.dirname(output_file))


    # create dataframe for csv of image points (attributes copied from df)
    shp_cols = list(shapefile.columns)
    csv_cols = ['id', 'img_path', 'image_points'].extend(shp_cols)
    df = pd.DataFrame(columns=csv_cols)
    all_rows = []


    # loop over boxes
    for idx, geom in enumerate(geoms):
        # get coordinates from row (-1 to get rid of last (== first))
        poly = np.array(geom[0]['coordinates'][0][:-1])
        valid_list = get_candidate_img_list(extents, poly)

        # append centroid
        xlim = (np.amin(poly[:,0]), np.amax(poly[:,0]))
        ylim = (np.amin(poly[:,1]), np.amax(poly[:,1]))
        center = np.array([(xlim[0]+xlim[1])/2, (ylim[0]+ylim[1])/2]).reshape(-1,2)
        poly = np.vstack([center, poly])

        # loop over valid list images, do reverse projection
        for img_i in valid_list:
            # get image coordinates of polygon vertices
            pts, exif_dict = gis2img(img_paths[img_i], poly, dsm, dsm_arr, cfg)

            # confirm validity (widht and height should be contained)
            h, w = float(exif_dict['Exif Image Height']), float(exif_dict['Exif Image Width']),
            if np.any(pts[:,0] < 0) or np.any(pts[:,0] >= w): continue
            if np.any(pts[:,1] < 0) or np.any(pts[:,1] >= h): continue

            # convert format for valid points
            pts = [tuple(pt) for pt in pts.tolist()[1:]] # skip centroid
            
            # get current row with 3 new attrs and all existing ones (remove 'geometry' -- UTM coords)
            row = {'id':row_ids[idx], 'img_path':img_names[img_i], 'geometry':str(pts).replace(' ', '')}
            # row.update(shapefile.iloc[idx])
            # del row['geometry']
            all_rows.append(row)
    df = pd.DataFrame(all_rows)
    df.to_csv(output_file)

    df_id = df.groupby('img_path')['id'].apply(list).reset_index(name='id')
    df_geometry = df.groupby('img_path')['geometry'].apply(list).reset_index(name='geometry')
    df_id = df_id.assign(geometry=df_geometry['geometry'])

    # plot
    for idx, row in df_id.iterrows():
        img_path = f"{img_dir}/{row['img_path']}"
        ids, geometries = row['id'], row['geometry']

        img = skimage.io.imread(img_path)
        plt.imshow(img)

        for id, geometry in zip(ids, geometries):
            pts = literal_eval(geometry)
            pts.append(pts[0])     
            pts = np.array(pts)
            # pts += [100,-100]
            plt.plot(pts[:,0],pts[:,1], label=id)
        plt.legend()
        plt.show()



# -------------------- MAIN -------------------- #
if __name__ == "__main__":
    # parse arguments
    args = parser.parse_args()    

    # load config
    cfg = read_cfg(args.config_file)

    # batch processing
    if args.batch:
        print("Starting Batch Projection with passed config file:")
        pprint(cfg)
        # run in desired direction
        if args.img2gis:
            batch_img2gis(cfg)
        else: 
            batch_gis2img(cfg)


    # single-file processing
    else:
        print("Starting Single Image Projection with passed config file:")
        pprint(cfg)

        # run in desired direction
        if args.img2gis:
            single_img2gis(cfg)
        else: 
            single_gis2img(cfg)

        

        
            

