""" 
Script that implements the core functionality for the projection
Comprises two functions:
    1. img2gis() for a given image and a given set of points (e.g. one polygon)
    2. gis2img() for a given shapefile and a given image
"""

from tqdm import tqdm
import rasterio
import numpy as np
import matplotlib.pyplot as plt
from utils import *
from pyproj import Proj
import pandas as pd
import skimage


def img2gis(img_path:str, 
            points:np.array([np.uint16]), 
            dsm:rasterio.DatasetReader, # output of rasterio.opnen()
            dsm_arr:np.array,
            cfg:dict,
            full_image:bool=False
            ) -> [np.array([np.float32]), tuple([int, str])]: 
    """ 
    Projection from Image Coordinates to Georeferenced Shapefile 
    Inputs:
        - img_path: relative/absolute path to image 
        - points: numpy array of N image pixel location points to be projected to UTM. First point must be the centroid. Array is shaped (N,2), ordered clockwise for correct shape in shapefile.  
        - dsm_path: path to digital surface model (DSM) for canopy projection or digital terrain model (DTM) for ground point projection. 
                This dsm just needs to correspond to the area where the image footprint is expected to be.
        - config: dictionary of config file information. Contains the following main information.
            * calibration_dict: dict containing camera calibration info for the given image. 
                    Must contain keys <'f' in pixels, 'cx' in pixels, 'cy' in pixels, 'sw' (sensor width in meters), 'sh' (sensor_height in meters)>. 
                    Any of the following can be included for more precise results, but are assumed to be 0 if not present:
                        Radial distortion <'k1','k2','k3','k4'>
                        Tangential distortion <'p1','p2'>
                        Skew <'b1','b2'>
            * csv_path (recommended): path to tab-separated csv with headers of aligned images if available, e.g. from Agisoft. 
                    Expected headers: <PhotoID,X,Y,Z,Omega,Phi,Kappa,r11,r12,r13,r21,r22,r23,r31,r32,r33> 
                    Provides more precise results than relying on only exif headers for these values, but if not given, exif info is what the function uses. 
            * show_box(optional): when True, shows patch contained by given points in given image (as a box using min-max u-v limits). 
        - full_image: if True, replaces points with full extent of image – useful for getting ground footprint
    Output:
        - curr_poly: numpy array of points projected to UTM coordinates
        - tuple of (zone_num, zone_letter): additional geo-location information
    """

    # extract required path & settings for processng $$
    cameras_csv_path = cfg['paths'].get('aligned_cameras_csv_path', None)

    show_box = cfg['img2gis'].get('show_box',False)
    coarse_del_Z = cfg['img2gis'].get('coarse_height_increment', 0.5) # starting increment for H in coarse to fine strategy of ray tracing/dsm interception
    fine_del_H = cfg['img2gis'].get('fine_height_increment', 0.5)     # starting increment for H in coarse to fine strategy of ray tracing/dsm interception
    coarse_tolerance = cfg['img2gis'].get('coarse_tolerance', 2)      # when we want to move to finer increments
    fine_tolerance = cfg['img2gis'].get('fine_tolerance', 0.2)        # how close we want to approximate (±value is tolerable) 
    bound = cfg['img2gis'].get('boundary', 5)                         # get dsm value as median in a small bxb region centered at px and py and curr_H:

    # extract calbiration information (e.g. given for DJI P1 in comments)
    calibration_dict = cfg['camera_calibration']
    f = calibration_dict['f']          # 11017.3
    cx = calibration_dict['cx']        # -6.35398
    cy = calibration_dict['cy']        # 22.0852
    sw = calibration_dict['sw']        # 36.0448 / 1000
    sh = calibration_dict['sh']        # 24.024 / 1000
    k1 = calibration_dict.get('k1', 0) # -0.00761618
    k2 = calibration_dict.get('k2', 0) # -0.245679 
    k3 = calibration_dict.get('k3', 0) # -0.520738
    k4 = calibration_dict.get('k4', 0) # 0.0
    p1 = calibration_dict.get('p1', 0) # -0.000425435
    p2 = calibration_dict.get('p2', 0) # 0.000371409
    b1 = calibration_dict.get('b1', 0) # -8.81372
    b2 = calibration_dict.get('b2', 0) # -0.708353

    # extract exif header
    exif_dict = extract_exif(img_path)

    # extract epsg
    try: dsm_epsg = int(str(dsm.crs)[5:])
    except: exit("Error: DSM does not seem to have CRS information. Chcek and try again")

    
    # get relevant information exif
    Z0 = float(exif_dict['Absolute Altitude']) 
    altitude = float(exif_dict['Relative Altitude'])              # increase -> box bigger
    image_height = float(exif_dict['Exif Image Height'])
    image_width = float(exif_dict['Exif Image Width'])
    focal = float(exif_dict['Focal Length'].split(' ')[0]) / 1000 # increase -> box bigger

    # GPS from LRF values directly + convert lat,long to UTM
    gps_lat = convert_gps(exif_dict["GPS Latitude"])
    gps_lon = convert_gps(exif_dict["GPS Longitude"])
    p = Proj(dsm_epsg, preserve_units=True)
    X0, Y0 = p(gps_lon, gps_lat)


    # rotation matrix m - NOTE: different drones/sensors have diffferent configurations for these, so may need to update this. 
    yaw = float(exif_dict['Gimbal Yaw Degree']) 
    roll = float(exif_dict['Gimbal Roll Degree']) 
    pitch = float(exif_dict['Gimbal Pitch Degree'])
    (kappa, omega, phi) = (yaw, pitch+90, roll)     # +90 needed due to DJI convention of gimbal vs camera view
    m = rot_matrix(kappa, omega, phi)

    # if csv path is given, get more precise estimates of m and Z0
    if cameras_csv_path is not None:
        df = pd.read_table(cameras_csv_path, skiprows=1)
        df.columns = ['PhotoID','X','Y','Z','Omega','Phi','Kappa','r11','r12','r13','r21','r22','r23','r31','r32','r33']

        img_name = img_path.split('/')[-1].split('.')[0] 
        row = df.loc[df['PhotoID'] == img_name]
        assert len(row) > 0, 'Error! Image not found in csv file.'

        omega = float(row['Omega'])
        phi = float(row['Phi'])
        kappa = float(row['Kappa'])
        m = rot_matrix(kappa, omega, phi)

        # X and Y used from exif, Z from csv
        Z0 = float(row['Z'])
        # X0 = float(row['X'])
        # Y0 = float(row['Y'])

    # replace points with image corners if needed
    if full_image:
        points = np.array([[image_width//2+cx,image_height//2+cy],[0,0],[image_width,0],[image_width,image_height],[0,image_height]])

    # show image
    if show_box:
        img = skimage.io.imread(img_path)
        xlim = (np.amin(points[:,0]), np.amax(points[:,0]))
        ylim = (np.amin(points[:,1]), np.amax(points[:,1]))
        plt.imshow(img[ylim[0]:ylim[1], xlim[0]:xlim[1]])
        plt.title("Object encricled by given points.")
        plt.show()

    
    # -------- Iterative 'Reverse' Projection with Undistortion ---------- #  
    # to go from  UV (2D pixel space coordinates) ---> undistorted XYZ (3D local image space coordinates)

    # make intrinsic matrix and inverse
    k = np.array([[f+b1, b2,  cx+image_width//2],
                  [   0,  f, cy+image_height//2],
                  [   0,  0,                 1]])
    k_inv = np.linalg.inv(k)

    # if projecting full image (rather than box), change points to be the four corners
    full_extent = np.array([[image_width//2+cx,image_height//2+cy],[0,0],[image_width,0],[image_width,image_height],[0,image_height]])

    # utility for iterating for best estimate of x,y (image local coords at z=1)
    def iterate_undistort(points):
        points = points.reshape(-1,2)
        points = np.hstack([points, np.ones((len(points),1))])
        xy_prime = k_inv.dot(points.T)[:2,:]

        r = np.sqrt((xy_prime[0]**2 + xy_prime[1]**2))  # r_prime
        # x_0 and y_0 init -> initial undistorted image coordinates
        x = xy_prime[0] - xy_prime[0] * (k1 * r**2 + k2 * r**4 + k3 * r**6 + k4 * r**8) - p1 * (r**2 + 2 * xy_prime[0]**2) - 2 * p2 * xy_prime[0] * xy_prime[1] 
        y = xy_prime[1] - xy_prime[1] * (k1 * r**2 + k2 * r**4 + k3 * r**6 + k4 * r**8) - p2 * (r**2 + 2 * xy_prime[1]**2) - 2 * p1 * xy_prime[0] * xy_prime[1] 
        # iterate to get x,y -> final undistorted image coordinates
        for _ in range(10):
            r = np.sqrt((x**2 + y**2))
            x = xy_prime[0] - x * (k1 * r**2 + k2 * r**4 + k3 * r**6 + k4 * r**8) - p1 * (r**2 + 2 * x**2) - 2 * p2 * x * y
            y = xy_prime[1] - y * (k1 * r**2 + k2 * r**4 + k3 * r**6 + k4 * r**8) - p2 * (r**2 + 2 * y**2) - 2 * p1 * x * y

        return x,y
    
    x,y = iterate_undistort(points)
    xf,yf = iterate_undistort(full_extent)

    # get limits, scale xy --> sensor space location of centroid & principal point (m: meters)
    xmin, xmax = min(xf), max(xf)
    ymin, ymax = min(yf), max(yf)
    xf_m = (xf-xmin)/(xmax-xmin) * sw  # f: indicates full extent
    yf_m = (yf-ymin)/(ymax-ymin) * sh
    x_m = (x-xmin)/(xmax-xmin) * sw
    y_m = (y-ymin)/(ymax-ymin) * sh


    # ---------- Depth estimation through ray tracing ------------ #
    # Take centroid ray from its sensor location to lens, forward project it incremently until it 'intersects' DSM 
    # The following diagram and variable naming is used.
    # First we compute d' from known f and computed xy' distance (distance of centroid to principal point in meters in sensor space location)
    # Then, we increment curr_H, compute d using similarity property, and use this d to project the point and find its XYZ location (3D image coordinates in image space)
    # We convert this to UTM, and use the UTM coordinates to index the DSM and obtain the actual H (vertical distance to camera) at the computed (XY) point
    # If this distance is close enough to the current estimate of the distance, we stop as the current d is close enough to the actual distance from the camera to the centroid.
    # Otherwise, repeat with a larger curr_H.
    #   This can be thought of as tracing the ray in the direction of the centroid until it intercepts the DSM, with the goal of finding the length of the ray that just meets the DSM.

    # +-xy'-+                    [sensor plane]
    #  \    | 
    # d'\   | f
    #    \  |
    #     \ |
    #      \|
    # ------+------              [lens]
    #       |\
    #       | \
    #       |  \
    #       |   \
    #curr_Zc|    \  d (depth)     NOTE: curr_Zc (depth of principal point) roughly equal to curr_H (height above point in 3D world space)
    #       |     \                     for nadir images. For oblique, need to compute curr_H from curr d and 3d world distance b/n camera center and polygon centroid
    #       |      \
    #       |       \
    #       +---xy---+           [plane at Z=curr_H in image local coords]

    # compute xy' (undist_lengths_m) and d' (d_dash) in meters
    undist_lengths_m = np.sqrt((x_m-xf_m[0])**2 + (y_m-yf_m[0])**2) 
    d_dash_m = np.sqrt(undist_lengths_m**2 + focal**2)


    # start 
    del_Z = coarse_del_Z
    curr_Zc = focal # starting estimate for Z (focal length in meters)
    while True:
        # increment curr_H and compute curr_D
        curr_Zc += del_Z
        curr_d = curr_Zc * d_dash_m[0] / focal # [0] refers to first point (centroid)

        # image local 3D coordinates for centroid (meters) : X = x*Z (Z:depth)
        Xc = -x[0]*curr_d  # [0] refers to first point (centroid)
        Yc =  y[0]*curr_d  

        # convert local to global coordinate space (UTM) using rotation matrix m 
        points = np.vstack([Xc,Yc,curr_Zc]).T
        curr_poly = points.dot(m.T)
        curr_poly = curr_poly[:,:2] / curr_poly[:,2:3]
        curr_poly = -curr_d * curr_poly + np.array([X0, Y0])

        # calculate elevation above point in 3D world space:
        dist_m_sq = (X0-curr_poly[0][0])**2 + (Y0-curr_poly[0][1]) # UTM distance between centroid and camera center
        curr_H = np.sqrt(curr_d**2 - dist_m_sq)

        # calculate pixel coordinate in dsm raster
        py, px = dsm.index(curr_poly[0][0], curr_poly[0][1])

        # get dsm value as median in a small bxb region centered at px and py and curr_H:
        slice_ = dsm_arr[py-bound:py+bound,px-bound:px+bound]
        dsm_value = np.median(slice_[slice_ > 0])
        if np.isnan(dsm_value):
            print(slice_, np.unique(slice_))
            exit("Error! No valid values in projected location. Please make sure DSM is for the correct location")
        
        # get distance to camera according to DSM at projected ray's location (curr_H') 
        # and compute difference with curr_H. If diff < tolerance, the current ray has succesfully intercepted the DSM 
        computed_del_H = Z0 - dsm_value # Z0: absolute altitude of drone over sea level. DSM also has absolute altitudes. Both in m
        diff = computed_del_H - curr_H

        if diff < 0:
            # error if we overshot the DSM (and landed up inside) - need to retry 
            print(f"Computed H => {computed_del_H}, Actual => {curr_H}")
            exit("Error! Too coarse. Reduce height increment del_H or increase tolerance, and try again")     
        elif diff < fine_tolerance:
            # successful intercept
            point_depth = curr_Zc * d_dash_m[0] / focal 
            # print(f"Found a good estimate of depth: {point_depth}")
            break
        elif diff < coarse_tolerance: 
            # reduce the height incrmeent (start finer approximation)
            del_Z = fine_del_H 

    
    # -------------- Actual Projection -------------- #
    # Now that we have the point depth of the centroid, we need to compute the depths of all actual points
    # using the estimated value of curr_H. xy*Z = XY (local coordinate space)
    # NOTE: alternatively, we can reuse the projection idea from above, but it may not be as robust and will be slower (by a bit). 

    corrected_point_depths = curr_Zc * d_dash_m / focal
    X = -x*corrected_point_depths
    Y =  y*corrected_point_depths

    # append focal length
    f_col = np.array([curr_Zc]*len(X))    # f should be curr_H here (in 3D image coordinate system, all are assumed to be on the same plane (but natrually at different depths from the camera))
    points = np.vstack([X,Y, f_col]).T

    # calculate new polygon points (in global coord system)
    curr_poly = points.dot(m.T)
    curr_poly = curr_poly[:,:2] / curr_poly[:,2:3]
    curr_poly = -corrected_point_depths.reshape(-1,1) * curr_poly + np.array([X0, Y0])

    # compute distance to image center (only used for batch shapefile creation)
    dist = np.sqrt((X0-curr_poly[0][0])**2 + (Y0-curr_poly[0][1])**2)

    # return
    return curr_poly, dist






def gis2img(img_path:str, 
            poly:np.array([np.uint16]), 
            dsm,
            dsm_arr,
            cfg:dict,
            ) -> np.array([np.uint16]): 
    """ 
    Projection from Image Coordinates to Georeferenced Shapefile 
    Inputs:
        - img_path: relative/absolute path to image to which we want to project
        - poly: numpy array of N UTM coordiantes. First point must be the centroid. Array is shaped (N,2), ordered clockwise.  
        - dsm_path: path to digital surface model (DSM) for canopy projection or digital terrain model (DTM) for ground point projection. 
                This dsm just needs to correspond to the area where the image footprint is expected to be.
        - config: dictionary of config file information. Contains the following main information.
            * calibration_dict: dict containing camera calibration info for the given image. 
                    Must contain keys <'f' in pixels, 'cx' in pixels, 'cy' in pixels, 'sw' (sensor width in meters), 'sh' (sensor_height in meters)>. 
                    Any of the following can be included for more precise results, but are assumed to be 0 if not present:
                        Radial distortion <'k1','k2','k3','k4'>
                        Tangential distortion <'p1','p2'>
                        Skew <'b1','b2'>
            * csv_path (recommended): path to tab-separated csv with headers of aligned images if available, e.g. from Agisoft. 
                    Expected headers: <PhotoID,X,Y,Z,Omega,Phi,Kappa,r11,r12,r13,r21,r22,r23,r31,r32,r33> 
                    Provides more precise results than relying on only exif headers for these values, but if not given, exif info is what the function uses. 

    Output:
        - img_points: numpy array of points projected to image pixel coordinates
        - exif_dict: for post processing
    """

    # extract required path & settings for processng $$
    cameras_csv_path = cfg['paths'].get('aligned_cameras_csv_path', None)

    bound = cfg['gis2img']['boundary']

    # extract calbiration information (e.g. given for DJI P1 in comments)
    calibration_dict = cfg['camera_calibration']
    f = calibration_dict['f']          # 11017.3
    cx = calibration_dict['cx']        # -6.35398
    cy = calibration_dict['cy']        # 22.0852
    sw = calibration_dict['sw']        # 36.0448 / 1000
    sh = calibration_dict['sh']        # 24.024 / 1000
    k1 = calibration_dict.get('k1', 0) # -0.00761618
    k2 = calibration_dict.get('k2', 0) # -0.245679 
    k3 = calibration_dict.get('k3', 0) # -0.520738
    k4 = calibration_dict.get('k4', 0) # 0.0
    p1 = calibration_dict.get('p1', 0) # -0.000425435
    p2 = calibration_dict.get('p2', 0) # 0.000371409
    b1 = calibration_dict.get('b1', 0) # -8.81372
    b2 = calibration_dict.get('b2', 0) # -0.708353
    
    # extract exif header
    exif_dict = extract_exif(img_path)

    
    # extract epsg
    try: dsm_epsg = int(str(dsm.crs)[5:])
    except: exit("Error: DSM does not seem to have CRS information. Chcek and try again")
    
    
    # get relevant information exif
    Z0 = float(exif_dict['Absolute Altitude']) 
    altitude = float(exif_dict['Relative Altitude'])              # increase -> box bigger
    image_height = float(exif_dict['Exif Image Height'])
    image_width = float(exif_dict['Exif Image Width'])
    focal = float(exif_dict['Focal Length'].split(' ')[0]) / 1000 # increase -> box bigger

    # GPS from LRF values directly + convert lat,long to UTM
    gps_lat = convert_gps(exif_dict["GPS Latitude"])
    gps_lon = convert_gps(exif_dict["GPS Longitude"])
    p = Proj(dsm_epsg, preserve_units=True)
    X0, Y0 = p(gps_lon, gps_lat)

    # rotation matrix m - NOTE: different drones/sensors have diffferent configurations for these, so may need to update this. 
    # TODO: fix this computation
    yaw = float(exif_dict['Gimbal Yaw Degree']) 
    roll = float(exif_dict['Gimbal Roll Degree']) 
    pitch = float(exif_dict['Gimbal Pitch Degree'])
    (kappa, omega, phi) = (yaw, pitch+90, roll)     # +90 needed due to DJI convention of gimbal vs camera view
    m = rot_matrix(kappa, omega, phi)

    # if csv path is given, get more precise estimates of m and Z0
    if cameras_csv_path is not None:
        df = pd.read_table(cameras_csv_path, skiprows=1)
        df.columns = ['PhotoID','X','Y','Z','Omega','Phi','Kappa','r11','r12','r13','r21','r22','r23','r31','r32','r33']

        img_name = img_path.split('/')[-1].split('.')[0] 
        row = df.loc[df['PhotoID'] == img_name]
        assert len(row) > 0, 'Error! Image not found in csv file.'

        omega = float(row['Omega'])
        phi = float(row['Phi'])
        kappa = float(row['Kappa'])
        m = rot_matrix(kappa, omega, phi)

        # X and Y used from exif, Z from csv
        # X0 = float(row['X'])
        # Y0 = float(row['Y'])
        Z0 = float(row['Z'])
    
    # ----------------- DEPTH ESTIMATION -------------------

    # determine curr_H for centroid (first value in poly) -- altitude difference between the point and camera in 3D world space
    py, px = dsm.index(poly[0][0], poly[0][1])
    slice_ = dsm_arr[py-bound:py+bound,px-bound:px+bound]
    centroid_dsm_value = np.median(slice_[slice_ > 0])
    if np.isnan(centroid_dsm_value):
            print(slice_, np.unique(slice_))
            exit("Error! No valid values in projected location. Please make sure DSM is for the correct location")
    curr_H = Z0 - centroid_dsm_value
    # print(curr_H)

    # get distance to center of camera for each point
    XY_dist_m_sq = (poly[:,0]-X0)**2 + (poly[:,1]-Y0)**2
    # print(poly, '-', np.array([X0,Y0]), '=',  '=', XY_dist_m_sq)

    # get depth of each point 
    point_depths_m = np.sqrt(curr_H**2 + XY_dist_m_sq)
    # print(point_depths_m, curr_H)


    # -------------------- PROJECTION ---------------------
    # convert UTM to 3d image local coordinate space

    poly = (poly - np.array([X0,Y0])) / (-point_depths_m).reshape(-1,1)
    poly = np.hstack([poly, np.ones((len(poly), 1))])
    pts = poly.dot(m.T)      # note: for m, inv == transpose NOTE: maybe remove .T
    # print(pts)

    # make intrinsic matrix and inverse
    k = np.array([[f+b1, b2,  cx+image_width//2],
                  [   0,  f, cy+image_height//2],
                  [   0,  0,                 1]])
    k_inv = np.linalg.inv(k)

    # principal point + four corners of image 
    full_extent = np.array([[image_width//2+cx,image_height//2+cy],[0,0],[image_width,0],[image_width,image_height],[0,image_height]])
    
    # utility for iterating for best estimate of x,y (image local coords at z=1)
    # used to get max entents in sensor space for full extent of image
    def iterate_undistort(points):
        points = points.reshape(-1,2)
        points = np.hstack([points, np.ones((len(points),1))])
        xy_prime = k_inv.dot(points.T)[:2,:]

        r = np.sqrt((xy_prime[0]**2 + xy_prime[1]**2))  # r_prime
        # x_0 and y_0 init -> initial undistorted image coordinates
        x = xy_prime[0] - xy_prime[0] * (k1 * r**2 + k2 * r**4 + k3 * r**6 + k4 * r**8) - p1 * (r**2 + 2 * xy_prime[0]**2) - 2 * p2 * xy_prime[0] * xy_prime[1] 
        y = xy_prime[1] - xy_prime[1] * (k1 * r**2 + k2 * r**4 + k3 * r**6 + k4 * r**8) - p2 * (r**2 + 2 * xy_prime[1]**2) - 2 * p1 * xy_prime[0] * xy_prime[1] 
        # iterate to get x,y -> final undistorted image coordinates
        for _ in range(10):
            r = np.sqrt((x**2 + y**2))
            x = xy_prime[0] - x * (k1 * r**2 + k2 * r**4 + k3 * r**6 + k4 * r**8) - p1 * (r**2 + 2 * x**2) - 2 * p2 * x * y
            y = xy_prime[1] - y * (k1 * r**2 + k2 * r**4 + k3 * r**6 + k4 * r**8) - p2 * (r**2 + 2 * y**2) - 2 * p1 * x * y

        return x,y
    
    def distort(x,y):
        r = np.sqrt((x**2 + y**2))
        x_prime = x + x*(k1 * r**2 + k2 * r**4 + k3 * r**6 + k4 * r**8) - p1 * (r**2 + 2 * x**2) - 2 * p2 * x * y
        y_prime = y + y*(k1 * r**2 + k2 * r**4 + k3 * r**6 + k4 * r**8) - p2 * (r**2 + 2 * y**2) - 2 * p1 * x * y
        print(x, x_prime)
        print(y, y_prime)

        print(x_prime.shape, y_prime.shape)

        xy_prime = np.vstack([x_prime, y_prime, np.ones((len(x_prime)))]).T # 5x3 


        # print(points.shape)
        # exit()

    
    xf,yf = iterate_undistort(full_extent)
    xmin, xmax = min(xf), max(xf)
    ymin, ymax = min(yf), max(yf)

    # 3d image local space (XYZ) --> sensor space location in m
    # print(distort(pts[:,0], pts[:,1]))
    pts = pts*-focal   # backproject 
    pts = pts[:,:2]    # remove last col (-focal)
    pts = pts / np.array([sw, sh])
    # print(pts), exit()
    # print(pts)
    pts *= np.array([image_width+cx, -(image_height+cy)])
    pts += np.array([image_width/2, image_height/2]) 
    pts = pts.astype(np.uint64)

    # sensor space location --> undistorted points (u'v')
    # first convert to correct sensor space location by unscaling between limits
    # x = (pts[:,0]/sw) * (xmax-xmin) + xmin
    # y = (pts[:,1]/sh) * (ymax-ymin) + ymin

    # undist_pts = pts/np.array([sw, sh]) * np.array([(xmax-xmin),(ymax-ymin)]) + np.array([xmin,ymin])
    

    # undist_pts *= np.array([-(image_width+cx), -(image_height+cy)])
    # undist_pts += np.array([image_width/2, image_height/2]) 
    # undist_pts = undist_pts.astype(np.uint64)

    return pts, exif_dict



    # ------------------- DISTORTION -------------------------
    # undistorted points (u'v') --> distorted points (uv) [output]



def get_candidate_img_list(extents, corners):
    ''' returns list of idx in extents where corners is entirely contained within extents'''
    
    def contained(lbox, sbox):
        '''returns True if sbox contained in lbox'''
        s_xmin, s_ymin = np.amin(sbox, axis=0)
        s_xmax, s_ymax = np.amax(sbox, axis=0)

        l_xmin, l_ymin = np.amin(lbox, axis=0)
        l_xmax, l_ymax = np.amax(lbox, axis=0)

        return (s_xmin >= l_xmin and s_ymin >= l_ymin and s_xmax <= l_xmax and s_ymax <= l_ymax)
        
    valid = []
    for i in range(len(extents)):
        if contained(extents[i],corners):
            valid.append(i)
    return valid
        

def get_extent_of_all_images(img_paths, dsm, dsm_arr, cfg):
    ''' assuming 0 height for projection '''
    extents = []
    for path in tqdm(img_paths):
        extents.append(img2gis(path, None, dsm, dsm_arr, cfg, full_image=True)[0])
    return extents




if __name__ == "__main__":
    pass

    
