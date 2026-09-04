# process output from CellProfiler 
import os

import cv2
import numpy as np
import pandas as pd
import tifffile


def find_positive_nuclei(out_path, data_fn, gene_name, use_tiles=False, 
    fn_str=None, img_nuclei_fn=None, img_gene_fn=None, 
    threshold=.2, annotation_radius=3, annotation_thickness=5): 
    """identifies nuclei positive for an RNA"""
    var_name = "Intensity_MeanIntensity_" + gene_name + "_enhanced"
    
    data = pd.read_csv(data_fn)
    data2 = data[data[var_name] > threshold]
    valid_loc = data2[["Location_CenterMassIntensity_X_Nuclei", 
                       "Location_CenterMassIntensity_Y_Nuclei", 
                       "ImageNumber"]]

    if (img_nuclei_fn is not None) and (not use_tiles): 
        img_nuclei = cv2.imread(img_nuclei_fn, 0)
        img_nuclei_out = np.dstack([img_nuclei, img_nuclei, img_nuclei])
        for x, y in zip(valid_loc.iloc[:,0], valid_loc.iloc[:,1]): 
            cv2.circle(
                img_nuclei_out, 
                (int(x), int(y)), 
                annotation_radius, 
                (255, 0, 0), 
                thickness=annotation_thickness, 
                lineType=cv2.LINE_AA)
        out_fn = os.path.join(out_path, 
            fn_str + "selection_DAPI_thresh" + str(int(threshold * 100)) + ".tif")
        tifffile.imwrite(out_fn, img_nuclei_out)
    
    if (img_gene_fn is not None) and (not use_tiles): 
        img_gene = cv2.imread(img_gene_fn, 0)
        img_gene_out = np.dstack([img_gene, img_gene, img_gene])
        for x, y in zip(valid_loc.iloc[:,0], valid_loc.iloc[:,1]): 
            cv2.circle(
                img_gene_out, 
                (int(x), int(y)), 
                annotation_radius, 
                (255, 0, 0), 
                thickness=annotation_thickness, 
                lineType=cv2.LINE_AA)
        out_fn2 = os.path.join(out_path, 
            fn_str + "selection_" + gene_name + "_thresh" + str(int(threshold * 100)) + ".tif")
        tifffile.imwrite(out_fn2, img_gene_out)

    # write csv for imagej import 
    if use_tiles: 
        ff = "_tiles.csv"
    else: 
        ff = ".csv"
        
    out_fn3 = os.path.join(out_path, 
        fn_str + "selection_pts_thresh" + str(int(threshold * 100)) + ff)
    valid_loc.to_csv(out_fn3, sep=',', header=True, index=True)
    
    
def make_tiles(img_fn, out_path, max_h=2000, max_w=2000): 
    """create tiles from image for valid input to cellprofiler""" 
    img = cv2.imread(img_fn, 0)

    (iH, iW) = img.shape
    
    n_tiles_h = int(iH / max_h) if (iH % max_h == 0) else int(iH / max_h) + 1
    n_tiles_w = int(iW / max_w) if (iW % max_w == 0) else int(iW / max_w) + 1
    
    tiles = []
    metadata = pd.DataFrame(columns=["tile", "h_start", "h_end", "w_start", "w_end"])
    for h in range(n_tiles_h): 
        for w in range(n_tiles_w): 
            h_start = h * max_h 
            h_end = min((h + 1) * max_h, iH)
            w_start = w * max_w
            w_end = min((w + 1) * max_w, iW)
            sub = img[h_start:h_end, w_start:w_end]
            tiles.append(sub)
            metadata.loc[len(metadata.index)] = [len(tiles)-1, h_start, h_end, 
                                                 w_start, w_end] 
                                                 
    # write outputs 
    fn_part = os.path.split(img_fn)[-1].split(".")[0]
    
    mdata_out_fn = os.path.join(out_path, fn_part + "_tiles.csv")
    metadata.to_csv(mdata_out_fn, sep=',', header=True, index=False)
    
    for i in range(len(tiles)): 
        out_f = os.path.join(out_path, 
                             fn_part + "_tile_" + str(i).zfill(2) + ".tif")
        tifffile.imwrite(out_f, tiles[i])
        print("written:", fn_part + "_tile_" + str(i).zfill(2) + ".tif")
        
    return len(tiles)
   

def combine_tiles(tile_data_fn, selection_fn, output_fn, 
    x_str = "Location_CenterMassIntensity_X_Nuclei", 
    y_str = "Location_CenterMassIntensity_Y_Nuclei", 
    t_str = "ImageNumber"): 
    """combine points in tiles into one stitched image""" 
    tile_data = pd.read_csv(tile_data_fn) 
    sel_pts = pd.DataFrame(columns=["idx", x_str, y_str, t_str])
    
    # selections
    sel = pd.read_csv(selection_fn)
    tile_id = np.unique(sel[t_str]) - 1
    tile_id = tile_id.tolist()
    
    try: 
        for ii, t, x, y in zip(sel.iloc[:,0], sel[t_str], sel[x_str], sel[y_str]): 
            x1 = x + tile_data.loc[t-1].w_start
            y1 = y + tile_data.loc[t-1].h_start
            sel_pts.loc[len(sel_pts.index)] = [ii, x1, y1, t]
    except pd.errors.EmptyDataError: 
        print("File contains no selection points.")
        
    # save data 
    sel_pts.to_csv(output_fn, sep=',', header=True, index=False)
    
    return True
    
