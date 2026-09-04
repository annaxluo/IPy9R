# process output from CellProfiler 
import os

import numpy as np
import pandas as pd
import tifffile


def make_tiles(img_fn, out_path, out_fn_str, 
    round_label, ch_label, used_stacks, zplane_labels, 
    max_h=2000, max_w=2000): 
    """create tiles from image for valid input to cellprofiler""" 
    img = tifffile.imread(img_fn)

    (iZ, iH, iW) = img.shape
    n_tiles_h = int(iH / max_h) if (iH % max_h == 0) else int(iH / max_h) + 1
    n_tiles_w = int(iW / max_w) if (iW % max_w == 0) else int(iW / max_w) + 1
    
    # get size in micron
    tif_tags = {}
    with tifffile.TiffFile(img_fn) as tif: 
        for tag in tif.pages[0].tags.values(): 
            name, val = tag.name, tag.value
            tif_tags[name] = val
    x_res = tif_tags["XResolution"]
    y_res = tif_tags["YResolution"]
    scale_factor_h = y_res[1] / y_res[0]
    scale_factor_w = x_res[1] / x_res[0]
    
    # multiple z-planes
    for (z_i, zp_lab) in zip(used_stacks, zplane_labels): 
        tiles = []
        metadata = pd.DataFrame(columns=["tile", "h_start_px", "h_end_px", 
            "w_start_px", "w_end_px", "h_start", "h_end", 
            "w_start", "w_end"])
    
        for h in range(n_tiles_h): 
            for w in range(n_tiles_w): 
                h_start_px = h * max_h 
                h_end_px = min((h + 1) * max_h, iH)
                w_start_px = w * max_w
                w_end_px = min((w + 1) * max_w, iW)
                sub = img[z_i, h_start_px:h_end_px, w_start_px:w_end_px]
                tiles.append(sub)
                # in micron 
                h_start = h_start_px * scale_factor_h
                h_end = h_end_px * scale_factor_h
                w_start = w_start_px * scale_factor_w
                w_end = w_end_px * scale_factor_w
                metadata.loc[len(metadata.index)] = [len(tiles)-1, 
                                                     h_start_px, h_end_px, 
                                                     w_start_px, w_end_px, 
                                                     h_start, h_end, 
                                                     w_start, w_end] 
                                                     
        # write outputs     
        mdata_out_fn = os.path.join(out_path, 
            f"{out_fn_str}-r{round_label}-c{ch_label}-z{zp_lab}.csv")
        metadata.to_csv(mdata_out_fn, sep=',', header=True, index=False)
        
        for i in range(len(tiles)): 
            out_f = os.path.join(out_path, 
                                 f"{out_fn_str}-f{i}-r{round_label}-c{ch_label}-z{zp_lab}.tiff")
            tifffile.imwrite(out_f, tiles[i])
            print("written:", out_fn_str + "_tile_" + str(i).zfill(2) + ".tif")
        
    return len(tiles) # return number of tiles per z-plane
   
# # change this: 
# def combine_tiles(coordinates, tiles_list, 
    # x_str = "w_start_px", 
    # y_str = "h_start_px", 
    # t_str = "tile"): 
    # """combine tiles in tiles into one stitched image""" 
    # return True
        
def make_starfish_data(fov_per_tile, tile_info_fn, primary_dir, nuclei_dir, z_step): 
    """makes a structured data object for the Starfish pipeline""" 
    import csv
    import shutil
    
    # tile info
    tile_data = pd.read_csv(tile_info_fn) 
    n_tiles = tile_data.shape[0]
    # make fovs
    fovs = [fov_per_tile] * n_tiles
    # coordinates: zc values are arbitrary
    coordinates_of_fovs = []
    tile_data = tile_data.reset_index()
    for index, row in tile_data.iterrows(): 
        coords = {
            'xc_min': row['w_start'],
            'xc_max': row['w_end'],
            'yc_min': row['h_start'],
            'yc_max': row['h_end'],
            'zc_min': z_step * index,
            'zc_max': z_step * (index + 1) 
            }
        coordinates_of_fovs.append(coords)
    
    # write coordinates file for primary and nuclei in their respective directories
    with open(os.path.join(primary_dir, "coordinates.csv"), "w") as fh:
        csv_writer = csv.DictWriter(
            fh,
            [
                'fov', 'round', 'ch', 'zplane',
                'xc_min', 'yc_min', 'zc_min', 'xc_max', 'yc_max', 'zc_max',
            ]
        )
        csv_writer.writeheader()
        
        for fov_id, (fov_info, coordinate_of_fov) in enumerate(zip(fovs, coordinates_of_fovs)):
            for round_label, ch_label, zplane_label in fov_info:
                tile_coordinates = coordinate_of_fov.copy()
                tile_coordinates.update({
                    'fov': fov_id,
                    'round': round_label,
                    'ch': ch_label,
                    'zplane': zplane_label,
                })
                csv_writer.writerow(tile_coordinates)

    shutil.copyfile(
        os.path.join(primary_dir, "coordinates.csv"), os.path.join(nuclei_dir, "coordinates.csv"))
        
def convert_to_spacetx(primary_dir, primary_out, nuclei_dir, nuclei_out): 
    """converts structured Starfish data to SpaceTx format""" 
    from slicedimage import ImageFormat
    from starfish.experiment.builder import format_structured_dataset
    
    format_structured_dataset(
        primary_dir,
        os.path.join(primary_dir, "coordinates.csv"),
        primary_out,
        ImageFormat.TIFF,
    )
    
    format_structured_dataset(
        nuclei_dir,
        os.path.join(nuclei_dir, "coordinates.csv"),
        nuclei_out,
        ImageFormat.TIFF,
    )
    
    # modify experiment.json
    with open(os.path.join(primary_out, "experiment.json"), "r+") as fh:
        contents = fh.readlines()
        print("original experiment.json\n")
        print("".join(contents))
        contents[3] = ",".join([contents[3].strip("\n"),"\n"])
        contents.insert(4, '\t"nuclei": "../nuclei/nuclei.json"\n')  # new_string should end in a newline
        fh.seek(0)  # readlines consumes the iterator, so we need to start over
        fh.writelines(contents)  # No need to truncate as we are increasing filesize
        fh.seek(0)
        print("\nmodified experiment.json\n")
        print(fh.read())
        
        
def modify_codebook(primary_out, mappings_list): 
    """modify the codebook""" 
    import json
    
    with open(os.path.join(primary_out, "codebook.json"), "r") as fh:
        codebook_old = json.load(fh)

    codebook = {
        "version": codebook_old['version'], 
        "mappings": mappings_list
        }
        
    with open(os.path.join(primary_out, "codebook.json"), "w") as fh:
        json.dump(codebook, fh)


def spots_count(img_fn, mask_fn_list, feature_table_fn, contour_thickness=3): 
    """extract puncta within a mask"""
    from skimage import measure
    from skimage.draw import circle_perimeter, disk
    
    img = tifffile.imread(img_fn)
    
    spots_data = pd.read_csv(feature_table_fn)
    spots_in_mask = pd.DataFrame(columns=["X", "Y", "radius", "mask_id"])
    img_out = np.dstack([img, img, img])
    for mask_fn in mask_fn_list: 
        mask_id = os.path.split(mask_fn)[1].split("_")[0]
        mask = tifffile.imread(mask_fn)

        # threshold mask 
        mask = (mask > 0)
        # visualize barrel contours
        contours = measure.find_contours(mask)
        mask_shape = mask.shape
        rc = [disk(c, contour_thickness, shape=mask_shape) for c in contours[0]]
        for rc_ in rc: 
            img_out[rc_[0], rc_[1], :] = [0, 255, 0]
        
        # find selection points within the mask 
        in_mask_ = pd.DataFrame(columns=["X", "Y", "radius", "mask_id"])
        for x, y, rr in zip(spots_data["x_original"], spots_data["y_original"], spots_data["radius"]): 
            if mask[y, x]: 
                in_mask_.loc[len(in_mask_.index)] = [x, y, rr, mask_id]
                
        # draw selections 
        for x_, y_, r_ in zip(in_mask_["X"], in_mask_["Y"], in_mask_["radius"]):
            row, col = circle_perimeter(y_, x_, int(r_), shape=mask_shape)
            img_out[row, col, :] = [255, 0, 0]
                    
        spots_in_mask = pd.concat([spots_in_mask, in_mask_])
        
        print("processed mask", mask_id)
        
    return img_out, spots_in_mask
    
    
def spots_count_stack(img_fn, mask_fn_list, feature_table_fn, contour_thickness=3): 
    """count spots from a stack
       - img_fn: tif, a stack 
    """
    from skimage import measure
    from skimage.draw import circle_perimeter, disk
    
    img = tifffile.imread(img_fn)
    
    spots_data = pd.read_csv(feature_table_fn)
    spots_in_mask = pd.DataFrame(columns=["X", "Y", "Z", "intensity", "radius", "mask_id"])
    img_out = np.zeros(img.shape, dtype="uint8")
    for mask_fn in mask_fn_list: 
        mask_id = os.path.split(mask_fn)[1].split("_")[0]
        mask = tifffile.imread(mask_fn)

        # threshold mask 
        mask = (mask > 0)
        # visualize barrel contours
        contours = measure.find_contours(mask)
        mask_shape = mask.shape
        rc = [disk(c, contour_thickness, shape=mask_shape) for c in contours[0]]
        for rc_ in rc: 
            img_out[:, rc_[0], rc_[1]] = 255
        
        # find selection points within the mask 
        in_mask_ = pd.DataFrame(columns=["X", "Y", "Z", "intensity", "radius", "mask_id"])
        for (x, y, z, rr, intensity) in zip(spots_data["x_original"], spots_data["y_original"], 
            spots_data["z"], spots_data["radius"], spots_data["intensity"]): 
            if mask[y, x]: 
                in_mask_.loc[len(in_mask_.index)] = [x, y, z, intensity, rr, mask_id]
                
        # draw selections 
        for (x_, y_, z_, r_) in zip(in_mask_["X"], in_mask_["Y"], in_mask_["Z"], in_mask_["radius"]):
            row, col = circle_perimeter(y_, x_, int(r_), shape=mask_shape)
            img_out[z_, row, col] = 255
                    
        spots_in_mask = pd.concat([spots_in_mask, in_mask_])
        
        print("processed mask", mask_id)
        
    return img_out, spots_in_mask
    
    
def summarize_count(spots_in_mask, groups, mask_fn_list): 
    """create summary of spots count stats"""
    # list of mask ids
    mask_ids = []
    for mask_fn in mask_fn_list: 
        mask_ids.append(os.path.split(mask_fn)[1].split("_")[0])
    
    summary = pd.DataFrame(columns=["mask_id", "spot_count", "group"])
    
    for mask_id in mask_ids: 
        count = sum(spots_in_mask["mask_id"] == mask_id)
        m_group = [grp for grp, m_id in groups.items() if mask_id in m_id]
        summary.loc[len(summary.index)] = [mask_id, count, m_group[0]]
        
    return summary
    

def summarize_count_stack(spots_in_mask, groups, mask_fn_list): 
    """create summary of spots count stats"""
    # list of mask ids
    mask_ids = []
    for mask_fn in mask_fn_list: 
        mask_ids.append(os.path.split(mask_fn)[1].split("_")[0])
    
    summary = pd.DataFrame(columns=["mask_id", "group", "z", "spot_count", 
        "mean_intensity", "std_intensity", "mean_radius", "std_radius"])
    
    z_idx = np.unique(spots_in_mask['Z'])
    
    for mask_id in mask_ids: 
        for z_ in z_idx: 
            valid_idx = (spots_in_mask["mask_id"] == mask_id) & (spots_in_mask["Z"] == z_)
            count = sum(valid_idx)
            m_group = [grp for grp, m_id in groups.items() if mask_id in m_id]
            summary.loc[len(summary.index)] = [mask_id, m_group[0], z_, count, 
                np.mean(spots_in_mask["intensity"][valid_idx]), 
                np.std(spots_in_mask["intensity"][valid_idx]), 
                np.mean(spots_in_mask["radius"][valid_idx]), 
                np.std(spots_in_mask["radius"][valid_idx])]
        
    return summary


def spot_density(summary, img_fn, mask_fn_list): 
    """compute density of cells within the masked region""" 
    # get resolution information from tif file 
    tif_tags = {}
    with tifffile.TiffFile(img_fn) as tif: 
        for tag in tif.pages[0].tags.values(): 
            name, val = tag.name, tag.value
            tif_tags[name] = val
    if len(tif_tags) == 0: 
        print("cannot open tiff file.")
        return None
        
    x_res = tif_tags["XResolution"]
    y_res = tif_tags["YResolution"]
    resolution_factor = (x_res[1]/x_res[0])*(y_res[1]/y_res[0]) # micron2 per pixel2
    
    summary_out = summary.copy()
    summary_out["Area"] = ""
    summary_out["Density"] = ""
    for mask_fn in mask_fn_list: 
        mask_id = os.path.split(mask_fn)[1].split("_")[0]
        row_id = summary_out.loc[summary_out['mask_id'] == mask_id].index[0]
        
        mask = tifffile.imread(mask_fn)
        # threshold mask 
        mask = (mask > 0)
        # compute mask area
        m_area = np.count_nonzero(mask) * resolution_factor # in micron^2
        m_density = summary_out.loc[row_id, ["spot_count"]][0] / m_area
        summary_out.loc[row_id, ["Area", "Density"]] = [m_area, m_density]
    
    return summary_out
    
    
def spot_density_stack(summary, img_fn, mask_fn_list): 
    """compute density of cells within the masked region""" 
    # get resolution information from tif file 
    tif_tags = {}
    with tifffile.TiffFile(img_fn) as tif: 
        for tag in tif.pages[0].tags.values(): 
            name, val = tag.name, tag.value
            tif_tags[name] = val
    if len(tif_tags) == 0: 
        print("cannot open tiff file.")
        return None
        
    x_res = tif_tags["XResolution"]
    y_res = tif_tags["YResolution"]
    resolution_factor = (x_res[1]/x_res[0])*(y_res[1]/y_res[0]) # micron2 per pixel2
    
    summary_out = summary.copy()
    summary_out["Area"] = ""
    summary_out["Density"] = ""
    for mask_fn in mask_fn_list: 
        mask_id = os.path.split(mask_fn)[1].split("_")[0]
        row_id = summary_out.loc[summary_out['mask_id'] == mask_id].index
        
        mask = tifffile.imread(mask_fn)
        # threshold mask 
        mask = (mask > 0)
        # compute mask area
        m_area = np.count_nonzero(mask) * resolution_factor # in micron^2
        m_density = summary_out.loc[row_id, "spot_count"] / m_area
        summary_out.loc[row_id, "Area"] = m_area
        summary_out.loc[row_id, "Density"] = m_density
    
    return summary_out
