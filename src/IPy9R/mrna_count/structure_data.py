# src/IPy9R/mrna_count/structure_data.py

# structure image tiles for starfish input
import argparse
import csv
import json
import os
import shutil

import numpy as np
import tifffile
from slicedimage import ImageFormat
from starfish import Experiment
from starfish.experiment.builder import format_structured_dataset


def run(data_path_base,
        roi_list, 
        used_channels,
        channel_gene_names, 
        used_z_list,         
        selection_rect, 
        image_id_list, 
        z_thickness, 
        roi_folder_prefix="confocal_63x_"
        ) -> None: 
    """
    Structure image tiles for starfish input.  
    """
    # validate input lengths 
    n = len(roi_list)

    if not (len(image_id_list) == len(used_z_list) == len(selection_rect) == n):
        raise ValueError(
            "roi_list, image_id_list, used_z_list, and selection_rect must have the same length."
        )

    # 1. prepare paths -----------------------------------------------------
    # output path
    output_path_base_list = [os.path.join(data_path_base, roi_folder_prefix + s)
        for s in roi_list]
    output_path = os.path.join(data_path_base, "SF_analysis") 
    
    # 2. define fovs and coordinates ----------------------------------------
    # fovs: only 1 round, 1 channel, and 1 z-plane
    fovs = [[(0, ch, 0) for ch in range(len(used_channels))]] * len(roi_list)

    coordinates_of_fovs = []

    for (roi_d, im_, z_, rect_) in zip(output_path_base_list, image_id_list, used_z_list, selection_rect): 
        img_fn = os.path.join(roi_d, "CP_inputs", 
            used_channels[0] + "-" + im_ + ("-000" + str(z_) if z_ is not None else "") + ".tif")
        with tifffile.TiffFile(img_fn) as tif: 
            x_res = tif.pages[0].tags['XResolution'].value
            w_ = tif.pages[0].tags['ImageWidth'].value if rect_ is None else rect_[2]
            w_micron = w_ / (x_res[0] / x_res[1])
            y_res = tif.pages[0].tags['YResolution'].value
            h_ = tif.pages[0].tags['ImageLength'].value if rect_ is None else rect_[3]
            h_micron = h_ / (y_res[0] / y_res[1])
        coords_ = {
            'xc_min': 0.0, 
            'xc_max': w_micron, 
            'yc_min': 0.0, 
            'yc_max': h_micron, 
            'zc_min': 0, 
            'zc_max': z_thickness 
            }
        coordinates_of_fovs.append(coords_)

    # 3. Update filenames -------------------------------------------------
    primary_dir = os.path.join(output_path, "primary_dir") 
    os.makedirs(primary_dir, exist_ok=True)
            
    # copy and rename files  
    original_image_size_list = [] # original image size
    selection_ul_micron = [] # position of selection rect UL in micron
    for fov_id, (fov_info, roi_d, im_, z_, rect_) in enumerate(zip(fovs, output_path_base_list, 
        image_id_list, used_z_list, selection_rect)): 
        # source image fns
        src_fns = [os.path.join(roi_d, "CP_inputs", 
            ch + "-" + im_ + ("-000" + str(z_) if z_ is not None else "") + ".tif") 
            for ch in used_channels] 
        # destination image fns
        des_fns = []
        for round_label, ch_label, zplane_label in fov_info: 
            fn_ = os.path.join(primary_dir, 
                f"primary-f{fov_id}-r{round_label}-c{ch_label}-z{zplane_label}.tiff")
            des_fns.append(fn_)
        for (src_, des_) in zip(src_fns, des_fns): 
            if rect_ is not None:
                src_img = tifffile.imread(src_)
                des_img = src_img[rect_[1]:rect_[1]+rect_[3], rect_[0]:rect_[0]+rect_[2]]
                tifffile.imwrite(des_, des_img)
            else: 
                shutil.copyfile(src_, des_)
        # original image size
        with tifffile.TiffFile(src_fns[0]) as tif: 
            w_ = tif.pages[0].tags['ImageWidth'].value 
            h_ = tif.pages[0].tags['ImageLength'].value 
            original_image_size_list.append([h_, w_]) # row, col            
        
            # selection rectangle (x,y) in micron
            x_res = tif.pages[0].tags['XResolution'].value
            x_ = 0 if rect_ is None else rect_[0]
            x_micron = x_ / (x_res[0] / x_res[1])
            y_res = tif.pages[0].tags['YResolution'].value
            y_ = 0 if rect_ is None else rect_[1]
            y_micron = y_ / (y_res[0] / y_res[1])
            selection_ul_micron.append([x_micron, y_micron])
    

    # write coordinate files 
    coordinates_fn = os.path.join(primary_dir, "coordinates.csv") 
    with open(coordinates_fn, "w") as fh: 
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


    # 4. Convert structured data into SpaceTx formmat------------------------------------
    primary_out = os.path.join(output_path, "primary")  
    os.makedirs(primary_out, exist_ok=True)

    format_structured_dataset(
        primary_dir,
        os.path.join(primary_dir, "coordinates.csv"),
        primary_out,
        ImageFormat.TIFF,
    )

    # updata codebook
    cb_fn = os.path.join(primary_out, "codebook.json")

    mappings_list = []
    for (fov_, gene_name) in zip(fovs[0], channel_gene_names): 
        mapping_ = {
            "codeword": [
                {"c": fov_[1], "r": fov_[0], "v": 1.0}
            ], 
            "target": gene_name}
        mappings_list.append(mapping_)
        
    codebook_new = {
        "version": '0.0.0', 
        "mappings": mappings_list
        }
    
    with open(cb_fn, "w") as fh:
        json.dump(codebook_new, fh)

    # 5. write data structure records -------------------------------------
    exp = Experiment.from_json(os.path.join(primary_out, "experiment.json"))

    fov_names = [fov.name for fov in exp.fovs()]

    fov_info = {k:None for k in fov_names}
    for ii in range(len(fov_names)): 
        if used_z_list[ii] is not None:
            channel_images = {k:f"{used_channels[k]}-{image_id_list[ii]}-000{used_z_list[ii]}.tif" 
                for k in range(len(used_channels))}
        else: 
            channel_images = {k:f"{used_channels[k]}-{image_id_list[ii]}.tif" 
                for k in range(len(used_channels))}
                
        ret_ = {
            "fov_name": fov_names[ii], 
            "used_roi": roi_list[ii], 
            "channels": {k:v for k, v in zip(range(len(used_channels)), used_channels)}, 
            "channel_genes": {k:v for k, v in zip(range(len(channel_gene_names)), channel_gene_names)}, 
            "used_z_stack": used_z_list[ii], 
            "original_image_size": original_image_size_list[ii], 
            "selection_rect": selection_rect[ii], 
            "selection_ul_micron": selection_ul_micron[ii], 
            "channel_images": channel_images
            }
        fov_info[fov_names[ii]] = ret_
            
    data_struct = {
        "data_path_base": data_path_base, 
        "experiment_json": "primary/experiment.json", 
        "primary_dir": "primary_dir", 
        "primary_out": "primary", 
        "num_fovs": len(fov_names), 
        "num_channels": len(used_channels), 
        "num_rois": np.unique(roi_list).shape[0], 
        "fov_info": fov_info}
        
    data_struct_fn = os.path.join(output_path, "data_structure.json")    
    with open(data_struct_fn, "w") as out_f:
        json.dump(data_struct, out_f, indent=4)
        
        
def main(): 
    parser = argparse.ArgumentParser(description="Make data structure for Starfish.")
    parser.add_argument("--config", required=True)

    args = parser.parse_args()

    from IPy9R.config import get_step_config, load_config

    config = load_config(args.config)
    step_config = get_step_config(config, "mrna_count")

    run(
        data_path_base=step_config["data_path_base"],
        roi_list=step_config["roi_list"], 
        used_channels=step_config["used_channels"], 
        channel_gene_names=step_config["channel_gene_names"],  
        used_z_list=step_config["used_z_list"], 
        selection_rect=step_config["selection_rect"], 
        image_id_list=step_config["image_id_list"], 
        z_thickness=step_config["z_thickness"], 
        roi_folder_prefix=step_config.get("roi_folder_prefix", "confocal_63x_")
    )
    
if __name__ == "__main__":
    main()