# src/IPy9R/mrna_count/detect_mrna_spots.py

# preprocess image tiles for starfish processing: intensity normalization, 
# scaling and clipping; followed by detection of mRNA puncta (spots) from image tiles 
import argparse
import json
import os

import numpy as np
import tifffile
from skimage.draw import circle_perimeter
from starfish import Experiment
from starfish.spots import DecodeSpots, FindSpots
from starfish.types import TraceBuildingStrategies

from IPy9R.mrna_count.quantification_utils import otsu_threshold, preprocess_image


def run(data_path_base,
        roi_list, 
        used_channels,
        channel_gene_names, 
        min_sigma=3,
        max_sigma=10,
        num_sigma=1,
        search_radius=10,
        otsu_n_classes=3
        ) -> None: 
    """
    Preprocess image tiles for starfish processing: intensity normalization, 
    scaling and clipping; followed by detection of mRNA puncta (spots) from image tiles 
    """
    # input validation
    if len(channel_gene_names) != len(used_channels):
        raise ValueError("channel_gene_names and used_channels must have the same length")

    # 1. load data -----------------------------------------------
    input_path = os.path.join(data_path_base, "SF_analysis") 
    experiment = Experiment.from_json(os.path.join(input_path, "primary", "experiment.json"))
    #print(experiment.fovs())

    output_path_base = os.path.join(input_path, "SF_outputs")
    os.makedirs(output_path_base, exist_ok=True)
        
    preprocess_output_path = os.path.join(output_path_base, "preprocess_outputs")
    os.makedirs(preprocess_output_path, exist_ok=True)

    threshold_output_path = os.path.join(output_path_base, "threshold_outputs")
    os.makedirs(threshold_output_path, exist_ok=True)
       
    detection_output_path = os.path.join(output_path_base, "detection_outputs")
    os.makedirs(detection_output_path, exist_ok=True)
        
    # 2. preprocess and detect puncta from each image -----------------------
    # fov information
    data_struct_fn = os.path.join(input_path, "data_structure.json")
    with open(data_struct_fn, "r") as f: 
        data_struct = json.load(f)

    for fov in experiment.fovs(): 
        fov_name = fov.name
        img = fov.get_image("primary")
        
        # validate channels
        if img.num_chs != len(channel_gene_names):
            raise ValueError(
                f"{fov_name}: Image has {img.num_chs} channels, "
                f"but channel_gene_names has {len(channel_gene_names)} entries."
            )
        
        # preprocess image
        print("preprocessing image:" + fov_name)
        (processed_img, saved_fns_dict) = preprocess_image(img, 
            save_filtered=True, save_norm=True, 
            output_path=preprocess_output_path, fov_id=fov_name)

        # get otsu thresholds
        print("computing otsu thresholds")
        blob_thresholds = otsu_threshold(processed_img, n_classes=otsu_n_classes, 
            save_thresholded=True, output_path=threshold_output_path, fov_id=fov_name)

        # detect puncta ------------------------
        #thr_keys = [f"threshold_{i}" for i in range(len(blob_thresholds))]
        #detected_outputs_ = [None for x in range(len(blob_thresholds))]
        ch_keys = [f"channel_{j}" for j in range(img.num_chs)]
        detected_outputs = {f"channel_{ch}":{f"threshold_{th}":None for th in range(len(blob_thresholds))} 
            for ch in range(img.num_chs)}
        
        for ii in range(len(blob_thresholds)):     
            thr = blob_thresholds[ii]
            print("detecting puncta using threshold", ii)
            bd = FindSpots.BlobDetector(
                min_sigma=min_sigma,
                max_sigma=max_sigma,
                num_sigma=num_sigma, 
                threshold=thr,
                measurement_type='mean',
                exclude_border=False, 
                detector_method='blob_log', 
                is_volume=False, 
            )
            bd_spots = bd.run(image_stack=processed_img)
            decoder = DecodeSpots.PerRoundMaxChannel(
                codebook=experiment.codebook,
                anchor_round=0,
                search_radius=search_radius,
                trace_building_strategy=TraceBuildingStrategies.NEAREST_NEIGHBOR)
            bd_decoded = decoder.run(spots=bd_spots)
            decode_mask = bd_decoded['target'] != 'nan'
            bd_feats = decode_mask.to_features_dataframe()
            bd_feats['threshold'] = thr
            #bd_table = build_spot_traces_exact_match(bd_spots)
            #bd_feats = bd_table.to_features_dataframe() 
            # visualize ---------------------------------------------
            # parameters for selection rect adjustment        
            rect_ = data_struct['fov_info'][fov_name]['selection_rect']
            if rect_ is not None: 
                mask_shape = data_struct['fov_info'][fov_name]['original_image_size']
                sel_ul_micron = data_struct['fov_info'][fov_name]['selection_ul_micron']
                # adjust feature table
                bd_feats['x'] = bd_feats['x'] + rect_[0]
                bd_feats['xc'] = bd_feats['xc'] + sel_ul_micron[0]
                bd_feats['y'] = bd_feats['y'] + rect_[1]
                bd_feats['yc'] = bd_feats['yc'] + sel_ul_micron[1]     
            else: 
                mask_shape = processed_img.tile_shape

        
            for ch in range(img.num_chs): 
                target_ = channel_gene_names[ch]
                mask_out = np.zeros(mask_shape, dtype=np.uint8) 
                # plot spots
                bd_feats_sub = bd_feats.loc[bd_feats['target']==target_]
                for (xx, yy, rr) in zip(bd_feats_sub['x'], bd_feats_sub['y'], bd_feats_sub['radius']):
                    row, col = circle_perimeter(yy, xx, int(rr), shape=mask_shape)
                    mask_out[row, col] = 1
                
                out_fn =  os.path.join(detection_output_path, f"thresh{ii}-{fov_name}-c{ch}-r0-z0.tiff")
                tifffile.imwrite(out_fn, mask_out)
                # save the feature table 
                bd_feats_fn = os.path.join(detection_output_path, f"thresh{ii}-{fov_name}-c{ch}-r0-z0.csv")
                bd_feats_sub.to_csv(bd_feats_fn)
                val_ = {"detected_images": os.path.basename(out_fn), "detected_outputs": os.path.basename(bd_feats_fn)}
                detected_outputs[f"channel_{ch}"][f"threshold_{ii}"] = val_

        # save detection information    
        output_data_fn = os.path.join(output_path_base, f"{fov_name}-detect-outputs.json")
        # input data
        with open(os.path.join(input_path, "primary", f"primary-{fov_name}.json"), 'r') as f:
            data_ = json.load(f)
            input_data = {"input_images": [x['file'] for x in data_['tiles']]}
        
        output_data0 = {'otsu_classes': otsu_n_classes, 
            'blob_thresholds': blob_thresholds.tolist(), 
            'histograms': saved_fns_dict['histograms']}
        
        output_data1 = {k:None for k in ch_keys}
        for ch in range(img.num_chs): 
            dat_1 = {k:input_data[k][ch] for k in input_data.keys()}
            dat_2 = {k:saved_fns_dict[k][ch] for k in list(saved_fns_dict.keys())[:-1]}
            output_data1[ch_keys[ch]] = {**dat_1, **dat_2, **detected_outputs[ch_keys[ch]]}
        
        with open(output_data_fn, 'w') as f:
            json.dump({**output_data0, **output_data1}, f, indent=4)


    # update data structure
    data_struct['num_thresholds'] = len(blob_thresholds)
    with open(data_struct_fn, "w") as f: 
        json.dump(data_struct, f, indent=4)
        


def main(): 
    parser = argparse.ArgumentParser(description="Detect mRNA spots.")
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
        min_sigma=step_config.get("min_sigma", 3),
        max_sigma=step_config.get("max_sigma", 10),
        num_sigma=step_config.get("num_sigma", 1),
        search_radius=step_config.get("search_radius", 10),
        otsu_n_classes=step_config.get("otsu_n_classes", 3)
    )
    
if __name__ == "__main__":
    main()
