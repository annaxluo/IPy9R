# src/IPy9R/mrna_count/summarize_spots_barrels.py

# Apply barrel masks and summarize spot statistics in each barrel 
import argparse
import json
import os

import numpy as np
import pandas as pd

from IPy9R.mrna_count.quantification_utils import load_mask, spots_count


def run(data_path_base,
        roi_list, 
        roi_area_thresh=0.4, 
        mask_suffix="_transformed.tif",
        roi_folder_prefix="confocal_63x_"
        ) -> None: 
    """
    Apply barrel masks and summarize spot statistics in each barrel
    """

    # 1. process paths ---------------------------------------------------------
    output_path_base = os.path.join(data_path_base, "SF_analysis", "SF_outputs")
    detection_path = os.path.join(data_path_base, "SF_analysis", "SF_outputs")

    output_path = os.path.join(output_path_base, "mask_outputs")
    os.makedirs(output_path, exist_ok=True)        

    # 2. read roi info and spot detection results----------------------------------------- 
    roi_summary_fn = os.path.join(data_path_base, "roi_summary.csv")
    roi_summary_df = pd.read_csv(roi_summary_fn, index_col=0)

    input_path = os.path.join(data_path_base, "SF_analysis") 
    data_struct_fn = os.path.join(input_path, "data_structure.json")
    with open(data_struct_fn, "r") as f: 
        data_struct = json.load(f)

    # 3. Process each ROI ----------------------------------------------------
    unique_roi_list = np.unique(roi_list)

    for used_roi in unique_roi_list: 
        print("processing ROI:", used_roi)
        # masks used
        roi_path_base = os.path.join(data_path_base, roi_folder_prefix + used_roi)
        mask_path = os.path.join(roi_path_base, "output_masks") 

        used_mask_ids = roi_summary_df.loc[
            (roi_summary_df['best_roi']==used_roi) & 
            (roi_summary_df['ratio'] > roi_area_thresh)]['mask_id'].tolist()
            
        if len(used_mask_ids) > 0: 
            # detection results 
            used_fovs = [f_n for f_n in data_struct['fov_info'].keys() 
                if (data_struct['fov_info'][f_n]['used_roi'] == used_roi)]
            
            for fov_name in used_fovs: 
                fov_info = data_struct['fov_info'][fov_name]
                with open(os.path.join(detection_path, f"{fov_info['fov_name']}-detect-outputs.json"), 'r') as f: 
                    detection_info = json.load(f)
                    
                # process each channel and threshold
                ch_keys = [ch for ch in range(len(fov_info['channels']))]
                th_keys = [k for k in range(len(detection_info['blob_thresholds']))]
                mask_detect_info = {f"channel_{ch}":None for ch in ch_keys}

                for ch in ch_keys: 
                    mask_detect_info[f"channel_{ch}"] = {f"threshold_{k}":None for k in th_keys}
                    for th in th_keys: 
                        bd_feats_fn = os.path.join(detection_path, "detection_outputs", 
                            f"thresh{th}-{fov_info['fov_name']}-c{ch}-r0-z0.csv")
                        bd_feats = pd.read_csv(bd_feats_fn, index_col=0)
                        mask_list = {k:None for k in used_mask_ids}
                        for mask_id_ in used_mask_ids: 
                            mask_ = load_mask(os.path.join(mask_path, mask_id_ + mask_suffix))
                            mask_list[mask_id_] = mask_
                        feats_out_fn =  os.path.join(output_path, f"thresh{th}-{fov_info['fov_name']}-c{ch}-r0-z0.csv")
                        image_out_fn = os.path.join(output_path, f"thresh{th}-{fov_info['fov_name']}-c{ch}-r0-z0.tif")
                        spots_count(bd_feats, mask_list, feats_out_fn, image_out_fn)
                        # update detection info 
                        mask_detect_info[f"channel_{ch}"][f"threshold_{th}"] = {
                            "mask_detected_images": os.path.basename(image_out_fn), 
                            "mask_detected_outputs": os.path.basename(feats_out_fn)}
                        print(f"processed {fov_name}, channel {ch}, threshold {th}")
                        
                # update detection info
                for ch in ch_keys: 
                    for th in th_keys: 
                        detection_info[f"channel_{ch}"][f"threshold_{th}"].update(
                            mask_detect_info[f"channel_{ch}"][f"threshold_{th}"])
                            
                with open(os.path.join(detection_path, f"{fov_info['fov_name']}-detect-outputs.json"), 'w') as f: 
                    json.dump(detection_info, f, indent=4)
                    
def main(): 
    parser = argparse.ArgumentParser(description="Summarize spots per barrel.")
    parser.add_argument("--config", required=True)

    args = parser.parse_args()

    from IPy9R.config import get_step_config, load_config

    config = load_config(args.config)
    step_config = get_step_config(config, "mrna_count")

    run(
        data_path_base=step_config["data_path_base"],
        roi_list=step_config["roi_list"], 
        roi_area_thresh=step_config.get("roi_area_thresh", 0.4), 
        mask_suffix=step_config.get("mask_suffix", "_transformed.tif"), 
        roi_folder_prefix=step_config.get("roi_folder_prefix", "confocal_63x_")
    )
    
if __name__ == "__main__":
    main()