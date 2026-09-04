# src/IPy9R/mrna_count/summarize_spots_experiment.py

# Combine results across ROIs and summarize spot statistics for an experiment. 
import argparse
import json
import os

import numpy as np
import pandas as pd


def run(data_path_base,
        roi_list, 
        spared_deprived_groups_fn,
        roi_area_thresh=0.4
        ) -> None: 
    """
    Combine results across ROIs and summarize spot statistics for an experiment.
    """

    # 1. process paths ---------------------------------------------------------
    output_path_base = os.path.join(data_path_base, "SF_analysis", "SF_outputs")
    mask_detection_path = os.path.join(data_path_base, "SF_analysis", "SF_outputs", "mask_outputs")

    with open(spared_deprived_groups_fn, "r") as f: 
        spared_deprived_groups = json.load(f)

    # 2. read roi info and spot detection results----------------------------------------- 
    roi_summary_fn = os.path.join(data_path_base, "roi_summary.csv")
    roi_summary_df = pd.read_csv(roi_summary_fn, index_col=0)

    input_path = os.path.join(data_path_base, "SF_analysis") 
    data_struct_fn = os.path.join(input_path, "data_structure.json")
    with open(data_struct_fn, "r") as f: 
        data_struct = json.load(f)
    # check required fields
    if "num_thresholds" not in data_struct:
        raise KeyError(
            "data_structure.json does not contain 'num_thresholds'. "
            "Run detect_mrna_spots.py before summarize_spots_experiment.py."
        )
    
    # 3. combine ROIs ------------------------------------------------------
    roi_out_df = roi_summary_df.copy()
    num_channels = data_struct['num_channels']
    num_thresholds = data_struct['num_thresholds']

    ll = [f"counts_t{i}_ch{j}" for i in range(num_thresholds) for j in range(num_channels)]
    roi_out_df2 = pd.DataFrame(columns=['mask_id'] + ll)

    unique_roi_list = np.unique(roi_list)

    for used_roi in unique_roi_list: 
        # roi info
        used_mask_ids = roi_summary_df.loc[
            (roi_summary_df['best_roi']==used_roi) & 
            (roi_summary_df['ratio'] > roi_area_thresh)]['mask_id'].tolist()
        ret = pd.DataFrame(columns=['mask_id'] + ll)
        
        # fov info
        used_fovs = [f_n for f_n in data_struct['fov_info'].keys() 
            if (data_struct['fov_info'][f_n]['used_roi'] == used_roi)]
            
        for mask_id in used_mask_ids: 
            row_ = [mask_id]
            for thr in range(num_thresholds): 
                for ch in range(num_channels): 
                    mask_detection_ = []
                    for fov_name in used_fovs:
                        mask_detected_fn = os.path.join(mask_detection_path, 
                            f"thresh{thr}-{fov_name}-c{ch}-r0-z0.csv")
                        if not os.path.exists(mask_detected_fn):
                            raise FileNotFoundError(
                                f"Missing mask detection file: {mask_detected_fn}. "
                                "Run summarize_spots_barrels.py before summarize_spots_experiment.py."
                            )
                            
                        mask_detection_.append(pd.read_csv(mask_detected_fn, index_col=0))
                    # count puncta
                    mask_detection = pd.concat(mask_detection_, axis=0, ignore_index=True)
                    num_spots = sum(mask_detection['mask_id']==mask_id)
                    row_.append(num_spots)
            ret.loc[len(ret.index)] = row_
        roi_out_df2 = pd.concat([roi_out_df2, ret])
        

    roi_out_df = roi_out_df.join(roi_out_df2.set_index('mask_id'), on='mask_id')
    for thr in range(num_thresholds): 
        for ch in range(num_channels): 
            roi_out_df[f"ratio_t{thr}_ch{ch}"] = \
                roi_out_df[f"counts_t{thr}_ch{ch}"] / roi_out_df['63x_area']

    spared_deprived = ["deprived" if s in spared_deprived_groups['deprived'] else "spared" 
        for s in roi_out_df['mask_id']]
    roi_out_df['spared_deprived'] = spared_deprived

    # print outputs
    for thr in range(num_thresholds): 
        for ch in range(num_channels): 
            c_str = f"ratio_t{thr}_ch{ch}"
            d_m = np.nanmean(roi_out_df.loc[roi_out_df['spared_deprived']=="deprived"][c_str])
            s_m = np.nanmean(roi_out_df.loc[roi_out_df['spared_deprived']=="spared"][c_str])
            print(f"threshold {thr}, channel {ch}")
            print(f"mean density for deprived barrels: {d_m: .5f}, spared barrels: {s_m: .5f}")

    out_fn = os.path.join(output_path_base, "experiment_counts.csv")
    roi_out_df.to_csv(out_fn)
    
def main(): 
    parser = argparse.ArgumentParser(description="Summarize spots per experiment.")
    parser.add_argument("--config", required=True)

    args = parser.parse_args()

    from IPy9R.config import get_step_config, load_config

    config = load_config(args.config)
    step_config = get_step_config(config, "mrna_count")
        
    run(
        data_path_base=step_config["data_path_base"],
        roi_list=step_config["roi_list"], 
        spared_deprived_groups_fn=step_config["spared_deprived_groups_fn"],
        roi_area_thresh=step_config.get("roi_area_thresh", 0.4)
    )
    
if __name__ == "__main__":
    main()
    

