# src/IPy9R/alignment/run_template_matching.py

# Perform multi-scale template matching on `image` and `template`, 
# and identify the optimal scaling factor. And then scale the `image` using the 
# optimal scaling, and compute the perspective transformation parameters between scaled 
# `image` and `template`. Write a validation image using the optimal 
# transformation parameters for visual inspection.

import argparse
import os
import sys
import time

import cv2
import numpy as np
import tifffile

from IPy9R.alignment.align_images_utils import TemplateMatchingScale


def run(param_fn, 
        image_fn, 
        template_fn,  
        val_output_dir, 
        init_mask_coords, 
        scale_min,
        scale_max,         
        image_crop_ul=None, 
        image_crop_hw=None, 
        image_scale_factor=0.60, 
        template_crop_ul=None, 
        template_crop_hw=None, 
        template_scale_factor=0.10): 
    """
    Perform multi-scale template matching on `image` and `template`, 
    and identify the optimal scaling factor. Write a validation image using the optimal 
    transformation parameters for visual inspection. 
    """
    os.makedirs(val_output_dir, exist_ok=True)

    # template matching -----------------------------------------------------
    matcher = TemplateMatchingScale(param_fn)
    matcher.read_image(
        image_fn, 
        image_scale_factor,
        crop_ul=image_crop_ul, 
        crop_hw=image_crop_hw
    )
    matcher.read_template(
        template_fn, 
        template_scale_factor, 
        crop_ul=template_crop_ul, 
        crop_hw=template_crop_hw, 
        mask_coords=init_mask_coords)
    
    # initial scales
    scales = np.linspace(scale_min, scale_max, 20)
    
    t_start = time.time()
    
    try: 
        found, found_img = matcher.template_matching_scale(scales)   
        if found is None:
            raise ValueError("Template matching failed: no valid match found.")    
            
        opt_idx = int(np.argmin(np.abs(scales - found[2])))
        
        if opt_idx==0: 
            raise ValueError(
                f"Borderline scaling factor: {found[2]}. "
                "Decrease the lower limit of `scales`."
            )
        
        if opt_idx==(len(scales)-1): 
            raise ValueError(
                f"Borderline scaling factor: {found[2]}. "
                "Increase the upper limit of `scales`."
            )
            
    except cv2.error as e: 
        sys.exit(str(e))
        
    except ValueError as e:
        sys.exit(str(e))
        
    print("template matching 1, elapsed time: ", time.time() - t_start)


    # refine template matching --------------------------------------
    opt_idx = int(np.argmin(np.abs(scales - found[2])))
    scales2 = np.linspace(scales[max(0, opt_idx-1)], scales[opt_idx+1], 100)
    try: 
        found, found_img = matcher.template_matching_scale(scales2)   
        if found is None:
            raise ValueError("Template matching failed: no valid match found.")  
            
        opt_idx2 = int(np.argmin(np.abs(scales2 - found[2])))
        
        if opt_idx2==0: 
            raise ValueError(
                f"Borderline scaling factor: {found[2]}. "
                "Decrease the lower limit of `scales`."
            )
        
        if opt_idx2==(len(scales2)-1): 
            raise ValueError(
                f"Borderline scaling factor: {found[2]}. "
                "Increase the upper limit of `scales`."
            )
            
    except cv2.error as e: 
        sys.exit(str(e))
        
    except ValueError as e:
        sys.exit(str(e))
    
    print("template matching 2, elapsed time: ", time.time() - t_start)
        
    # compute perspective transformation parameters --------------------------------
    try: 
        ret = matcher.template_matching_scale_perspective()
        if not ret:
            raise ValueError("Perspective transform failed.")
    except cv2.error as e: 
        sys.exit(str(e))
    except ValueError as e:
        sys.exit(str(e))
        
    print("perspective transform, elapsed time: ", time.time() - t_start)

    # write validation images 
    out_fn0 = os.path.join(val_output_dir, "val_init_mask_templ.tif")
    tifffile.imwrite(out_fn0, matcher.initial_mask)

    init_mask_out = matcher.transform_mask()
    if init_mask_out is None:
        sys.exit("Initial mask transformation failed.")
    out_fn1 = os.path.join(val_output_dir, "val_init_mask_out.tif")
    tifffile.imwrite(out_fn1, init_mask_out)


def main(): 
    parser = argparse.ArgumentParser(description="Perform multi-scale template matching.")
    parser.add_argument("--config", required=True)

    args = parser.parse_args()

    from IPy9R.config import get_step_config, load_config

    config = load_config(args.config)
    step_config = get_step_config(config, "alignment")

    run(
        param_fn=step_config["param_fn"], 
        image_fn=step_config["image_fn"], 
        template_fn=step_config["template_fn"],  
        val_output_dir=step_config["val_output_dir"], 
        init_mask_coords=step_config["init_mask_coords"], 
        scale_min=step_config["scale_min"],
        scale_max=step_config["scale_max"],
        image_crop_ul=step_config["image_crop_ul"], 
        image_crop_hw=step_config["image_crop_hw"], 
        image_scale_factor=step_config["image_scale_factor"], 
        template_crop_ul=step_config["template_crop_ul"], 
        template_crop_hw=step_config["template_crop_hw"], 
        template_scale_factor=step_config["template_scale_factor"], 
    )
    
if __name__ == "__main__":
    main()