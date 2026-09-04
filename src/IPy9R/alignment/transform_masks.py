# src/IPy9R/alignment/transform_masks.py

# Apply optimal scaling factors on all barrel masks to transform them  from lower 
# magnification to higher magnification. 

import argparse
import os
import sys

import cv2
import tifffile

from IPy9R.alignment.align_images_utils import TemplateMatchingScale


def run(param_fn, 
        mask_dir, 
        used_masks,   
        boundary_mask_fn, 
        output_dir
        ) -> None: 
    """
    Apply optimal scaling factors on all barrel masks to transform them  from lower 
    magnification to higher magnification. 
    """
    os.makedirs(output_dir, exist_ok=True)
    # load matcher 
    matcher = TemplateMatchingScale(param_fn)
    
    mask_list = [fn for fn in os.listdir(mask_dir) if (".tif" in fn) and 
        fn.split("_")[0] in used_masks]
    
    # detect image boundary 
    boundary_mask = cv2.imread(boundary_mask_fn, 0) # need to be binary 
    if boundary_mask is None:
        sys.exit(f"Could not read boundary mask: {boundary_mask_fn}")
    
    # tranform masks 
    for mask_fn in mask_list: 
        mask_in = cv2.imread(os.path.join(mask_dir, mask_fn), 0) 
        if mask_in is None:
            sys.exit(f"Could not read mask: {os.path.join(mask_dir, mask_fn)}")
            
        try: 
            mask_out = matcher.transform_mask(mask_in)
        except cv2.error as e: 
            sys.exit(str(e))
        
        if mask_out is None: 
            sys.exit("Mask transformation failed.")
        
        if boundary_mask.shape != mask_out.shape:
            sys.exit(
                f"Boundary mask shape {boundary_mask.shape} does not match "
                f"transformed mask shape {mask_out.shape}."
            )
            
        # add image boundary 
        mask_out2 = cv2.bitwise_and(mask_out, mask_out, mask=boundary_mask)
        #mask_out2 = mask_out
        out_fn = os.path.join(output_dir, mask_fn.replace(".tif", "_transformed.tif"))
        tifffile.imwrite(out_fn, mask_out2)
        print("transformed mask", mask_fn)
        
        
def main(): 
    parser = argparse.ArgumentParser(description="Transform all masks.")
    parser.add_argument("--config", required=True)

    args = parser.parse_args()

    from IPy9R.config import get_step_config, load_config

    config = load_config(args.config)
    step_config = get_step_config(config, "transform_masks")

    run(
        param_fn=step_config["param_fn"],
        mask_dir=step_config["mask_dir"], 
        used_masks=step_config["used_masks"],   
        boundary_mask_fn=step_config["boundary_mask_fn"], 
        output_dir=step_config["output_dir"], 
    )
    
if __name__ == "__main__":
    main()