# src/IPy9R/cell_count/analyze_cellcount.py

# Process outputs from CellProfiler (nuclei segmentation)
import argparse
import os

import tifffile

from IPy9R.cell_count.analyze_image_utils import select_in_mask
from IPy9R.cell_count.cellprofiler_utils import combine_tiles


def run(data_path, stack_str, slide_id, nuclei_ch, mask_stack) -> None: 
    """Process nuclei segmentation outputs from CellProfiler."""
    data_fn = os.path.join(data_path, "CP_outputs", f"{stack_str}_Nuclei_obj.csv")
    tile_data_fn = os.path.join(data_path, "CP_inputs", 
        f"image_tiles_{stack_str}", f"{nuclei_ch}-{slide_id}-{stack_str}_tiles.csv") 
        
    # nuclei image
    nuclei_img_fn = os.path.join(data_path, "CP_inputs", 
        f"{nuclei_ch}-{slide_id}-{stack_str}.tif")
        
    # barrel masks 
    mask_path = os.path.join(data_path, "output_masks") 

    # 1. Combine nuclei objects from files--------------------------
    output_fn = os.path.join(data_path, "CP_outputs", 
        f"{stack_str}_Nuclei_obj_combinedTiles.csv")
        
    _ = combine_tiles(
        tile_data_fn=tile_data_fn, 
        selection_fn=data_fn, 
        output_fn=output_fn) # ret: True 
        
        
    # 2. Select nuclei objects within barrel masks------------------------------
    nuclei_fn = os.path.join(data_path, "CP_outputs", 
        f"{stack_str}_Nuclei_obj_combinedTiles.csv")
    
    # read masks
    mask_fn_list = [m for m in os.listdir(mask_path) if ("_transformed.tif" in m) and 
        (m.split("_")[0] in mask_stack)]
    mask_fn_list = [os.path.join(mask_path, m) for m in mask_fn_list]

    # count DAPI nuclei 
    nuclei_img_out, nuclei_in_mask = select_in_mask(
        in_img_fn=nuclei_img_fn, 
        mask_fn_list=mask_fn_list, 
        point_data_fn=nuclei_fn, 
        contour_thickness=1, 
        circle_radius=1, 
        x_str="Location_CenterMassIntensity_X_Nuclei", 
        y_str="Location_CenterMassIntensity_Y_Nuclei"
    )
        
    # save outputs
    out_fn1 = os.path.join(data_path, "CP_outputs", 
        f"{stack_str}_Nuclei_obj_inMask.tif")
    tifffile.imwrite(out_fn1, nuclei_img_out)

    out_fn2 = os.path.join(data_path, "CP_outputs", 
        f"{stack_str}_Nuclei_obj_inMask.csv")
    nuclei_in_mask.to_csv(out_fn2, sep=",", header=True, index=True)


def main(): 
    parser = argparse.ArgumentParser(description="Process CellProfiler outputs.")
    parser.add_argument("--config", required=True)

    args = parser.parse_args()

    from IPy9R.config import get_step_config, load_config

    config = load_config(args.config)
    step_config = get_step_config(config, "analyze_cellcount")

    run(
        data_path=step_config["data_path"],
        stack_str=step_config["stack_str"],
        slide_id=step_config["slide_id"],
        nuclei_ch=step_config["nuclei_ch"],
        mask_stack=step_config["mask_stack"],
    )
    
if __name__ == "__main__":
    main()
