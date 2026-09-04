# src/IPy9R/cell_count/compute_stats.py

# After mannually curating Nuclei selection in ImageJ, reprocess the nuclei segmentation, 
# compute microglia density and microglia-to-nuclei ratio. 
import argparse
import json
import os

from IPy9R.cell_count.analyze_image_utils import (
    cell_density,
    cell_to_nuclei_ratio,
    select_in_mask,
    summarize_count,
)


def run(data_path, stack_str, slide_id, nuclei_ch, rna_ch, mask_stack, spared_deprived_groups_fn) -> None: 
    """Process nuclei segmentation outputs from CellProfiler."""
    # process inputs 
    nuclei_img_fn = os.path.join(data_path, "CP_inputs", 
        f"{nuclei_ch}-{slide_id}-{stack_str}.tif")
        
    with open(spared_deprived_groups_fn, "r") as f: 
        spared_deprived_groups = json.load(f)
        
    mask_path = os.path.join(data_path, "output_masks")
    mask_fn_list = [m for m in os.listdir(mask_path) if ("_transformed.tif" in m) and 
        (m.split("_")[0] in mask_stack)]
    mask_fn_list = [os.path.join(mask_path, m) for m in mask_fn_list]
    
    rna_img_fn = os.path.join(data_path, "CP_inputs", 
        f"{rna_ch}-{slide_id}-{stack_str}.tif")
    
    curated_nuclei_fn = os.path.join(data_path, "CP_outputs", 
        f"{stack_str}_Nuclei_obj_inMask_v2.csv")
    if not os.path.exists(curated_nuclei_fn):
        raise FileNotFoundError(f"Curated nuclei file not found: {curated_nuclei_fn}")
        
    nuclei_img_out2, nuclei_in_mask2 = select_in_mask(
        in_img_fn=nuclei_img_fn, 
        mask_fn_list=mask_fn_list, 
        point_data_fn=curated_nuclei_fn, 
        contour_thickness=1, 
        circle_radius=1, 
        x_str="X", 
        y_str="Y"
    )
    
    # microglia in mask 
    curated_cell_fn = os.path.join(data_path, "CP_outputs", 
        f"CellCount_{rna_ch}-{slide_id}-{stack_str}.csv")
        
    cell_img_out, cell_in_mask = select_in_mask(
        in_img_fn=rna_img_fn, 
        mask_fn_list=mask_fn_list, 
        point_data_fn=curated_cell_fn, 
        contour_thickness=1, 
        circle_radius=1, 
        x_str="X", 
        y_str="Y"
    )

    # summary stats 
    summary_cell = summarize_count(cell_in_mask, spared_deprived_groups, mask_fn_list)

    summary_nuclei = summarize_count(nuclei_in_mask2, spared_deprived_groups, mask_fn_list)

    # compute cell density 
    cell_density_df = cell_density(summary_cell, rna_img_fn, mask_fn_list)

    out_fn3 = os.path.join(data_path, "CP_outputs", 
        f"{stack_str}_cell_inMask_summary.csv")
    cell_density_df.to_csv(out_fn3, sep=",", header=True, index=True)

    # compute cell to nuclei ratio 
    cell_ratio = cell_to_nuclei_ratio(cell_density_df, summary_nuclei)
    out_fn4 = os.path.join(data_path, "CP_outputs", 
        f"{stack_str}_cell_ratio_inMask_summary.csv")
    cell_ratio.to_csv(out_fn4, sep=",", header=True, index=True)

def main(): 
    parser = argparse.ArgumentParser(description="Compute microglia statistics.")
    parser.add_argument("--config", required=True)

    args = parser.parse_args()

    from IPy9R.config import get_step_config, load_config

    config = load_config(args.config)
    step_config = get_step_config(config, "compute_stats")

    run(
        data_path=step_config["data_path"],
        stack_str=step_config["stack_str"],
        slide_id=step_config["slide_id"],
        nuclei_ch=step_config["nuclei_ch"],
        rna_ch=step_config["rna_ch"],
        mask_stack=step_config["mask_stack"],
        spared_deprived_groups_fn=step_config["spared_deprived_groups_fn"],
    )
    
if __name__ == "__main__":
    main()
