# src/IPy9R/cell_count/preprocess_cellprofiler.py

# preprocess images for CellProfiler nuclei segmentation by creating image tiles.
import argparse
import os

from IPy9R.cell_count.cellprofiler_utils import make_tiles


def run(data_path, stack_str, slide_id, nuclei_ch, rna_ch, max_h=2000, max_w=2000):
    """Preprocess images for CellProfiler nuclei segmentation by creating image tiles."""    
    # img1: nuclei channel 
    img1_fn = os.path.join(data_path, "CP_inputs", f"{nuclei_ch}-{slide_id}-{stack_str}.tif")
    # img2: mRNA channel 
    img2_fn = os.path.join(data_path, "CP_inputs", f"{rna_ch}-{slide_id}-{stack_str}.tif")

    tile_output_path = os.path.join(data_path, "CP_inputs", f"image_tiles_{stack_str}")
    os.makedirs(tile_output_path, exist_ok=True)
    
    n_tiles_nuclei = make_tiles(img1_fn, tile_output_path, max_h, max_w)
    n_tiles_rna = make_tiles(img2_fn, tile_output_path, max_h, max_w)
    
    return {
        "n_tiles_nuclei": n_tiles_nuclei,
        "n_tiles_rna": n_tiles_rna,
        "tile_output_path": tile_output_path, 
        }


def main(): 
    parser = argparse.ArgumentParser(description="Create image tiles for CellProfiler.")
    
    # add arguments from config
    parser.add_argument("--config", required=True)

    args = parser.parse_args()
    
    from IPy9R.config import get_step_config, load_config

    config = load_config(args.config)
    step_config = get_step_config(config, "preprocess_cellprofiler")

    result = run(
        data_path=step_config["data_path"],
        stack_str=step_config["stack_str"],
        slide_id=step_config["slide_id"],
        nuclei_ch=step_config["nuclei_ch"],
        rna_ch=step_config["rna_ch"],
        max_h=step_config.get("max_h", 2000),
        max_w=step_config.get("max_w", 2000),
    )

    print(result)
    
    
if __name__ == "__main__":
    main()