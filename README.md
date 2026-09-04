# IPy9R

[![CI](https://github.com/annaxluo/IPy9R/actions/workflows/ci-light.yml/badge.svg)](https://github.com/annaxluo/IPy9R/actions/workflows/ci-light.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

`IPy9R` is a Python package that contains the image-processing workflows in my thesis project to quantify smFISH data. The package currently supports three main analyses:

1. **Image alignment**
- Estimate a model to project lower-magnification images (e.g., 10x) to higher-magnification images (e.g., 63x) using multi-scale template matching followed by perspective transformation. 
- Transform annotations of barrels in the mouse posteromedial barrel subfield (PMBSF) made in lower-magnification images to barrel mask for higher-magnification images.

2. **Cell counting**
- Preprocess images for `CellProfiler` by creating tiles.
- Combine `CellProfiler` tile outputs, and aggregate counts of DAPI-positive nuclei and specific cell types (e.g., microglia) inside each barrel of the PMBSF.
- Compute cell density and cell-to-nuclei ratios for each barrel. 

3. **mRNA spot counting**
- Make data structure for confocal image data compatible for `Starfish`.
- Detect mRNA puncta using `BlobDetector` in `Starfish`.
- Count spots inside individual barrels.
- Summarize spot counts across ROIs and experiments.

The package uses a `src/` layout and is configured through YAML files in the `configs/` directory.

---

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/annxluo/IPy9R.git
cd IPy9R
```

### 2. Create a virtual environment

Using Conda:

```bash
conda create -n ipy9r python=3.10
conda activate ipy9r
```

### 3. Install the package

For lightweight development installation:

```bash
pip install -e .
```

Some workflow-specific dependencies are optional because packages such as OpenCV, Starfish, and SlicedImage can be heavy.

To install alignment dependencies:

```bash
pip install -e ".[alignment]"
```

To install cell-counting dependencies:

```bash
pip install -e ".[cell-count]"
```

To install mRNA-counting dependencies:

```bash
pip install -e ".[mrna-count]"
```

To install all optional workflow dependencies:

```bash
pip install -e ".[all]"
```

For development tools:

```bash
pip install -e ".[dev]"
```

---

## Configuration files

The package can be run from YAML configuration files.

The main example config files are:

```text
configs/alignment.yaml
configs/cell_count.yaml
configs/mrna_count.yaml
```

Modify inputs in these configuration files. 

---

# Task 1: Image alignment

To evaluate gene expression in the mouse posteromedial barrel subfield (PMBSF), we identified a method to method to delineate barrel boundaries using autofluorescence signal in the GFP channel in 10x images acquired on a widefield microscope. 

However, green autoflourescence was weaker at the higher resolutions (e.g., 20x, 63x) at which the smFISH signals were detected. To enable identification of barrel boundaries at higher resolutions, I designed an image process pipeline that includes the following steps: 
- Takes in fluorescence images of the same channel (e.g., DAPI) acquired at low and high resolution, and estimate the best scaling factor for the low-resolution image via **multi-scale template matching**. 
- Scales the low-resolution image, and estimate the **perspective transformation** matrix that project the low-resolution image to high-resolution image. 
- Applies the optimal scaling factor and perspective transformation matrix (defined as an "Alignment Parameter Set") to barrel annotations (binary masks of individual barrels) made a low resolution to estimate barrel masks at high resolution. 

<p align="center">
  <img src="docs/assets/multi-scale-alignment.jpg" alt="Projecting barrel masks from lower-magnification to higher-magnification" width="700">
</p>

<p align="center">
  rojecting barrel masks from lower-magnification to higher-magnification. 
</p>


The relevant config file is:

```text
configs/alignment.yaml
```

---

## 1.1 Update `alignment.yaml`

Example structure:

```yaml
alignment:
  param_fn: "/path/to/confocal_20x/alignment_params.json"

  image_fn: "/path/to/confocal_20x/CP_inputs/C3-slide-stack.tif"

  template_fn: "/path/to/keyence/input_images/keyence_image.tif"

  val_output_dir: "/path/to/confocal_20x/output_masks/val"

  init_mask_coords: [1551, 3363, 500, 400]

  scale_min: 0.8
  scale_max: 1.2

  image_crop_ul: null
  image_crop_hw: null
  image_scale_factor: 0.60

  template_crop_ul: null
  template_crop_hw: null
  template_scale_factor: 0.12

transform_masks:
  param_fn: "/path/to/confocal_20x/alignment_params.json"

  mask_dir: "/path/to/keyence/input_masks"

  used_masks:
    - delta
    - A1
    - A2
    - B1
    - B2
    - B3
    - C1
    - C2
    - C3
    - D1
    - D2
    - D3
    - E1
    - E2

  boundary_mask_fn: "/path/to/confocal_20x/output_masks/image_boundary_mask.tif"

  output_dir: "/path/to/confocal_20x/output_masks"
```

Important fields:

| Field | Meaning |
|---|---|
| `alignment.param_fn` | Output JSON file where alignment parameters are saved. |
| `alignment.image_fn` | High-magnification target image. |
| `alignment.template_fn` | Lower-magnification reference/template image. |
| `alignment.val_output_dir` | Directory for validation images. |
| `alignment.init_mask_coords` | Initial template mask coordinates as `[Y, X, H, W]`. Use `null` for default. |
| `alignment.scale_min`, `alignment.scale_max` | Search range for scale matching. |
| `transform_masks.mask_dir` | Directory containing original barrel masks. |
| `transform_masks.used_masks` | Mask IDs to transform. |
| `transform_masks.boundary_mask_fn` | Binary mask defining the valid target image boundary. |
| `transform_masks.output_dir` | Directory where transformed masks are written. |

---

## 1.2 Run template matching

```bash
python -m IPy9R.alignment.run_template_matching --config configs/my_alignment.yaml
```

This step writes alignment parameters to:

```text
alignment.param_fn
```

and validation images to:

```text
alignment.val_output_dir
```

Inspect the validation outputs before proceeding.

---

## 1.3 Transform barrel masks

```bash
python -m IPy9R.alignment.transform_masks --config configs/my_alignment.yaml
```

This step reads the masks listed under `transform_masks.used_masks`, transforms them using the saved alignment parameters, applies the image boundary mask, and writes outputs such as:

```text
B2_transformed.tif
C3_transformed.tif
delta_transformed.tif
```

to:

```text
transform_masks.output_dir
```

---

# Task 2: Cell counting

The cell-counting workflow has three main steps:

1. Preprocess images into tiles for CellProfiler.
2. Analyze CellProfiler nuclei outputs.
3. Compute cell statistics inside barrel masks.

The relevant config file is:

```text
configs/cell_count.yaml
```

---

## 2.1 Update `cell_count.yaml`

Example structure:

```yaml
shared:
  data_path: /data_base
  stack_str: stack01
  slide_id: slide01
  nuclei_ch: C3
  rna_ch: C1

  mask_stack:
    - B2
    - B3
    - C2
    - C3
    - C4
    - C5
    - C6
    - D1
    - D2
    - D3
    - D4
    - D5
    - D6
    - D7
    - delta
    - E1
    - E2
    - E3

preprocess_cellprofiler:
  max_h: 2000
  max_w: 2000

analyze_cellcount: {}

compute_stats:
  spared_deprived_groups_fn: "configs/spared_deprived_barrels.json"
```

Important fields:

| Field | Meaning |
|---|---|
| `shared.data_path` | Base path for one cell-counting dataset. |
| `shared.stack_str` | Stack identifier, for example `stack01`. |
| `shared.slide_id` | Slide/image identifier. |
| `shared.nuclei_ch` | DAPI/nuclei channel ID. |
| `shared.rna_ch` | RNA/cell-marker channel ID. |
| `shared.mask_stack` | Barrel mask IDs to include. |
| `preprocess_cellprofiler.max_h` | Maximum tile height for CellProfiler inputs. |
| `preprocess_cellprofiler.max_w` | Maximum tile width for CellProfiler inputs. |
| `compute_stats.spared_deprived_groups_fn` | JSON file assigning masks to spared/deprived groups. |

Expected input structure under `data_path`:

```text
<data_path>/
├── CP_inputs/
│   ├── C3-<slide_id>-<stack_str>.tif
│   └── C1-<slide_id>-<stack_str>.tif
├── CP_outputs/
└── output_masks/
    ├── B2_transformed.tif
    ├── B3_transformed.tif
    └── ...
```

---

## 2.2 Preprocess images for CellProfiler

```bash
python -m IPy9R.cell_count.preprocess_cellprofiler --config configs/my_cell_count.yaml
```

This creates image tiles under:

```text
<data_path>/CP_inputs/image_tiles_<stack_str>/
```

Use these tiles as CellProfiler inputs.

---

## 2.3 Run CellProfiler externally

After preprocessing, run the appropriate CellProfiler pipeline outside this package.

The expected CellProfiler output file is:

```text
<data_path>/CP_outputs/<stack_str>_Nuclei_obj.csv
```

The tile metadata file should be:

```text
<data_path>/CP_inputs/image_tiles_<stack_str>/<nuclei_ch>-<slide_id>-<stack_str>_tiles.csv
```

---

## 2.4 Analyze CellProfiler nuclei outputs

```bash
python -m IPy9R.cell_count.analyze_cellcount --config configs/my_cell_count.yaml
```

This step:

1. Combines tile-level CellProfiler outputs.
2. Selects nuclei inside transformed masks.
3. Writes visualization and CSV outputs.

Expected outputs include:

```text
<data_path>/CP_outputs/<stack_str>_Nuclei_obj_combinedTiles.csv
<data_path>/CP_outputs/<stack_str>_Nuclei_obj_inMask.tif
<data_path>/CP_outputs/<stack_str>_Nuclei_obj_inMask.csv
```

---

## 2.5 Manually curate nuclei selections if needed

If manual curation is needed, edit/export the curated nuclei file as:

```text
<data_path>/CP_outputs/<stack_str>_Nuclei_obj_inMask_v2.csv
```

This file is used by the statistics step.

---

## 2.6 Compute cell statistics

```bash
python -m IPy9R.cell_count.compute_stats --config configs/my_cell_count.yaml
```

This step expects:

```text
<data_path>/CP_outputs/<stack_str>_Nuclei_obj_inMask_v2.csv
<data_path>/CP_outputs/CellCount_<rna_ch>-<slide_id>-<stack_str>.csv
```

Outputs include:

```text
<data_path>/CP_outputs/<stack_str>_cell_inMask_summary.csv
<data_path>/CP_outputs/<stack_str>_cell_ratio_inMask_summary.csv
```

---

# Task 3: mRNA spot counting

The mRNA-counting workflow has four main steps:

1. Structure data for Starfish.
2. Detect mRNA spots.
3. Summarize spots inside barrel masks for each ROI.
4. Summarize spots across the experiment.

The relevant config file is:

```text
configs/mrna_count.yaml
```

---

## 3.1 Update `mrna_count.yaml`

Example structure:

```yaml
mrna_count:
  data_path_base: "/path/to/experiment"

  roi_folder_prefix: "confocal_63x_"

  roi_list:
    - ROI1
    - ROI2
    - ROI3
    - ROI4
    - ROI5

  used_channels:
    - C2
    - C4

  channel_gene_names:
    - geneA
    - geneB

  image_id_list:
    - "image1"
    - "image2"
    - "image3"
    - "image4"
    - "image5"

  used_z_list:
    - 3
    - 3
    - 3
    - 3
    - 3

  selection_rect:
    - null
    - null
    - null
    - null
    - null

  z_thickness: 1.0

  min_sigma: 3
  max_sigma: 10
  num_sigma: 1
  search_radius: 10
  otsu_n_classes: 3

  roi_area_thresh: 0.4

  mask_suffix: "_transformed.tif"

  spared_deprived_groups_fn: "configs/spared_deprived_barrels.json"
```

Important fields:

| Field | Meaning |
|---|---|
| `data_path_base` | Base experiment folder. |
| `roi_folder_prefix` | Prefix used for ROI folders, for example `confocal_63x_`. |
| `roi_list` | ROI identifiers to process. |
| `used_channels` | Image channels used for mRNA counting. |
| `channel_gene_names` | Gene names corresponding to `used_channels`. |
| `image_id_list` | Image identifiers for each ROI. |
| `used_z_list` | Z-plane index for each ROI. |
| `selection_rect` | Optional crop rectangle for each ROI. Use `null` for no crop. |
| `z_thickness` | Z thickness used in Starfish coordinate metadata. |
| `min_sigma`, `max_sigma`, `num_sigma` | Blob detector scale parameters. |
| `search_radius` | Spot detector search radius. |
| `otsu_n_classes` | Number of Otsu classes used to derive thresholds. |
| `roi_area_thresh` | Minimum ROI overlap ratio for including a barrel. |
| `mask_suffix` | Suffix for transformed mask files. |
| `spared_deprived_groups_fn` | JSON file defining spared/deprived barrel groups. |

Expected input structure:

```text
<data_path_base>/
├── roi_summary.csv
├── confocal_63x_ROI1/
│   ├── CP_inputs/
│   │   ├── C2-image1-0003.tif
│   │   └── C4-image1-0003.tif
│   └── output_masks/
│       ├── B2_transformed.tif
│       └── ...
├── confocal_63x_ROI2/
│   └── ...
└── SF_analysis/
```

---

## 3.2 Structure data for Starfish

```bash
python -m IPy9R.mrna_count.structure_data --config configs/my_mrna_count.yaml
```

This step creates structured Starfish input data under:

```text
<data_path_base>/SF_analysis/
```

It writes files such as:

```text
<data_path_base>/SF_analysis/data_structure.json
<data_path_base>/SF_analysis/primary/
<data_path_base>/SF_analysis/primary_dir/
```

---

## 3.3 Detect mRNA spots

```bash
python -m IPy9R.mrna_count.detect_mrna_spots --config configs/my_mrna_count.yaml
```

This step loads the Starfish experiment, preprocesses images, computes thresholds, detects spots, and writes detection outputs.

Expected output folders:

```text
<data_path_base>/SF_analysis/SF_outputs/preprocess_outputs/
<data_path_base>/SF_analysis/SF_outputs/threshold_outputs/
<data_path_base>/SF_analysis/SF_outputs/detection_outputs/
```

Example detection outputs:

```text
thresh0-f0-c0-r0-z0.csv
thresh0-f0-c0-r0-z0.tiff
thresh1-f0-c1-r0-z0.csv
```

The step also updates:

```text
<data_path_base>/SF_analysis/data_structure.json
```

with the number of thresholds.

---

## 3.4 Summarize spots inside barrel masks

```bash
python -m IPy9R.mrna_count.summarize_spots_barrels --config configs/my_mrna_count.yaml
```

This step applies transformed barrel masks to detected spots and writes mask-level detection outputs.

Expected output folder:

```text
<data_path_base>/SF_analysis/SF_outputs/mask_outputs/
```

Example outputs:

```text
thresh0-f0-c0-r0-z0.csv
thresh0-f0-c0-r0-z0.tif
```

---

## 3.5 Summarize spots across the experiment

```bash
python -m IPy9R.mrna_count.summarize_spots_experiment --config configs/my_mrna_count.yaml
```

This step combines ROI-level and barrel-level spot counts across the experiment.

Expected output:

```text
<data_path_base>/SF_analysis/SF_outputs/experiment_counts.csv
```

---

# Spared/deprived barrel groups

Both the cell-counting and mRNA-counting workflows can use a JSON file defining spared and deprived barrel groups.

Example:

```json
{
  "spared": ["B2", "B3", "C2", "C3"],
  "deprived": ["D1", "D2", "D3"]
}
```

The default path used in the example configs is:

```text
configs/spared_deprived_barrels.json
```

Update this file to match your experiment.

---

# Development and checks

Install development tools:

```bash
pip install -e ".[dev]"
```

Run Ruff:

```bash
ruff check src
```

Run syntax/compile check:

```bash
python -m compileall -q src
```

Validate package metadata:

```bash
validate-pyproject pyproject.toml
```

Build source distribution and wheel:

```bash
python -m build --sdist --wheel
```

Run all lightweight local checks:

```bash
validate-pyproject pyproject.toml
python -m compileall -q src
ruff check src
python -m build --sdist --wheel
```

Generated folders such as `build/`, `dist/`, `*.egg-info/`, and `__pycache__/` should not be committed to GitHub.

---

# Continuous integration

This repository uses GitHub Actions for lightweight CI. The workflow checks package metadata, source compilation, Ruff linting, and package building.

The CI badge at the top of this README shows the current status of the workflow.

---

# License

This project is licensed under the GNU General Public License v3.0.

See the GPL v3 license terms here:

https://www.gnu.org/licenses/gpl-3.0.en.html
```

Also make sure your repository has an actual license file:

```text
LICENSE
```

For GPL v3, you can add it from GitHub:

```text
Add file > Create new file > LICENSE
```

Then choose:

```text
GNU General Public License v3.0
```