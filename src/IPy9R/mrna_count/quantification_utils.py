# utils
import os

import numpy as np
from starfish import ImageStack
from starfish.types import Axes


def preprocess_image(image: ImageStack, 
    save_filtered=False, save_norm=False, output_path=None, fov_id=None): 
    """preprocess an image"""
    from starfish.core.types import Levels
    from starfish.image import Filter

    # bandpass filter --------------------------------------------
    # pre-filter clip to remove low-intensity background signal
    clip1 = Filter.Clip(p_min=50, p_max=100)
    bandpass = Filter.GaussianHighPass(sigma=9, is_volume=False, 
        level_method=Levels.SCALE_BY_IMAGE) # remove autofluorescence
    # post-filter clip to eliminate all but the highest-intensity peaks
    clip2 = Filter.Clip(p_min=95, p_max=100, is_volume=False)

    filtered = clip1.run(image, in_place=False, verbose=True, n_processes=8)
    bandpass.run(filtered, in_place=True, verbose=True, n_processes=8)
    clip2.run(filtered, in_place=True, verbose=True, n_processes=8)
    
    if save_filtered and (output_path is not None) and (fov_id is not None): 
        saved_filtered_fns = save_stack(filtered, "filtered", output_path, fov_id)
    else:
        saved_filtered_fns = None
        
    # normalize intensity distribution across channels
    mh_c = Filter.MatchHistograms({Axes.CH})
    scaled_c = mh_c.run(filtered, in_place=False, verbose=True, n_processes=8)
    
    if save_norm and (output_path is not None) and (fov_id is not None): 
        saved_norm_fns = save_stack(scaled_c, "normalized", output_path, fov_id)
    else:  
        saved_norm_fns = None
        
    # plot histograms
    if (output_path is not None): 
        saved_hist_fn = compare_histogram(image, filtered, scaled_c, 
            output_fn=os.path.join(output_path, f"{fov_id}-histograms.png"), 
            stack1_str="original", 
            stack2_str="filtered", 
            stack3_str="normalized")
    else: 
        saved_hist_fn = None
        
    # return normalized images and saved fns
    saved_fns_dict = {"filtered_images": saved_filtered_fns, 
        "normalized_images": saved_norm_fns, 
        "histograms": saved_hist_fn}
        
    return (scaled_c, saved_fns_dict)
    
    
def save_stack(image: ImageStack, save_type, output_path, fov_id): 
    """save images in an ImageStack"""    
    import tifffile
    
    num_rounds = image.shape['r']
    num_channels = image.shape['c']
    num_z_stacks = image.shape['z']
    
    saved_fn = []
    
    for r in range(num_rounds): 
        for c in range(num_channels): 
            for z in range(num_z_stacks): 
                # export images 
                out_fn = os.path.join(output_path, 
                    f"{save_type}-{fov_id}-c{c}-r{r}-z{z}.tiff")
                tifffile.imwrite(out_fn, 
                    image.sel({Axes.CH: c, Axes.ROUND: r, Axes.ZPLANE: z}).xarray.squeeze())
                saved_fn.append(os.path.basename(out_fn))
                
    return saved_fn
    
    
def compare_histogram(stack1: ImageStack, stack2: ImageStack, stack3: ImageStack, 
    output_fn, stack1_str="original", stack2_str="filtered", stack3_str="normalized"): 
    """compare histograms for each slice in the image stack"""
    import matplotlib
    from matplotlib import pyplot as plt
    from starfish.util.plot import intensity_histogram
    
    num_channels = stack1.shape['c']
        
    # visualize intensity histograms
    matplotlib.rcParams["figure.dpi"] = 450
    matplotlib.rcParams['axes.titlesize'] = 10
    matplotlib.rcParams['xtick.labelsize'] = 8
    matplotlib.rcParams['ytick.labelsize'] = 8
    f, axs = plt.subplots(3, num_channels, sharey=True)
    if num_channels == 1: 
        axs = np.expand_dims(axs, axis=1)
    f.suptitle('Intensity Histogram (log)')
    f.set_figheight(2 * num_channels)
    f.set_figwidth(4)
    f.tight_layout()

    # Plot intensity distribution of entire as a histogram with 50 bins
    for c in range(num_channels): 
        intensity_histogram(stack1, sel={Axes.CH: c}, log=True, bins=50, ax=axs[0][c],
            title=f"{stack1_str}:ch_{c}")
            
    for c in range(num_channels): 
        intensity_histogram(stack2, sel={Axes.CH: c}, log=True, bins=50, ax=axs[1][c],
            title=f"{stack2_str}:ch_{c}")

    for c in range(num_channels): 
        intensity_histogram(stack3, sel={Axes.CH: c}, log=True, bins=50, ax=axs[2][c],
            title=f"{stack3_str}:ch_{c}")

    # save output
    plt.savefig(output_fn)    
    
    return os.path.basename(output_fn)


def otsu_threshold(image: ImageStack, n_classes: int, 
    save_thresholded=False, output_path=None, fov_id=None): 
    """determine threshold for blob detection"""
    import tifffile
    from skimage.filters import threshold_multiotsu
    
    # for each channel
    num_channels = image.shape['c']
    
    # blob_thresholds = []
    # for c in range(num_channels): 
        # img_ = image.sel({Axes.CH: c}).xarray.squeeze().to_numpy()
        # thresh_ = threshold_multiotsu(img_, classes=n_classes)        
        # blob_thresholds.append(thresh_)    
        # regions = np.digitize(img_, bins=thresh_)
        # if save_thresholded and (output_path is not None) and (fov_id is not None): 
            # out_fn = os.path.join(output_path, 
                # f"otsu-{fov_id}-r0-c{c}-z0.tiff")
            # tifffile.imwrite(out_fn, regions.astype(np.uint8))
    
    # uniform thresholds for all channels
    img_list = []
    for c in range(num_channels): 
        img_list.append(image.sel({Axes.CH: c}).xarray.squeeze().to_numpy())
    img_concat = np.concatenate(img_list, axis=0)

    # thresholding 
    blob_thresholds = threshold_multiotsu(img_concat, classes=n_classes)   

    # save outputs
    for c in range(num_channels):
        if save_thresholded and (output_path is not None) and (fov_id is not None): 
            regions = np.digitize(img_list[c], bins=blob_thresholds)    
            out_fn = os.path.join(output_path, f"otsu-{fov_id}-c{c}-r0-z0.tiff")
            tifffile.imwrite(out_fn, regions.astype(np.uint8))
    
    # return thresholds
    return blob_thresholds
    
    
def load_mask(mask_fn): 
    """read a mask and process it""" 
    import cv2
    
    mask_in = cv2.imread(mask_fn, 0) 
    if mask_in is None:
        raise FileNotFoundError(f"Could not read mask file: {mask_fn}")
        
    # thresholding 
    _, mask_out = cv2.threshold(mask_in, 127, 255, cv2.THRESH_BINARY) # T = 0 before
    
    return mask_out    
    
    
def spots_count(bd_feats, mask_list, feats_out_fn, image_out_fn, contour_thickness=3): 
    """count puncta within a mask"""
    import tifffile
    from skimage import measure
    from skimage.draw import circle_perimeter, disk
    
    feats_ = bd_feats.copy()
    def __in_mask(x, y):
        x = int(round(x))
        y = int(round(y))
    
        for mask_id, mask_ in mask_list.items():
            h, w = mask_.shape[:2]
            if x < 0 or x >= w or y < 0 or y >= h:
                continue            
            if mask_[y, x]:
                return mask_id
        return "nan" 
                
    ret_ = feats_.apply(lambda row: __in_mask(row['x'], row['y']), axis=1)
    
    if feats_.shape[0] > 0: 
        feats_['in_mask'] = ret_ != "nan"
        feats_['mask_id'] = ret_
    else: 
        import pandas as pd
        feats_ = pd.concat([feats_, pd.DataFrame(columns=['in_mask', 'mask_id'])])
        
    # draw spots within a mask
    mask_out = np.zeros_like(mask_list[list(mask_list.keys())[0]], dtype=np.uint8)
    mask_shape = mask_out.shape
    # plot spots
    feats_sub = feats_.loc[feats_['in_mask']]
    for (xx, yy, rr) in zip(feats_sub['x'], feats_sub['y'], feats_sub['radius']):
        row, col = circle_perimeter(yy, xx, int(rr), shape=mask_shape)
        mask_out[row, col] = 1
    # plot mask boundary
    # visualize barrel contours
    #for i, (mask_id, mask_) in enumerate(mask_list.items()):
    for mask_ in mask_list.values():
        contours = measure.find_contours(mask_ > 0)
        for ctr in contours: 
            rc = [disk(c, contour_thickness, shape=mask_shape) for c in ctr]
            for rc_ in rc: 
                mask_out[rc_[0], rc_[1]] = 1
        
    # save outputs
    feats_sub = feats_sub.reset_index(drop=True)
    feats_sub.to_csv(feats_out_fn)
    
    tifffile.imwrite(image_out_fn, mask_out)
