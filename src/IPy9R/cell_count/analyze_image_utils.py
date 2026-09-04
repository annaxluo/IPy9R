# analyze images with masks and cellprofiler outputs 
import os

import cv2
import numpy as np
import pandas as pd
import tifffile


def cell_count(img_fn, mask_fn_list, data_fn, contour_thickness=3, circle_radius=20): 
    """count the number of cells within a mask region""" 
    img = cv2.imread(img_fn, 0)
    
    sel_data = pd.read_csv(data_fn)
    sel_in_mask = pd.DataFrame(columns=["X", "Y", "mask_id"])
    img_out = np.dstack([img, img, img])
    for mask_fn in mask_fn_list: 
        mask_id = os.path.split(mask_fn)[1].split("_")[0]
        mask = cv2.imread(mask_fn, 0)

        # threshold mask 
        _, mask = cv2.threshold(mask, 127, 255, cv2.THRESH_BINARY)
        # find selection points within the mask 
        in_mask_ = pd.DataFrame(columns=["X", "Y", "mask_id"])
        for x, y in zip(sel_data["X"], sel_data["Y"]): 
            if mask[int(y), int(x)]: 
                in_mask_.loc[len(in_mask_.index)] = [x, y, mask_id]
                # visualize 
                contours, hierarchy = cv2.findContours(mask, 
                    cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
                cv2.drawContours(img_out, contours, -1, (0,255,0), 
                    contour_thickness)
                # draw selections 
                for x1, y1 in zip(in_mask_["X"], in_mask_["Y"]): 
                    cv2.circle(img_out, (int(x1), int(y1)), circle_radius, (255, 0, 0), -1)
                    
        sel_in_mask = pd.concat([sel_in_mask, in_mask_])
        
        print("processed mask", mask_id)
        
    return img_out, sel_in_mask
    
    
def DAPI_count(nuclei_img_fn, mask_fn_list, nuclei_data_fn, 
    contour_thickness=1, circle_radius=1, 
    x_str="Location_CenterMassIntensity_X_Nuclei", 
    y_str="Location_CenterMassIntensity_Y_Nuclei"): 
    """count the number of DAPI nuclei in each mask"""
    img = cv2.imread(nuclei_img_fn, 0)
    
    nuclei_dat = pd.read_csv(nuclei_data_fn) 
    nuclei_in_mask = pd.DataFrame(columns=["X", "Y", "mask_id"])
    img_out = np.dstack([img, img, img])
    
    for mask_fn in mask_fn_list: 
        mask_id = os.path.split(mask_fn)[1].split("_")[0]
        mask = cv2.imread(mask_fn, 0)
        
        # threshold mask 
        _, mask = cv2.threshold(mask, 127, 255, cv2.THRESH_BINARY)
        # find selection points within the mask 
        in_mask_ = pd.DataFrame(columns=["X", "Y", "mask_id"])
        
        for x, y in zip(nuclei_dat[x_str], 
                        nuclei_dat[y_str]): 
            x = round(x)
            y = round(y)
            if mask[int(y), int(x)]: 
                in_mask_.loc[len(in_mask_.index)] = [x, y, mask_id]
                # visualize 
                contours, hierarchy = cv2.findContours(mask, 
                    cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
                cv2.drawContours(img_out, contours, -1, (0,255,0), 
                    contour_thickness)
                # draw selections 
                for x1, y1 in zip(in_mask_["X"], in_mask_["Y"]): 
                    cv2.circle(img_out, (int(x1), int(y1)), circle_radius, (255, 0, 0), -1)
        
        nuclei_in_mask = pd.concat([nuclei_in_mask, in_mask_])
        print("processed mask", mask_id)
    
    return img_out, nuclei_in_mask    
    
    
def select_in_mask(in_img_fn, mask_fn_list, point_data_fn, 
    contour_thickness=1, circle_radius=1, 
    x_str="Location_CenterMassIntensity_X_Nuclei", 
    y_str="Location_CenterMassIntensity_Y_Nuclei"): 
    """select points within each mask"""
    
    img = cv2.imread(in_img_fn, 0)
    
    point_dat = pd.read_csv(point_data_fn) 
    point_in_mask = pd.DataFrame(columns=[x_str, y_str, "mask_id"])
    img_out = np.dstack([img, img, img])
    
    for mask_fn in mask_fn_list: 
        mask_id = os.path.split(mask_fn)[1].split("_")[0]
        mask = cv2.imread(mask_fn, 0)
        
        # threshold mask 
        _, mask = cv2.threshold(mask, 127, 255, cv2.THRESH_BINARY)
        # find selection points within the mask 
        in_mask_ = pd.DataFrame(columns=[x_str, y_str, "mask_id"])
        
        for x, y in zip(point_dat[x_str], 
                        point_dat[y_str]): 
            x = round(x)
            y = round(y)
            if mask[int(y), int(x)]: 
                in_mask_.loc[len(in_mask_.index)] = [x, y, mask_id]
                # visualize 
                contours, hierarchy = cv2.findContours(mask, 
                    cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)
                cv2.drawContours(img_out, contours, -1, (0,255,0), 
                    contour_thickness)
                # draw selections 
                for x1, y1 in zip(in_mask_[x_str], in_mask_[y_str]): 
                    cv2.circle(img_out, (int(x1), int(y1)), circle_radius, (255, 0, 0), -1)
        
        point_in_mask = pd.concat([point_in_mask, in_mask_])
        print("processed mask", mask_id)
    
    return img_out, point_in_mask    
    
def summarize_count(sel_in_mask, groups, mask_fn_list): 
    """create summary of cell count stats"""
    # list of mask ids
    mask_ids = []
    for mask_fn in mask_fn_list: 
        mask_ids.append(os.path.split(mask_fn)[1].split("_")[0])
    
    summary = pd.DataFrame(columns=["mask_id", "count", "group"])
    
    for mask_id in mask_ids: 
        count = sum(sel_in_mask["mask_id"] == mask_id)
        m_group = [grp for grp, m_id in groups.items() if mask_id in m_id]
        if not m_group:
            raise ValueError(f"Mask {mask_id} is not assigned to any group")
        summary.loc[len(summary.index)] = [mask_id, count, m_group[0]]
        
    return summary


def cell_density(summary, img_fn, mask_fn_list): 
    """compute density of cells within the masked region""" 
    # get resolution information from tif file 
    tif_tags = {}
    with tifffile.TiffFile(img_fn) as tif: 
        for tag in tif.pages[0].tags.values(): 
            name, val = tag.name, tag.value
            tif_tags[name] = val
    if len(tif_tags) == 0: 
        print("cannot open tiff file.")
        return None
        
    x_res = tif_tags["XResolution"]
    y_res = tif_tags["YResolution"]
    resolution_factor = (x_res[1]/x_res[0])*(y_res[1]/y_res[0]) # micron2 per pixel2
    
    summary_out = summary.copy()
    summary_out["Area"] = ""
    summary_out["Density"] = ""
    for mask_fn in mask_fn_list: 
        mask_id = os.path.split(mask_fn)[1].split("_")[0]
        row_id = summary_out.loc[summary_out['mask_id'] == mask_id].index[0]
        
        mask = cv2.imread(mask_fn, 0)

        # threshold mask 
        _, mask = cv2.threshold(mask, 127, 255, cv2.THRESH_BINARY)
        m_area = cv2.countNonZero(mask) * resolution_factor # in micron^2
        m_density = summary_out.loc[row_id, ["count"]][0] / m_area
        summary_out.loc[row_id, ["Area", "Density"]] = [m_area, m_density]
    
    return summary_out
    

def cell_to_nuclei_ratio(summary_cell, summary_nuclei): 
    summary_out = summary_cell.copy()
    assert (summary_out["mask_id"] == summary_nuclei["mask_id"]).all()
    summary_out["nuclei_count"] = summary_nuclei["count"]
    
    summary_out["ratio"] = summary_out["count"] / summary_nuclei["count"]
    
    return summary_out
    
