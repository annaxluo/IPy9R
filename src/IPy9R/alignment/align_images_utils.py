# align images 
import json
import os

import cv2
import numpy as np
import tifffile


class TemplateMatchingScale: 
    """match template to image at multiple scales""" 
    
    def __init__(self, parameter_fn): 
    
        # file names 
        self.image_fn = None 
        self.template_fn = None 
        self.initial_mask_coords = None # [YX, HW] of initial mask, in pixel
    
        self.image_shape = None # shape of original image 
        self.image_crop_ul = None # upperleft corner of crop
        self.image_crop_hw = None # height and width of crop 
        self.image_scale_factor = None       
        self.image = None # needed?? 
        self.image_resized = None # needed?? 
        
        self.template_shape = None # shape of original template 
        self.template_crop_ul = None
        self.template_crop_hw = None
        self.template_scale_factor = None
        self.template = None # needed?? 
        self.template_resized = None # needed?? 
        
        self.initial_mask = None # needed?? 
        self.initial_mask_resized = None # needed?? 
        
        # best scale factor for template matching  
        self.tm_image_resized_shape = None # shape for mask transform 
        self.tm_image_ul = None # found: [x, y]
        self.tm_image_scale_factor = None # found: scaling factor 
        self.tm_perspective_M = None # matrix for perspective transform 
        self.tm_image_resized = None # needed?? 
        
        # file to store tuned parameters 
        self.parameters_fn = parameter_fn 
        
        # try to read parameters 
        self.__read_params()
        
    
    def __read_params(self): 
        """read parameters for template matching"""
        if not os.path.exists(self.parameters_fn): 
            return
            

        if os.path.exists(self.parameters_fn): 
            with open(self.parameters_fn, 'r') as f:
                # Reading from json file
                json_object = json.load(f)
                
                self.image_fn = json_object["image_fn"]
                self.template_fn = json_object["template_fn"]
                self.initial_mask_coords = json_object["initial_mask_coords"]
                
                self.image_shape = json_object["image_shape"]
                self.image_crop_ul = json_object["image_crop_ul"]
                self.image_crop_hw = json_object["image_crop_hw"]
                self.image_scale_factor = json_object["image_scale_factor"]
                
                self.template_shape = json_object["template_shape"]
                self.template_crop_ul = json_object["template_crop_ul"]
                self.template_crop_hw = json_object["template_crop_hw"]
                self.template_scale_factor = json_object["template_scale_factor"]
                
                self.tm_image_resized_shape = json_object["tm_image_resized_shape"]
                self.tm_image_ul = json_object["tm_image_ul"]
                self.tm_image_scale_factor = json_object["tm_image_scale_factor"]
                self.tm_perspective_M = json_object["tm_perspective_M"]


    def __write_param(self): 
        """write parameters to files""" 
        out = {
            "image_fn": self.image_fn, 
            "template_fn": self.template_fn, 
            "initial_mask_coords": self.initial_mask_coords, 
            "image_shape": self.image_shape, 
            "image_crop_ul": self.image_crop_ul,
            "image_crop_hw": self.image_crop_hw,
            "image_scale_factor": self.image_scale_factor,
            "template_shape": self.template_shape, 
            "template_crop_ul": self.template_crop_ul,             
            "template_crop_hw": self.template_crop_hw, 
            "template_scale_factor": self.template_scale_factor, 
            "tm_image_resized_shape": self.tm_image_resized_shape, 
            "tm_image_ul": self.tm_image_ul, 
            "tm_image_scale_factor": self.tm_image_scale_factor, 
            "tm_perspective_M": self.tm_perspective_M
        }
   
        # Serializing json
        json_object = json.dumps(out, indent=4)
 
        # Writing to sample.json
        with open(self.parameters_fn, "w") as f:
            f.write(json_object)
    
    
    def __adjust_crop_size(self, crop_wh, scale_factor): 
        """auto-adjust crop size according to scale_factor to ensure integer size  
        """
        while (scale_factor * crop_wh[0]) - int((scale_factor * crop_wh[0])) > 0: 
            crop_wh[0] = int(int((scale_factor * crop_wh[0])) / scale_factor)
        
        while (scale_factor * crop_wh[1]) - int((scale_factor * crop_wh[1])) > 0: 
            crop_wh[1] = int(int((scale_factor * crop_wh[1])) / scale_factor)
        
        return crop_wh
        
        
    def read_image(self, image_fn, scale_factor, crop_ul=None, crop_hw=None): 
        """read image"""                                                        
        img = cv2.imread(image_fn, 0) # read in grayscale
        
        if img is None:
            raise FileNotFoundError(f"Could not read image: {image_fn}")

        # crop a rectangular region 
        if crop_ul is None or crop_hw is None: 
            crop_ul = [0, 0]
            crop_hw = list(img.shape) 
            
        crop_hw = self.__adjust_crop_size(crop_hw, scale_factor)
        
        img_cropped = img[crop_ul[0]:crop_ul[0]+crop_hw[0], 
            crop_ul[1]:crop_ul[1]+crop_hw[1]]
      
        # resize image
        img_resized = cv2.resize(img_cropped, 
            None, 
            fx = scale_factor, 
            fy = scale_factor, 
            interpolation = cv2.INTER_AREA)
            
        # for instance using to compute transformation parameters  
        self.image = img 
        self.image_resized = img_resized 
        
        self.image_fn = image_fn
        self.image_shape = img.shape 
        self.image_crop_ul = crop_ul 
        self.image_crop_hw = crop_hw
        self.image_scale_factor = scale_factor    
        
        self.__write_param()

        
    def read_template(self, template_fn, scale_factor, crop_ul=None, crop_hw=None, 
        mask_coords=None, default_mask_size=.1): 
        """read template""" 
        
        templ = cv2.imread(template_fn, 0) # read in grayscale
        if templ is None:
            raise FileNotFoundError(f"Could not read template: {template_fn}")

        # crop a rectangular region 
        if crop_ul is None or crop_hw is None: 
            crop_ul = [0, 0]
            crop_hw = list(templ.shape) 
            
        crop_hw = self.__adjust_crop_size(crop_hw, scale_factor)
        
        templ_cropped = templ[crop_ul[0]:crop_ul[0]+crop_hw[0], 
            crop_ul[1]:crop_ul[1]+crop_hw[1]]
      
        # resize image
        templ_resized = cv2.resize(templ_cropped, 
            None, 
            fx = scale_factor, 
            fy = scale_factor, 
            interpolation = cv2.INTER_AREA)
            
        # for instance using to compute transformation parameters  
        self.template = templ
        self.template_resized = templ_resized   
        self.template_fn = template_fn         
        
        self.template_shape = templ.shape
        self.template_crop_ul = crop_ul 
        self.template_crop_hw = crop_hw   
        self.template_scale_factor = scale_factor
        
        # read initial mask 
        mask = mask_resized = None
        if mask_coords is not None: 
            mask, mask_resized = self.__make_initial_mask(mask_coords)
        else: # use an initial mask centering the template 
            mask, mask_resized = self.__make_initial_mask_default(default_mask_size)
                
        self.initial_mask = mask 
        self.initial_mask_resized = mask_resized
        
        self.initial_mask_coords = mask_coords
        
        self.__write_param()
        
        
    def process_mask(self, mask): 
        """read a mask and process it according to template matching parameters""" 
        
        # thresholding 
        _, mask = cv2.threshold(mask, 127, 255, cv2.THRESH_BINARY) # T = 0 before
        
        crop_ul = self.template_crop_ul
        crop_hw = self.template_crop_hw
        scale_factor = self.template_scale_factor
        
        mask_cropped = mask[crop_ul[0]:crop_ul[0]+crop_hw[0], 
            crop_ul[1]:crop_ul[1]+crop_hw[1]]
        mask_resized = cv2.resize(mask_cropped, 
            None, 
            fx = scale_factor, 
            fy = scale_factor, 
            interpolation = cv2.INTER_NEAREST)
        
        return mask_resized    
        
    
    def __make_initial_mask(self, mask_coords): 
        """makes an initial mask"""
        mask = np.zeros(self.template_shape, dtype="uint8")
        mask_hw = (mask_coords[2], mask_coords[3])
        mask_ul = (mask_coords[0], mask_coords[1])
        mask[mask_ul[0]-mask_hw[0]:mask_ul[0]+mask_hw[0], 
             mask_ul[1]-mask_hw[1]:mask_ul[1]+mask_hw[1]] = 255
        # crop and resize 
        mask_resized = self.process_mask(mask)
        
        return mask, mask_resized 
        
        
    def __make_initial_mask_default(self, default_mask_size=.1): 
        """makes the default initial mask"""
        mask = np.zeros(self.template_shape, dtype="uint8")
        mask_hw = (int(mask.shape[0]*default_mask_size/2), 
                   int(mask.shape[1]*default_mask_size/2))
        mask_ul = (int(mask.shape[0] / 2), int(mask.shape[1] /2))
        mask[mask_ul[0]-mask_hw[0]:mask_ul[0]+mask_hw[0], 
             mask_ul[1]-mask_hw[1]:mask_ul[1]+mask_hw[1]] = 255
        # crop and resize 
        mask_resized = self.process_mask(mask)

        return mask, mask_resized 
        

    def template_matching_scale(self, scales): 
        """template matching at multiple scales"""
        # check required parameters 
        flag = False 
        if self.template_resized is None: 
            print("read template first.")
            flag = True 
        
        if self.image_resized is None: 
            print("read image first.")
            flag = True 
            
        if flag: 
            return None, None
        
        template_edged = cv2.Canny(self.template_resized, 50, 200)
        (tH, tW) = self.template_resized.shape[:2]
        
        found = None 
        found_img = None
        for scale in scales: 
            # match template to image 
            image_resized = cv2.resize(self.image_resized, 
                None, 
                fx=scale, 
                fy=scale, 
                interpolation = cv2.INTER_AREA)
                
            edged = cv2.Canny(image_resized, 50, 200)
            result = cv2.matchTemplate(edged, 
                template_edged, 
                cv2.TM_CCOEFF, 
                mask = self.initial_mask_resized)
            (_, maxVal, _, maxLoc) = cv2.minMaxLoc(result)

            if found is None or maxVal > found[0]:
                found = (maxVal, maxLoc, scale)
                # for the instance to compute transformation parameters 
                self.tm_image_resized = image_resized
                clone = np.dstack([image_resized, image_resized, image_resized])
                found_img = cv2.rectangle(clone, (maxLoc[0], maxLoc[1]), 
                    (maxLoc[0] + tW, maxLoc[1] + tH), (0, 0, 255), 2)
                
            if image_resized.shape[0] < tH or image_resized.shape[1] < tW: 
                break        
        
        # add parameters 
        self.tm_image_resized_shape = self.tm_image_resized.shape 
        self.tm_image_ul = found[1]
        self.tm_image_scale_factor = found[2]
        
        # TODO: write template matching parameters to json 
        self.__write_param()
        
        return found, found_img 
        
 
    def template_matching_scale_perspective(self): 
        """Template matching using homography"""
        # check required parameters 
        if self.template_resized is None: 
            print("read template first and perfrom template matching first.")
            return False       
        
        if (self.tm_image_scale_factor is None) or (self.tm_image_resized is None): 
            print('Perform template matching first.')
            return False
        
        templ = self.template_resized      
        img = self.tm_image_resized
        
        # Initiate SIFT detector
        MIN_MATCH_COUNT = 10
        sift = cv2.SIFT_create()
        
        # find the keypoints and descriptors with SIFT
        kp1, des1 = sift.detectAndCompute(templ, None)
        kp2, des2 = sift.detectAndCompute(img, None)
        if des1 is None or des2 is None:
            print("Could not compute SIFT descriptors.")
            return False
        
        FLANN_INDEX_KDTREE = 1
        index_params = dict(algorithm = FLANN_INDEX_KDTREE, trees = 5)
        search_params = dict(checks = 50)
        flann = cv2.FlannBasedMatcher(index_params, search_params)
        matches = flann.knnMatch(des1, des2, k=2)
        
        # store all the good matches as per Lowe's ratio test.
        good = []
        for pair in matches:
            if len(pair) < 2:
                continue
            m, n = pair
            if m.distance < 0.7*n.distance:
                good.append(m)
        
        if len(good) > MIN_MATCH_COUNT:
            src_pts = np.float32([ kp1[m.queryIdx].pt for m in good ]).reshape(-1,1,2)
            dst_pts = np.float32([ kp2[m.trainIdx].pt for m in good ]).reshape(-1,1,2)
            M, mask = cv2.findHomography(src_pts, dst_pts, cv2.RANSAC, 5.0)
            if M is None:
                print("Homography estimation failed.")
                return False
            # for the instance to compute transformation parameters 
            self.tm_perspective_M = M.tolist()
        else:
            print("Not enough matches are found - {}/{}".format(len(good), MIN_MATCH_COUNT))
            return False           
        
        # TODO: write template matching parameters to json 
        self.__write_param()
        
        return True
        
        
    def transform_mask(self, mask_in=None): 
        """perform perspective transform and rescaling of a mask on the 'template' to 
           produce a mask for the image. 
        """  
        # check required parameters 
        if self.image_scale_factor is None: 
            print("read image first")
            return None
            
        if (self.tm_image_scale_factor is None) or (self.tm_image_resized_shape is None): 
            print("perform template matching first.")
            return None

        if self.tm_perspective_M is None: 
            print("perform perspective transform first.")
            return None
        
        # perspective transform of input mask
        (img_H, img_W) = self.tm_image_resized_shape
        
        if mask_in is None: 
            if self.initial_mask is None: 
                mask_in, _ = self.__make_initial_mask_default()
            else: 
                mask_in = self.initial_mask
        
        # resize mask_in
        mask_resized = self.process_mask(mask_in)
        
        mask_out = cv2.warpPerspective(mask_resized, 
            np.array(self.tm_perspective_M), 
            (img_W, img_H))
            
        # rescale back to fit original image 
        r = 1 / (self.image_scale_factor * self.tm_image_scale_factor)
        mask_out1 = cv2.resize(mask_out, 
            None, 
            fx=r, 
            fy=r, 
            interpolation = cv2.INTER_NEAREST)
            
        # TODO: resize to target shape, ie self.image_crop_hw??
        # pad image 
        top = self.image_crop_ul[0]
        left = self.image_crop_ul[1]
        bottom = self.image_shape[0] - top - mask_out1.shape[0]
        right = self.image_shape[1] - left - mask_out1.shape[1]
        
        mask_out2 = cv2.copyMakeBorder(
            mask_out1, 
            top, bottom, left, right, 
            cv2.BORDER_CONSTANT, 
            None, 
            value = 0)
        
        return mask_out2
        
    
    def map_inverse(self, image_in): 
        """inversely maps image to template"""
        # check required parameters 
        if self.image_scale_factor is None: 
            print("read image first")
            return None
            
        if (self.tm_image_scale_factor is None) or (self.tm_image_resized_shape is None): 
            print("perform template matching first.")
            return None

        if self.tm_perspective_M is None: 
            print("perform perspective transform first.")
            return None
        
        # perspective transform of input mask
        (img_H, img_W) = self.tm_image_resized_shape
        
        
    def vis(self, out_fn_str): 
        """write matched or transformed mask for validation"""
        # transformed template and mask 
        templ = self.template_resized 
        templ_out = self.transform_mask(templ)
        init_mask_out = self.transform_mask()
        
        out_fn1 = out_fn_str + '_tranformed_template.tif'
        tifffile.imwrite(out_fn1, templ_out)
        out_fn2 = out_fn_str + '_tranformed_template_mask.tif'
        tifffile.imwrite(out_fn2, init_mask_out)
        
        # mask can be used to overlay with confocal images

