function[roi] = ExtractROI_SPM(roinii,roistr,filelist,binary)
%% ------------------- Script Description ------------------------------
% Script for extracting values within a 3D map (beta or t map). Uses
% functions from SPM. 
% *Check output values and orientation especially if files are not
% generated using SPM.
%
% Inputs:
% 1) roinii: A string pointing to the nii file of the ROI matrix
% 2) roistr: A string containing the name of the ROI. This will also be the
%            name that the ROI will be saved as.
% 3) filelist: An array with pointers to each of the nii files where the
%              values would be extracted from. Typically a t or beta
%              map.
% 4) binary: default = 1; for non-binary image input 0.
% Outputs a single structure with the following fields:
% 1) name : name of the image 
% 2) mat : matrix value of the source image
% 3) spatialsnr : rough estimate of the spatial snr of the image based on
%                 the mean and standard deviation of values (for reference
%                 only, depends on the nature of your  data)
% 4) roi_val : values of the source image in the specified ROI
% 5) roi_zscore : z-scored values of roi-val (z-score within ROI)
% 6) roi_mean : mean value within the ROI
%
% MNI coordinates for each ROI will also be generated. MNI conversion taken
% from Xu Cui's function(cor2mni.m), and has not been checked.
%
% Completed by JH - 5/12/2017
%
% Notes:
% Function has been tested using spm_vol & spm_read_vols from SPM8. SPM12
% has not changed the 2 functions, so it should work similarly. 
% This has NOT been tested with other versions of SPM.
% 
% Version history:
% Edits 14/3/2018
% Modify output structure for easier viewing. 
% Added field to save the descrip of the extracted file
%
% Edits 22/3/2018
% Modified binary checking to ensure that the ROI is binary.
% Modified script to take in non-binary image with input checking
%
% Edits 27/9/2018
% Commented out the line saving each mat for each extracted volume to save
% space. Can uncomment if require those mats (Line 120).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
if nargin == 3
    binary = 1;
end
roi = struct();

% Create ROI information
roi2load = roinii;
r = spm_vol(roi2load);
rmat = spm_read_vols(r);
r_size = size(rmat);

% Check and ensure ROI is binary
if binary == 1
    min_val = min(min(min(rmat)));
    max_val = max(max(max(rmat)));
    
    chkbin = find(rmat~=0 & rmat~=1);
    
    if min_val ~= 0 || max_val ~= 1 || ~isempty(chkbin)
        error('ROI not binarised')
    end
end

r_pos = find(rmat ~= 0);
num_vx = size(r_pos,1);

roi.info.name = roistr;
roi.info.mat = rmat;
roi.info.size = num_vx;

% Get coordinates of the masked voxels
[x y z] = ind2sub(r_size,r_pos);
roi.info.coord = [x y z];

% convert voxel coordinate to MNI -- Taken from Xu Cui, need to check again
T = r.mat;
vx_mni = T * [x(:,1) y(:,1) z(:,1) ones(num_vx,1)]';
vx_mni = vx_mni';
vx_mni(:,4) = [];

roi.info.mni = vx_mni;

for i = 1 : length(filelist)
    disp(sprintf('Extracting roi value from file %2d',i));
    img2load = filelist{i};
    img = spm_vol(img2load);
    imgmat = spm_read_vols(img);
    img_size = size(imgmat);
    
    % Get rough spatial SNR for each of the image map
    snr = mean(imgmat(~isnan(imgmat))) / std(imgmat(~isnan(imgmat)));
    
    % Check img dimension with roi dimensions
    imgT = img.mat;
    
    if sum(sum(imgT ~= T))~=0 || sum(r_size~= img_size)~=0
        error('Image dimension mismatch')
    end
    
    if binary ~= 1
        w_imgmat = imgmat .* rmat;
        imgmat = w_imgmat;
    end
    
    img_val = imgmat(r_pos);
    img_z =zscore(img_val); 
    img_mean = mean(img_val);
    
    roi.output.name{i} = img2load;
    roi.output.descrip{i} = img.descrip;
    %roi.output.mat{i} = imgmat;
    roi.output.spatialsnr{i} = snr;
    roi.output.roi_val{i} = img_val;
    roi.output.roi_zscore{i} = img_z;
    roi.output.roi_mean{i} = img_mean;
    
%     roi.output(i,1).name = img2load;
%     roi.output(i,1).mat = imgmat;
%     roi.output(i,1).spatialsnr = snr;
%     roi.output(i,1).roi_val = img_val;
%     roi.output(i,1).roi_zscore = img_z;
%     roi.output(i,1).roi_mean = img_mean;
end
disp(sprintf('Completed for ROI: %s',roistr));
end
        
    
    
