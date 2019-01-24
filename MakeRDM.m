function [rdm] = MakeRDM(roi,dist)
%% --------------------- Script Description -----------------------------
% Script for creating a dissimilarity matrix using the specified distance
% metric taking in the output from ExtractROI.m.
%
% Takes in the following inputs:
% 1) roi - data structure from ExtractROI.m with the following req fields:
%   i) descrip - 1*M cell string array with M being the number of conds
%  ii) roi_val - 1*M cell array each containing a mat for the values in
%                each ROI in the respective conditions
% iii) info.name - name of the ROI that values was extracted from
%
% 2) dist - distace metric to use. Takes the following string -
%           'correlation','spearman','euclidean' etc (refer to pdist2.m)
%
% Generate a structure with the representational dissimilarity matrix and
% the condition labels
%
% Completed by JH 7/12/2017
%
% Version history
% Edits 15/3/2018
% Include the name of the ROI in the output for easier referencing
%%
rdm = struct();
rdm.dist = dist;

cond_name = roi.output.descrip;
temp = roi.output.roi_val;

numvx = roi.info.size;
numcond = length(temp);

all_cond_mat = zeros(numcond,numvx);

for i = 1 : numcond
    all_cond_mat(i,:) = temp{i};
end

mat = pdist(all_cond_mat,dist);
rdm.mat = squareform(mat);
rdm.cond = cond_name;
rdm.roi = roi.info.name;

end
