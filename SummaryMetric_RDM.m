function [rdm_sm] = SummaryMetric_RDM(rdm,ds)
%% --------------------- Script Description -----------------------------
% Script for extracting summary metric from RDM taking in output from
% MakeRDM.m. Extracts the within and between category distances from RDMs.
% Includes both Item and Category level summary.
%
% Takes in the following inputs:
% 1) rdm - structure created from MakeRDM.m with the following req fields:
%   i) mat - the representational dissimilarity matrix
%  ii) cond - the name of each condition in the rdm
% iii) roi - name of roi rdm was created from
%
% 2) ds - structure with the following required field:
%   i) convect - M*N matrix containing the conditions to compare, with M
%                being the number of items and N being the number
%                of comparisons. Comparisons, in this context, refers to
%                definintion of within and between category comparisons.
%                Conditions can take on any integers, with conditions
%                sharing the same label being considered as 'within' and
%                conditions with different labels as 'between' (for X
%                conditions, each condition will have X - 1 'between'
%                category distance). Ignores '0'.
%  ii) conname - 1*N cell string array with name of the comparisons.
%
% Produces the output structure rdm_sm with the following fields saved
% under the respective comparisons:
% 1) catg - contains the category averages
%       i) within : X * 1 vector showing within category distance for all
%                   categories.
%     iii) between : X * 1 vector averaging all between category
%                        distances
%       v) info : X * 1 vector of Between_avg - Within category
%                     distances
%
% 2) trial - Contains trial level distances
%       i) within: M * 1 vector containing the trial's distance to other
%                  items in the same category.
%      ii) between: M * 1 vector averaging all between category
%                       distances
%     iii) info_avg : M * 1 matrix of Between_avg - Within
%
% 3) avg - contains average all all categories
%       i) within : single value for within category distances
%      ii) between: single value for between category distances
%     iii) info : single value for Between - Within
%
% Completed by JH  28/3/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
rdm_sm = struct();
rdm_sm.roi = rdm.roi;
numitem = size(rdm.mat,1);
[M, N] = size(ds.convect);
numCond = max(ds.convect);

if M ~= numitem
    error('Number of items in convect does not match RDM');
end

for i = 1 : N
    cond = ds.convect(:,i);
    rdm_sm.conname{i} = ds.conname{i};
    
    full_itemw = zeros(numitem,1);
    full_itemb = zeros(numitem,1);
    
    numCond = max(cond);
    
    for j = 1 : numCond
        % Identify conditions
        pos1 = find(cond == j);
        pos2 = find(cond > 0 & cond ~= j);
        
        % Get Within-Category sub-matrix
        with_mat = rdm.mat(pos1, pos1);
        item_w = sum(with_mat,2) / (size(with_mat,2)-1);
        full_itemw(pos1,1) = item_w;
        
        rdm_sm.con{i}.catg.within(j,1) = mean(item_w);
        
        % Between-Category sub-matrix
        btw_mat = rdm.mat(pos1, pos2);
        item_b = sum(btw_mat,2) / size(btw_mat,2);
        full_itemb(pos1,1) = item_b;
        
        rdm_sm.con{i}.catg.between(j,1) = mean(item_b);
        
        rdm_sm.con{i}.catg.info(j,1) = ...
            rdm_sm.con{i}.catg.between(j,1) - rdm_sm.con{i}.catg.within(j,1);
    end
    rdm_sm.con{i}.item.within = full_itemw;
    rdm_sm.con{i}.item.between = full_itemb;
    rdm_sm.con{i}.item.info = ...
        rdm_sm.con{i}.item.between - rdm_sm.con{i}.item.within;
    
    % Get average across all categories
    rdm_sm.con{i}.avg.within = mean(rdm_sm.con{i}.catg.within);
    rdm_sm.con{i}.avg.between = mean(rdm_sm.con{i}.catg.between);
    rdm_sm.con{i}.avg.info = rdm_sm.con{i}.avg.between - rdm_sm.con{i}.avg.within;
end



