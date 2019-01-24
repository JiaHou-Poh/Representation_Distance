function [cv] = CrossVal_Dist(ds,dist,perms)
%% --------------------- Script Description -----------------------------
% Script for doing cross validation using dist metric. Creates an average
% template pattern using the training blocks and 'test' using the left out
% run. Testing would be done by subtracting distance of within to between.
% If there are >2 categories, average difference would be taken. If there
% are multiple test-items in each fold, distance will be computed on the
% pattern for individual items, and also on the pattern for the averaged
% pattern. 
% * Note: it is assumed that each item is tested only once.
%
% Takes in the following input: 
% 1) ds - data structure with the following fields:
%     i) runs: Column vector indicating which run each trial belongs to
%    ii) cond: Column vector indicating which cond each trial belongs to
%   iii) vx: M*N matrix of values (M = items, N = no. voxels/values)
%    iv) fold: M*N matrix assigning each run to train/test (M = no. items
%              N = no. iteration) ; 1 - Test, 2 - Train, 0 - Unused
%     v) condname : String array indicating name of each condition
% 
% 2) dist - distance metric to use. Takes the following string -
%           'correlation','spearman','euclidean' etc (refer to pdist2.m)
% 
% 3) perm - takes in an integer to indicate if permutation testing is done
%           0 to indicate no permutation test *To be added*
%
% A single output structure would be generated with the following fields:
% 1) cv.catg - which contains the category level analysis. Includes the
%             following:
%       i) within : M*N matrix showing the within category distances
%      ii) between : M*N matrix showing the between category distances
%     iii) info : M*N matrix showing Between - Within (distance measures,
%                   further means less similar)
%          ---- M corresponds to number of CV folds and N corresponds to
%               each category
%      iv) within_avg : 1*N vector showing average within category
%                       distances across all folds
%       v) between_avg: 1*N vector showing average between category
%                       distances across all folds
%      vi) info_avg: 1*N vector showing the difference of the averages
%                     across folds
%
%  2) cv.avg - Contain the average across all categories and all folds
%       i) within : Single value for within category distance across all
%                   categories and fold
%      ii) between: Single value for between category distance across all
%                   categories and fold
%     iii) info : Single value for the difference across all categories and
%                 fold
%  
%  3) cv.trial - Contains the trial level analysis. Includes the following:
%       i) within : single value containing the trial's distance to the
%                   average pattern of it's own category
%      ii) between: vector containing the trial's distance to the
%                   average pattern of each of the other categories
%     iii) info : vector containing the difference of between and within
%                 for each category
%      iv) info_avg : single value derived from the average of info
%
% Completed by JH - 7/12/2017
% Functions to include:
% 1) Generate null distribution for permutation testing
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
cv = struct();
cv.dist = dist;

numFold = size(ds.fold,2);
numCond = max(ds.cond);
numItem = size(ds.vx,1);
numVx = size(ds.vx,2);

cv.perm = perms;

% Identify conditions
con = struct();
for i = 1 : numCond
    con(i).pos = find(ds.cond == i);
end

% Check if items tested only once
chkTest = sum(ds.fold == 1,2);
% if length(unique(chkTest)) ~= 1
%     error('Each item must belong to ''test'' equal number of times')
% end
if max(chkTest)>1
    error('Each item can only belong to ''test'' once');
end

for i = 1 : numFold
    con_cv = struct();
    
    % Identify trials for testing & training
    cv_vect = ds.fold(:,i);
    test_trials = find(cv_vect == 1);
    train_trials = find(cv_vect == 2);
    
    avg_test_mat = zeros(numCond,numVx);
    avg_train_mat = zeros(numCond,numVx);
                
            
    % Identify trial number for each testing and training condition
    for j = 1 : numCond
        con_cv(j).testIdx = intersect(con(j).pos,test_trials);
        con_cv(j).trainIdx = intersect(con(j).pos,train_trials);
        
        % Extract value for each condition
        % Note that using pdist2, the distance is computed between rows
        % Each row is considered as a 'point'
        con_cv(j).testvx = ds.vx(con_cv(j).testIdx,:);
        con_cv(j).trainvx = ds.vx(con_cv(j).trainIdx,:);
        
        con_cv(j).testvx =  con_cv(j).testvx;
        con_cv(j).trainvx = con_cv(j).trainvx;
        
        % Create average pattern for each condition
        con_cv(j).testavg = mean(con_cv(j).testvx,1);
        con_cv(j).trainavg = mean(con_cv(j).trainvx,1);
        
        avg_test_mat(j,:) = con_cv(j).testavg;
        avg_train_mat(j,:) = con_cv(j).trainavg; 
    end
    
   
    % Do the comparison for each fold. Create confusion matrix for each
    % pair of vector. Diagonals will be the within category for each
    % category, and the off-diagonals will be the between. Matrix is not
    % symmetrical. In output row corresponds to training category, column
    % corresponds to testing category
    conf_mat = pdist2(avg_train_mat,avg_test_mat,dist);
    cv.catg.within(i,:) = diag(conf_mat);
    
    off_diag = conf_mat .* ~eye(size(conf_mat));
    cv.catg.between(i,:) = sum(off_diag,1) / (size(off_diag,1)-1); %mean without diagonal
    
    % This is distance measure, so between should be greater than within
    cv.catg.info(i,:) = cv.catg.between(i,:) - cv.catg.within(i,:);
    
    % Do the itemised correlation within and between category
    % Assumes each item only tested once
    all_col = 1 : numCond;
    
    for j = 1 : numCond
        trial_test_mat = con_cv(j).testvx;
        orig_idx = con_cv(j).testIdx; % the original index position of trial
        
        % Over here, output row corresponds to test items and column
        % coresponds to train category
        conf_mat = pdist2(trial_test_mat,avg_train_mat,dist);
        
        % Get within category
        trial_within = conf_mat(:,j); 
        cv.trial.within(orig_idx,1) = trial_within;
          
        % Get between category
        btw_idx = all_col(all_col~=j);
        
        trial_between = conf_mat(:,btw_idx);
        cv.trial.between(orig_idx,:) = trial_between;
        
        % Get difference
        cv.trial.info(orig_idx,:) = trial_between - trial_within;
        cv.trial.info_avg(orig_idx,1) = mean(cv.trial.info(orig_idx,:),2);
        
    end
    % Compute item distinctiveness measure - Quantified as distance
    % within weighted by the difference of between - within). If item
    % is not more similar to it's own category (dist within > dist
    % between), item would have negative weight.
    
    cv.trial.weight_distinct = ...
        cv.trial.within .* cv.trial.info;
end
cv.catg.within_avg = mean(cv.catg.within,1);
cv.catg.between_avg = mean(cv.catg.between,1);
cv.catg.info_avg = cv.catg.between_avg - cv.catg.within_avg;
cv.catg.info_avg_all = mean(cv.catg.info_avg,2);

cv.avg.within = mean(cv.catg.within_avg,2);
cv.avg.between = mean(cv.catg.between_avg,2);
cv.avg.info = cv.avg.between - cv.avg.within;


% Generate null distribution using permutation testing -- TO ADD ---
end
        
        
        
    
        
    
    
    
    
    