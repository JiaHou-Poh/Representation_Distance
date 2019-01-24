function [pd] = Pairwise_Dist(ds,dist)
%% --------------------- Script Description -----------------------------
% Script for doing pairwise similarity using dist metric. Compares items
% from the same condition. If there are >2 items belonging to the same
% condition, final output will be the averaged pairwise correlations
% between each pair (i.e. note for comparing each item to the condition
% average, use CrossVal_Dist)). For across-run repetitions, a separate
% distance value will also be generated for each within-run pair.
%
% Takes in the following input:
% 1) ds - data structure with the following fields
%     i) runs: Column vector indicating which run each trial belongs to
%    ii) cond: Column vector indicating which cond each trial belongs to
%   iii) vx: M*N matrix of values (M = items, N = no. voxels/values)
%    iv) condname: String array for name of each condition(optional)
%
% 2) dist - distance metric to use. Takes the following string -
%           'correlation','spearman','euclidean' etc (refer to pdist2.m)
%
% A single output structure would be generated with the following fields:
% 1) pd - Contains the average distance for all conditions. Includes:
%       i) within_run : vector containing the average distance between
%                items of the same condition within the same run
%      ii) run_avg : single value containing the average of all the within
%                   run distances
%     iii) all_runs : single value containing the average of all pairwise
%                    comparisons regardless of run
%
% Completed by JH 23/3/2018
% Version history
% Edits 26/3/2018
% Modified the within run correlation to input NaN if a specific condition
% doesnt exist within a specific run

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
pd = struct();
pd.dist = dist;

numCond = max(ds.cond);
numItem = size(ds.vx,1);
numVx = size(ds.vx,2);
numRun = max(ds.runs);

% Identify conditions
con = struct();
for i = 1 : numCond
    con(i).pos = find(ds.cond == i);
end

% Create similarity measure across all runs
for i = 1 : numCond
    con(i).fullmat = ds.vx(con(i).pos,:);
    full_rdm = squareform(pdist(con(i).fullmat,dist));
    
    % Take the mean without the off-diagonal - matrix is symmetrical
    idx = ones(size(full_rdm));
    ut_idx = triu(idx,1);
    all_corr = full_rdm(ut_idx==1);
    pd.all_runs(i,1) = mean(all_corr,1);
end

% Create similarity measure only within run
for i = 1 : numCond
    for j = 1 : numRun
        run_idx = find(ds.runs == j);
        cond_run_idx = intersect(con(i).pos, run_idx);
        
        if isempty(cond_run_idx)
            run_corr = NaN;
        else        
            con(i).run(j).fullmat = ds.vx(cond_run_idx,:);
            full_run_rdm = squareform(pdist(con(i).run(j).fullmat,dist));
            
            % Take the mean without the off-diagonal - matrix is symmetrical
            idx = ones(size(full_run_rdm));
            ut_idx = triu(idx,1);
            run_corr = full_run_rdm(ut_idx==1);
        end
        pd.within_run(i,j) = run_corr;
    end
end
pd.run_avg = nanmean(pd.within_run,2);

end
        
        
        
        
        
        
        
        
