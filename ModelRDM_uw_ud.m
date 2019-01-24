function [mod_rdm] = ModelRDM_uw_ud(cond,edge)
%% --------------------- Script Description -----------------------------
% Script for creating a model dissimilarity matrix (unweighted,undirected)
% that can be used as a predictor for comparison to other RDM. 
% Model created will assume connections ('1') within each row listed as 
% being in the same condition and no connections ('0') between the
% different conditions. 
% 
% To manually define connections, input edge can be use to add connections.
% 
% Takes in the with the following inputs:
% 1) cond - datate structure with the following fields: 
%     i) names - 1*M cell string array with M being the no. conditions
%    ii) vect - 1*M cell array each containing a vector for rows belonging 
%               to the respective conditions.
%   iii) matsize - Integer to indicate size of the RDM (square matrix)
% 
% 2) edge - N*2 matrix with N being the number of additional edges to
%           connect. Values within each row correspond to Row and Column to
%           define a connection.
%
% Created by JH 22/9/2018
%
%% 
if nargin < 2
    edge = [];
end
mod_rdm = struct();
matsize = cond.matsize;

d = ones(1,matsize);
mat = diag(d);

labels = cell(matsize,1);

for i = 1 : length(cond.names)
    v = cond.vect{i};
    
    % Get all within-condition edges
    [R, C] = meshgrid(v,v);
    tmp = cat(2,C,R);
    pos = reshape(tmp,[],2);
    
    mat(pos(:,1),pos(:,2)) = 1;
    labels(pos(:,1)) = cond.names(i);
end

for i = 1 : length(edge)
    pos = edge(i,:);
    
    mat(pos(:,1),pos(:,2)) = 1;
end

mod_rdm.labels = labels;
mod_rdm.mat = mat;

% Visualise rdm
imagesc(mat);
colormap(flipud(gray));

