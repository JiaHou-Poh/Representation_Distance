% Create data structure for Pattern analysis with Cross validation
% Requires structure ds with the following fields:
% ds - data structure with the following fields:
%     i) runs: Column vector indicating which run each trial belongs to
%    ii) cond: Column vector indicating which cond each trial belongs to
%   iii) vx: M*N matrix of values (M = items, N = no. voxels/values)
%    iv) fold: M*N matrix assigning each run to train/test (M = no. items
%              N = no. iteration) ; 1 - Test, 2 - Train, 0 - Unused
%     v) condname : String array indicating name of each condition
%
% This version also includes selection for Top X% of voxels

clear all
clear mex

mainFldr = '/Volumes/Tera2b/Experiments/ACTS2_MVPA/';
memFldr = [mainFldr 'Analysis/MemSplit/'];
roiFldr = [mainFldr,'Analysis/ROI_Output/'];

roiname = 'PPA';
memtype = {'all','hrem','lrem','rem','forg'};
%memtype={'all'};
imgtype = {'spmT'};

distmeasure = 'correlation';
featsel = 85; % Select the most strongly activated percentile of voxels from the 'All' condition


saveFldr = [mainFldr, 'Analysis/ROI_RSA/' roiname '/'];
if ~exist(saveFldr)
    mkdir(saveFldr)
end

numfolds = 4;

subjects =[4];
%subjects = [4 5 6 9 10 11 12 13 14 15 16 18 19 21 22 23 24 26 27 28 29 30 31 32:37 39:44]; %
sessions = {'Stim','Sham'};

for i = 1 : length(subjects)
    disp(sprintf('Processing subject %d\n',subjects(i)));
    subjno = subjects(i); % numeric
    id = sprintf('ACTS2%03d', subjno);
 
    for s = 1 : length(sessions)
        ses = sessions{s};
        rsa = struct();
        
        % Load the roi mat
        roimat = [roiFldr, roiname '/' id '_' ses '_' imgtype{1}, '_ROI.mat'];
        r = load(roimat);
        r = r.roi;
        
        inp = r.input;
        
        num_roi = length(r.(roiname));
        feat = struct();
        
        % Define threshold separatly for each session
        cutoff = zeros(1,num_roi);
        
        for nr = 1 : num_roi
            w_out = r.(roiname)(nr);
        
            for t = 1 : length(memtype)
                ds = struct();

                cond = inp.(memtype{t}).catg;
                %cond = inp.(memtype{t}).subcatg;
                ds.cond(:,1) = cond;
                ds.numcond = max(cond);

                %Check that length of the mat matches the output
                out = w_out.(memtype{t}).output;

                if length(out) ~= length(cond)
                    error('number of trials mismatch')
                end

                % Get the voxel values for each of the images
                for n = 1 : length(out)
                    ds.val(n,:) = out(n).roi_mean;
                    ds.vx(n,:) = out(n).roi_val;
                    ds.vxz(n,:) = out(n).roi_zscore;
                    ds.snr(n,:) = out(n).spatialsnr;
                end
                
                % Find the top X percentile of voxels from ds
                ds.vxmean = mean(ds.vx,1);
                if size(ds.vxmean,1)~=1
                    error('vector length incorrect')
                end
                
                if strcmp (memtype{t}, 'all')
                    cutoff(1,nr) = prctile(ds.vxmean,featsel);
                    feat.roi(nr).idx = find(ds.vxmean >= cutoff(1,nr));
                end
                
                % Keep relevant columns
                ds.vx = ds.vx(:,feat.roi(nr).idx);
                ds.vxz = ds.vxz(:,feat.roi(nr).idx);
                ds.val = mean(mean(ds.vx,2));
                ds.snr(:) = NaN;
                

                folds = zeros(length(cond),numfolds);
                testcount = zeros(numfolds,ds.numcond);
                traincount = zeros(numfolds,ds.numcond);

                wrun = inp.(memtype{t}).idx(:,1);
                ds.runs(:,1) = wrun;

                % Define the folds
                for f = 1 : numfolds
                    test_idx = find(wrun == f);
                    train_idx = find(wrun ~= f);
                    folds(test_idx,f) = 1;
                    folds(train_idx,f) = 2;

                    % Count number of items from each category in each fold
                    for c = 1 : ds.numcond
                        testcount(f,c) = length(find(cond(test_idx) == c));
                        traincount(f,c) = length(find(cond(train_idx) == c));
                    end
                end
                
                ds.fold = folds;
                ds.testcount = testcount;
                ds.traincount = traincount;

                % Run the Cross val function
                [cv] = CrossVal_Dist(ds,distmeasure,0);
                %[cv] = CrossVal_Dist(ds,'euclidean',0);
                %[cv] = CrossVal_Dist(ds,'spearman',0);
                
                % Assign to structure
                rsa.roi(nr).(memtype{t}).ds = ds;
                rsa.roi(nr).(memtype{t}).results = cv;
            end
        end
        % Save mat separately for each session
        savemat = [saveFldr, id ,'_' ses '_' imgtype{1}, '_TopVxPct' num2str(featsel) '_ROI.mat'];
        save(savemat,'rsa')
        disp(sprintf('Completed subject %d Session %s \n',subjects(i), ses));
    end
end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
