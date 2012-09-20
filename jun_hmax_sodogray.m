function [C2res_so,C2res_do, C2res_gray] = jun_hmax_sodogray(numTrain,numPhases)


tmpDir = sprintf('../results/101/hmax/%i_train',numTrain);
load(fullfile(tmpDir,sprintf('train_%i_train_cal101.mat', ...
             numTrain)));
load(fullfile(tmpDir,sprintf('test_%i_train_cal101.mat', ...
              numTrain)));        

outDir = sprintf('../results/101/hmax/%i_train/%i_phase',numTrain,numPhases);


READPATCHESFROMFILE = 1;
patchSizes = [4 8 12 16];
numPatchSizes = length(patchSizes);
numPatchesPerSize = 250;
numChannel = 8;
Nfeatures = numPatchesPerSize*numPatchSizes;


cI = jun_readAllImages_multiclass(train,test);
clear train;
clear test;

if isempty(cI{1})
    error(['No training images were loaded -- did you remember to' ...
	       ' change the path names?']);
end     
%       cItrain = [cI{1};cI{2}];
      %below the c1 prototypes are extracted from the images/ read from file
if ~READPATCHESFROMFILE
    %take more time to compute
    cPatches = extractRandC1Patches_normPhase(cI{1}, numPatchSizes, ...
            numPatchesPerSize, patchSizes,numChannel,numPhases); %fix: extracting from positive only
    
    save(fullfile(outDir,sprintf('dict_c1_%i_numTrain_%i_patches_%i_sizes.mat', ...
             numTrain, numPatchesPerSize, length(patchSizes))) ,'cPatches','-v7.3');
    
else
    fprintf('reading patches');
    cPatches_gray = load(fullfile(outDir,sprintf('dict_c1gray_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases, numPatchesPerSize, length(patchSizes))));
    cPatches_so = load(fullfile(outDir,sprintf('dict_c1so_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases,numPatchesPerSize, length(patchSizes))));      
    cPatches_do = load(fullfile(outDir,sprintf('dict_c1do_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases, numPatchesPerSize, length(patchSizes))));
    
    cPatches_so = cPatches_so.cPatches_so;
    cPatches_do = cPatches_do.cPatches_do;
    cPatches_gray = cPatches_gray.cPatches_gray;
end
       

%----Settings for Testing --------%
rot = [90 -45 0 45];
c1ScaleSS = [1:2:18];
RF_siz    = [7:2:39];
c1SpaceSS = [8:2:22];
minFS     = 7;
maxFS     = 39;
div = [4:-.05:3.2];
Div       = div;
%      %--- END Settings for Testing --------%
      
fprintf(1,'Initializing gabor filters -- full set...');
%creates the gabor filters use to extract the S1 layer
[fSiz,filter, filters,c1OL,numSimpleFilters] = init_color_gabor_phases(rot, RF_siz, Div,numChannel,numPhases);
[fSiz_gray,filters_gray,c1OL_gray,numSimpleFilters_gray] = init_gabor(rot, RF_siz, Div); 
fprintf(1,'done\n');

             
% % The actual C2 features are computed below for each one of the training/testing directories
tic
for i = 1:2,
%     C2res_so{i} = extractC2forcell_so(filters,fSiz,c1SpaceSS,c1ScaleSS,...
%         c1OL,cPatches_so,cI{i},numPatchSizes,numChannel,numPhases, numTrain,outDir,i);
% 
    [C2res_do{i},C2res_so{i},C2res_gray{i}] = extractC2forcell_sodogray_new(filters_gray,fSiz_gray,c1OL_gray,cPatches_gray,filter, filters,fSiz,c1SpaceSS,c1ScaleSS,...
        c1OL,cPatches_so,cPatches_do,cI{i},numPatchSizes,numChannel,numPhases, numTrain,outDir,i);
   
%     C2res_gray{i} = extractC2forcell_gray(filters_gray,fSiz_gray,c1SpaceSS,c1ScaleSS,...
%         c1OL_gray,cPatches_gray,cI{i},numPatchSizes, numTrain,outDir,i);
    toc  
    
end
totaltimespectextractingC2 = toc;



save(fullfile(outDir,sprintf('c2so_%i_numTrain_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain, numPatchesPerSize, length(patchSizes))), 'C2res_so','-v7.3');

save(fullfile(outDir,sprintf('c2do_%i_numTrain_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain, numPatchesPerSize, length(patchSizes))), 'C2res_do','-v7.3');

save(fullfile(outDir,sprintf('c2gray_%i_numTrain_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain, numPatchesPerSize, length(patchSizes))), 'C2res_gray','-v7.3');
        
return



