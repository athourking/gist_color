function [cPatches_so,cPatches_do,cPatches_gray] = jun_c1Patches_sodogray_8scenes(numTrain,mode)


        
tmpDir = sprintf('/gpfs/data/jz7/ColorGist/results/%i_train',numTrain);
train = load(fullfile(tmpDir,sprintf('train_%i_train.mat', ...
             numTrain)));
train = train.train;
test = load(fullfile(tmpDir,sprintf('test_%i_train.mat', ...
             numTrain)));         
test = test.test;
          
datapath = '/users/jz7/data/ColorGist/images/spatial_envelope_256x256_static_8outdoorcategories';

cI = jun_readAllImages_8scenes(train,test,datapath);       
if isempty(cI)
    error(['No training images were loaded -- did you remember to' ...
	       ' change the path names?']);
end     
%       cItrain = [cI{1};cI{2}];





%%
READPATCHESFROMFILE = 0;
patchSizes = [4 8 12 16];
numPatchSizes = length(patchSizes);
numPatchesPerSize = 250;
numChannel = 8;
Nfeatures = numPatchesPerSize*numPatchSizes;


      %below the c1 prototypes are extracted from the images/ read from file
if ~READPATCHESFROMFILE
    %take more time to compute
    cPatches_so = extractRandC1Patches_so_mode(cI{1}, numPatchSizes, ...
            numPatchesPerSize, patchSizes,numChannel,numTrain,mode); %fix: extracting from positive only
%     save(fullfile(outDir,sprintf('dict_c1so_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
%              numTrain,numPhases, numPatchesPerSize, length(patchSizes))) ,'cPatches_so','-v7.3');
   
    cPatches_do = extractRandC1Patches_normPhase_mode(cI{1}, numPatchSizes, ...
            numPatchesPerSize, patchSizes,numChannel,numTrain,mode); %fix: extracting from positive only 
%     save(fullfile(outDir,sprintf('dict_c1do_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
%             numTrain,numPhases, numPatchesPerSize, length(patchSizes))) ,'cPatches_do','-v7.3');
    
    cPatches_gray = extractRandC1Patches_gray_mode(cI{1}, numPatchSizes, ...
        numPatchesPerSize, patchSizes,numTrain,mode); 
%     save(fullfile(outDir,sprintf('dict_c1gray_%i_train_%i_patches_%i_sizes.mat', ...
%              numTrain, numPatchesPerSize, length(patchSizes))) ,'cPatches_gray','-v7.3');
 
    
else
    fprintf('reading patches');
    cPatches_so = load(sprintf('dict_c1_%i_split_%i_patches_%i_sizes.mat', ...
            split, numPatchesPerSize, length(patchSizes)) ,'cPatches_so');  
    cPatches_do = load(sprintf('dict_c1_%i_split_%i_patches_%i_sizes.mat', ...
            split, numPatchesPerSize, length(patchSizes)) ,'cPatches_do'); 
    cPatches_gray = load(sprintf('dict_c1_%i_split_%i_patches_%i_sizes.mat', ...
            split, numPatchesPerSize, length(patchSizes)) ,'cPatches_gray');

    cPatches_so = cPatches_do.cPatches_so;   
    cPatches_do = cPatches_do.cPatches_do;
    cPatches_gray = cPatches_gray.cPatches_gray;
    
end
   

         
return



