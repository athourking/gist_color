function C2res_so = jun_hmax_gist_8scenes(numTrain,mode)

tmpDir = sprintf('/gpfs/data/jz7/ColorGist/results/%i_train/1_split',numTrain);
train = load(fullfile(tmpDir,sprintf('train_%i_train_1_split.mat', ...
             numTrain)));
train = train.train;
test = load(fullfile(tmpDir,sprintf('test_%i_train_1_split.mat', ...
             numTrain)));         
test = test.test;
          
datapath = '/users/jz7/data/ColorGist/images/spatial_envelope_256x256_static_8outdoorcategories';

cI = jun_readAllImages_8scenes(train,test,datapath);       

patchSizes = [4 8 12 16];
numPatchSizes = length(patchSizes);
numPatchesPerSize = 250;
numChannel = 8;
Nfeatures = numPatchesPerSize*numPatchSizes;


switch mode
    case 0
    numPhases = 2;
    outDir = sprintf('/gpfs/data/jz7/ColorGist/results/%i_train/C2/%i_phases',numTrain,numPhases);

    cPatches_so = load(fullfile(outDir,sprintf('dict_c1so_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases, numPatchesPerSize, length(patchSizes))));
    cPatches_do = load(fullfile(outDir,sprintf('dict_c1do_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases, numPatchesPerSize, length(patchSizes))));      
    cPatches_gray = load(fullfile(outDir,sprintf('dict_c1gray_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases, numPatchesPerSize, length(patchSizes))));
    case 1
    numPhases = 1; 
    outDir = sprintf('/gpfs/data/jz7/ColorGist/results/%i_train/C2/%i_phases',numTrain,numPhases);

    cPatches_so = load(fullfile(outDir,sprintf('dict_c1so_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases, numPatchesPerSize, length(patchSizes))));
    cPatches_do = load(fullfile(outDir,sprintf('dict_c1do_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases, numPatchesPerSize, length(patchSizes))));      
    cPatches_gray = load(fullfile(outDir,sprintf('dict_c1gray_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases, numPatchesPerSize, length(patchSizes))));

    case 2
    numPhases = 4;   
    outDir = sprintf('/gpfs/data/jz7/ColorGist/results/%i_train/C2/%i_phases',numTrain,numPhases);

    cPatches_so = load(fullfile(outDir,sprintf('dict_c1so_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases, numPatchesPerSize, length(patchSizes))));
    cPatches_do = load(fullfile(outDir,sprintf('dict_c1do_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases, numPatchesPerSize, length(patchSizes))));      
    cPatches_gray = load(fullfile(outDir,sprintf('dict_c1gray_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases, numPatchesPerSize, length(patchSizes))));  
end


    cPatches_so = cPatches_so.cPatches;   
    cPatches_do = cPatches_do.cPatches;
    cPatches_gray = cPatches_gray.cPatches;
        
% % The actual C2 features are computed below for each one of the training/testing directories
tic
for i = 1:2,
%     i = 1;
    C2res_so{i} = extractC2forcell_so_mode(cPatches_so,cI{i},numPatchSizes,numChannel, numTrain,outDir,i,mode);
% 
%     [C2res_do{i},C2res_gray{i}] = extractC2forcell_dogray_mode(cPatches_gray,cPatches_do,cI{i},numPatchSizes,numChannel,outDir, i,mode);
       toc  
    
end
totaltimespectextractingC2 = toc;

save(fullfile(outDir,sprintf('c2so1_%i_phase_%i_patches_%i_sizes.mat', ...
             numPhases,numPatchesPerSize, length(patchSizes))), 'C2res_so','-v7.3');

% save(fullfile(outDir,sprintf('c2do_%i_phase_%i_patches_%i_sizes.mat', ...
%              numPhases,numPatchesPerSize, length(patchSizes))), 'C2res_do','-v7.3');
% 
% save(fullfile(outDir,sprintf('c2gray_%i_phase_%i_patches_%i_sizes.mat', ...
%               numPhases,numPatchesPerSize, length(patchSizes))), 'C2res_gray','-v7.3');
        
return



