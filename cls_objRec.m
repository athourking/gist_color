%
% This script trains a libSVM for classification:
% available at: http://www.csie.ntu.edu.tw/~cjlin/libsvm/

addpath(genpath('/gpfs/data/jz7/libsvm/libsvm-3.1/libsvm-3.1'));
datapath = '/gpfs/data/tserre/jzhang/general_color_image_processing_functions/Joost_Weijer/ColorDescriptors/soccer_data/soccer';

numPhases = 1;
numTest = 15;
numTrain = 25;
patchSizes = [4 8 12 16];
numPatchesPerSize = 250;


dirs = dir(datapath);
catNames = {dirs([dirs.isdir]).name};
catNames = setdiff(catNames, {'.', '..'});
[ans, inds] = sort(lower(catNames));
catNames = catNames(inds);
Nclasses = numel(catNames);


% training and testing labels
yTrain = [];
yTest = [];
for ii = 1:Nclasses
    yTrain =  [yTrain;ii*ones(numTrain,1)];
    yTest =  [yTest;ii*ones(numTest,1)];
end


for split = 1:3
    
    outDir = sprintf('/gpfs/data/tserre/jzhang/object_recognition/results/soccer/%i_split',split);
    
    c2so = load(fullfile(outDir,sprintf('eccv_c2so_%i_orients_%i_split_%i_patches_%i_sizes.mat', ...
        2, split, numPatchesPerSize, length(patchSizes))));
    c2so = c2so.C2res;
    
    c2do = load(fullfile(outDir,sprintf('eccv_c2do_%i_orients_%i_split_%i_patches_%i_sizes.mat', ...
        2, split,numPatchesPerSize, length(patchSizes))));
    c2do = c2do.C2res;
    
    c2gray = load(fullfile(outDir,sprintf('eccv_c2_%i_orients_%i_split_%i_patches_%i_sizes.mat', ...
        2 ,split,numPatchesPerSize, length(patchSizes))));
    c2gray = c2gray.C2res;
    
    
    
    XTrain.gray = [c2gray{1}]; XTest.gray =[c2gray{2}];
    XTrain.do = [c2do{1}]; XTest.do =[c2do{2}];
    XTrain.so = [c2so{1}]; XTest.so =[c2so{2}];

    XTrain.s1 = zscore(XTrain.so',0,1); XTest.s1 = zscore(XTest.so',0,1);
    XTrain.c1 = zscore(XTrain.do',0,1); XTest.c1 = zscore(XTest.do',0,1);
    XTrain.gray = zscore(XTrain.gray',0,1);XTest.gray = zscore(XTest.gray',0,1);
    
    XTrain.col = [XTrain.s1,XTrain.c1]; XTest.col =[XTest.s1,XTest.c1];
    XTrain.all = [XTrain.s1,XTrain.c1,XTrain.gray]; XTest.all =[XTest.s1,XTest.c1,XTest.gray];
    XTrain.shape = [XTrain.c1, XTrain.gray]; XTest.shape =[XTest.c1, XTest.gray];
    
    model.s1 = svmtrain_libsvm(yTrain, XTrain.s1, '-t 0 -c 1000 -b 1');
    model.c1 = svmtrain_libsvm(yTrain, XTrain.c1, '-t 0 -c 1000 -b 1');
    model.gray = svmtrain_libsvm(yTrain, XTrain.gray, '-t 0 -c 1000 -b 1');
    model.all = svmtrain_libsvm(yTrain, XTrain.all, '-t 0 -c 1000 -b 1');
    model.shape = svmtrain_libsvm(yTrain, XTrain.shape, '-t 0 -c 1000 -b 1');
    model.col = svmtrain_libsvm(yTrain, XTrain.col, '-t 0 -c 1000 -b 1');
    
    [~, accuracy.s1(:,split), pb.s1] = svmpredict(yTest, XTest.s1, model.s1,'-b 1');
    [~, accuracy.c1(:,split), pb.c1] = svmpredict(yTest, XTest.c1, model.c1,'-b 1');
    [~, accuracy.gray(:,split), pb.gray] = svmpredict(yTest, XTest.gray, model.gray,'-b 1');
    [~, accuracy.shape(:,split), ~] = svmpredict(yTest, XTest.shape, model.shape,'-b 1');
    [~, accuracy.all(:,split), ~] = svmpredict(yTest, XTest.all, model.all,'-b 1');
    [~, accuracy.col(:,split), ~] = svmpredict(yTest, XTest.col, model.col,'-b 1');
    
    
    
    % probability multiplication
    pb.all =  pb.s1.* pb.c1.* pb.gray;
    pb.shape =  pb.c1.* pb.gray;
    pb.col =  pb.s1.* pb.c1;
    
    
    for k = 1:length(yTest)
        [foo.all,predict_label.all(k)] = max(pb.all(k,:));
        [foo.shape,predict_label.shape(k)] = max(pb.shape(k,:));
        [foo.sodo,predict_label.col(k)] = max(pb.col(k,:));
    end
    
    for j = 1:Nclasses
        for i = 1:Nclasses
            Cm.all(i,j,split) = 100*sum((yTest==i) .* (predict_label.all'==j))/sum(yTest==i);
            Cm.shape(i,j,split) = 100*sum((yTest==i) .* (predict_label.shape'==j))/sum(yTest==i);
            Cm.col(i,j,split) = 100*sum((yTest==i) .* (predict_label.col'==j))/sum(yTest==i);
        end
    end
    
    
    perf.all(split) = mean(diag(Cm.all(:,:,split)));
    perf.shape(split) = mean(diag(Cm.shape(:,:,split)));
    perf.col(split) = mean(diag(Cm.col(:,:,split)));
    
end


