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
for ii = 1:7
    yTrain =  [yTrain;ii*ones(numTrain,1)];
    yTest =  [yTest;ii*ones(numTest,1)];
end


%% load C2 features
split = 2;
outDir = sprintf('/gpfs/data/tserre/jzhang/object_recognition/results/soccer/%i_split',split);


C2res = load(fullfile(outDir,sprintf('eccv_c2so_%i_orients_%i_split_%i_patches_%i_sizes.mat', ...
    2, split, numPatchesPerSize, length(patchSizes))));
C2res = C2res.C2res;
c2.so = [];
for ii = 1:7
    c2.so = [c2.so,C2res{1}(:,(ii-1)*numTrain+1:(ii-1)*numTrain+numTrain)];
    c2.so = [c2.so,C2res{2}(:,(ii-1)*numTest+1:(ii-1)*numTest+numTest)];
end
        
C2res = load(fullfile(outDir,sprintf('eccv_c2do_%i_orients_%i_split_%i_patches_%i_sizes.mat', ...
    2, split,numPatchesPerSize, length(patchSizes))));
C2res = C2res.C2res;
c2.do = [];
for ii = 1:7
    c2.do = [c2.do,C2res{1}(:,(ii-1)*numTrain+1:(ii-1)*numTrain+numTrain)];
    c2.do = [c2.do,C2res{2}(:,(ii-1)*numTest+1:(ii-1)*numTest+numTest)];
end

C2res = load(fullfile(outDir,sprintf('eccv_c2_%i_orients_%i_split_%i_patches_%i_sizes.mat', ...
    2 ,split,numPatchesPerSize, length(patchSizes))));
C2res = C2res.C2res;
c2.gray = [];
for ii = 1:7
    c2.gray = [c2.gray,C2res{1}(:,(ii-1)*numTrain+1:(ii-1)*numTrain+numTrain)];
    c2.gray = [c2.gray,C2res{2}(:,(ii-1)*numTest+1:(ii-1)*numTest+numTest)];
end


outDir = sprintf('/gpfs/data/jz7/Object_recognition/results/101/hmax/soccer/1_phase');

datasplits = load(fullfile(outDir,'datasplits.mat'));
datasplits = datasplits.datasplits;


%% classify
for split = 1:3
        
    tr = datasplits{split,1};
    ts = datasplits{split,2};
    trn = []; tst = [];
    for i = 1:7
        trn = [trn,40*(i-1)+tr(1+25*(i-1):25+25*(i-1))];
        tst = [tst,40*(i-1)+ts(1+15*(i-1):15+15*(i-1))];
    end
    
    
    XTrain.gray = [c2.gray(:,trn)]; XTest.gray =[c2.gray(:,tst)];
    XTrain.do = [c2.do(:,trn)]; XTest.do =[c2.do(:,tst)];
    XTrain.so = [c2.so(:,trn)]; XTest.so =[c2.so(:,tst)];

    
    XTrain.so = zscore(XTrain.so',0,1);
    XTest.so = zscore(XTest.so',0,1);
    XTrain.do = zscore(XTrain.do',0,1);
    XTest.do = zscore(XTest.do',0,1);
    XTrain.gray = zscore(XTrain.gray',0,1);
    XTest.gray = zscore(XTest.gray',0,1);
    
    XTrain.col = [XTrain.so,XTrain.do]; XTest.col =[XTest.so,XTest.do];
    XTrain.all = [XTrain.so,XTrain.do,XTrain.gray]; XTest.all =[XTest.so,XTest.do,XTest.gray];
    XTrain.shape = [XTrain.do, XTrain.gray]; XTest.shape =[XTest.do, XTest.gray];
    
    model.so = svmtrain_libsvm(yTrain, XTrain.so, '-t 0 -c 1000 -b 1');
    model.do = svmtrain_libsvm(yTrain, XTrain.do, '-t 0 -c 1000 -b 1');
    model.gray = svmtrain_libsvm(yTrain, XTrain.gray, '-t 0 -c 1000 -b 1');
    model.all = svmtrain_libsvm(yTrain, XTrain.all, '-t 0 -c 1000 -b 1');
    model.shape = svmtrain_libsvm(yTrain, XTrain.shape, '-t 0 -c 1000 -b 1');
    model.col = svmtrain_libsvm(yTrain, XTrain.col, '-t 0 -c 1000 -b 1');
    
    [~, accuracy.so(:,split), pb.so] = svmpredict(yTest, XTest.so, model.so,'-b 1');
    [~, accuracy.do(:,split), pb.do] = svmpredict(yTest, XTest.do, model.do,'-b 1');
    [~, accuracy.gray(:,split), pb.gray] = svmpredict(yTest, XTest.gray, model.gray,'-b 1');
    [~, accuracy.shape(:,split), ~] = svmpredict(yTest, XTest.shape, model.shape,'-b 1');
    [~, accuracy.all(:,split), ~] = svmpredict(yTest, XTest.all, model.all,'-b 1');
    [~, accuracy.col(:,split), ~] = svmpredict(yTest, XTest.col, model.col,'-b 1');
    
    
    
    % probability multiplication
    pb.all =  pb.so.* pb.do.* pb.gray;
    pb.shape =  pb.do.* pb.gray;
    pb.col =  pb.so.* pb.do;
    
    
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

meanreturn
