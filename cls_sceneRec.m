%
% This script trains a libSVM for classification:
% available at: http://www.csie.ntu.edu.tw/~cjlin/libsvm/

addpath(genpath('/gpfs/data/jz7/libsvm/libsvm-3.1/libsvm-3.1'));


% Parameters
HOMEIMAGES = '/gpfs/data/tserre/jzhang/ColorGist/images/spatial_envelope_256x256_static_8outdoorcategories';
categories = {'tallbuilding','insidecity','street','highway','coast','opencountry','mountain','forest'};
imageSize = 256; 
orientationsPerScale = [8 8 8 8];
numberBlocks = 4;
fc_prefilt = 4;
NtrainingPerClass = 100;

numChannel = 8;
numPhases = 2;
rot_s1 = [0 90];
rot =  0:22.5:157.5;
c1ScaleSS = 1:2:8;
RF_siz    = 7:6:39;
c1SpaceSS = 8:4:20;
minFS     = 7;
maxFS     = 39;
div = 4:-.05:3.2;
Div       = div(1:3:end);


Nscenes = 2680;
Nfeatures_s1 = length(rot)*length(RF_siz)*numberBlocks^2;
Nfeatures_c1 = length(rot)*(length(c1ScaleSS)-1)*numberBlocks^2;

tmpDir = sprintf('../results/0920');
C = load(fullfile(tmpDir,sprintf('C.mat')));
C = C.C;


outDir = sprintf('../results/0920');
Fg = load(fullfile(outDir,'F'));         
Fg = Fg.F;

Fso = load(fullfile(outDir,'F'));         
Fso = Fso.F;

Fdo = load(fullfile(outDir,'F'));         
Fdo = Fdo.F;


Nclasses = 8;
Cm = zeros(Nclasses, Nclasses,3);
perf = zeros(3,1);

s = [1 3 5 7 8];

for split = 1:10
    tmpDir = sprintf('../results/%i_train/%i_split',NtrainingPerClass,split);
    train_label = load(fullfile(tmpDir,sprintf('train_%i_train_%i_split.mat', ...
             NtrainingPerClass,split)));
    train_label = train_label.train;
    test_label = load(fullfile(tmpDir,sprintf('test_%i_train_%i_split.mat', ...
             NtrainingPerClass,split)));         
    test_label = test_label.test;
    
    yTrain = C(train_label); yTest = C(test_label);
    
    XTrain.do = Fdo(train_label,:); XTest.do = Fdo(test_label,:);
    XTrain.so = Fso(train_label,:); XTest.so = Fso(test_label,:);
    XTrain.gray = Fg(train_label,:); XTest.gray = Fg(test_label,:);
    XTrain.rgb = Frgb(train_label,:); XTest.rgb = Frgb(test_label,:);    
    
    XTrain.so = zscore(XTrain.so,0,1); XTest.so = zscore(XTest.so,0,1);
    XTrain.do = zscore(XTrain.do,0,1); XTest.do = zscore(XTest.do,0,1);
    XTrain.gray = zscore(XTrain.gray,0,1); XTest.gray = zscore(XTest.gray,0,1);
    XTrain.rgb = zscore(XTrain.rgb,0,1);  XTest.rgb = zscore(XTest.rgb,0,1);
    
    XTrain.all = [XTrain.so,XTrain.do,XTrain.gray]; XTest.all = [XTest.so,XTrain.do,XTrain.gray]; 
	XTrain.col = [XTrain.so,XTrain.do]; XTest.col = [XTest.so,XTrain.do];
    XTrain.shape = [XTrain.so,XTrain.do,XTrain.gray]; XTest.shape = [XTest.so,XTrain.do,XTrain.gray]; 
    
    
    %% libsvm
    model.s1 = svmtrain_libsvm(yTrain, XTrain.s1, '-t 2 -c 1000 -b 1');
    model.c1 = svmtrain_libsvm(yTrain, XTrain.c1, '-t 2 -c 1000 -b 1');
    model.gray = svmtrain_libsvm(yTrain, XTrain.gray, '-t 2 -c 1000 -b 1');
    model.rgb = svmtrain_libsvm(yTrain, XTrain.rgb, '-t 2 -c 1000 -b 1');
    model.all = svmtrain_libsvm(yTrain, XTrain.all, '-t 2 -c 1000 -b 1');
    model.shape = svmtrain_libsvm(yTrain, XTrain.shape, '-t 2 -c 1000 -b 1');
    model.col = svmtrain_libsvm(yTrain, XTrain.col, '-t 2 -c 1000 -b 1');
    

    [~, accuracy.s1(:,split), pb.s1] = svmpredict(yTest, XTest.s1, model.s1,'-b 1');    
    [~, accuracy.c1(:,split), pb.c1] = svmpredict(yTest, XTest.c1, model.c1,'-b 1');    
    [~, accuracy.gray(:,split), pb.gray] = svmpredict(yTest, XTest.gray, model.gray,'-b 1');    
    [~, accuracy.rgb(:,split), pb.rgb] = svmpredict(yTest, XTest.rgb, model.rgb,'-b 1');    
    [~, accuracy.shape(:,split), ~] = svmpredict(yTest, XTest.shape, model.shape,'-b 1');    
    [~, accuracy.all(:,split), ~] = svmpredict(yTest, XTest, model.all,'-b 1');    
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
            Cm.all(i,j,split) = 100*sum((yTest'==i) .* (predict_label.all'==j))/sum(yTest==i);
            Cm.shape(i,j,split) = 100*sum((yTest'==i) .* (predict_label.shape'==j))/sum(yTest==i);
             Cm.col(i,j,split) = 100*sum((yTest'==i) .* (predict_label.col'==j))/sum(yTest==i);
        end
    end
	
	
    perf.all(split) = mean(diag(Cm.all(:,:,split)));
    perf.shape(split) = mean(diag(Cm.shape(:,:,split)));
     perf.col(split) = mean(diag(Cm.col(:,:,split)));

end


%