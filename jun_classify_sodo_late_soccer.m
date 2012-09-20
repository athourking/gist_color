function [predict_label, accuracy, dec_values] = jun_classify_sodo_late_soccer(data_ind,datapath,numTrain,numPhase)

% addpath(genpath('/gpfs/data/jz7/Plot/haxby'));
% addpath(genpath('/gpfs/data/jz7/Plot/real2rgb'));
% addpath(genpath('/gpfs/data/jz7/tools'));
% addpath(genpath('/gpfs/data/jz7/Plot/barweb'));
addpath(genpath('/gpfs/data/jz7/libsvm/libsvm-3.1/libsvm-3.1'));
addpath(genpath('/gpfs/data/jz7/Object_recognition/standardmodelrelease/standardmodelmatlabrelease'));
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


yTrain = [];
yTest = [];
for ii = 1:7
     yTrain =  [yTrain,ii*ones(1,numTrain)];
     yTest =  [yTest,ii*ones(1,numTest)];
end


for split = 1:3

    outDir = sprintf('/gpfs/data/tserre/jzhang/object_recognition/results/soccer/%i_split',split);

    c2so = load(fullfile(outDir,sprintf('eccv_c2so_%i_orients_%i_split_%i_patches_%i_sizes.mat', ...
          2, split, numPatchesPerSize, length(patchSizes))));
    c2so = c2so.C2res;

    c2do = load(fullfile(outDir,sprintf('eccv_c2do_%i_orients_%i_split_%i_patches_%i_sizes.mat', ...
          4, split,numPatchesPerSize, length(patchSizes))));
    c2do = c2do.C2res;

    c2gray = load(fullfile(outDir,sprintf('eccv_c2_%i_orients_%i_split_%i_patches_%i_sizes.mat', ...
          4 ,split,numPatchesPerSize, length(patchSizes))));
    c2gray = c2gray.C2res;
    
%     tmpDir = sprintf('../results/101/hmax/soccer/%i_phase/%i_split',numPhases,split);
% 
%     tr_label = load(fullfile(tmpDir,sprintf('train_label_%i_split_cal101.mat', ...
%              split)));
%     tr_label = tr_label.train_label;
%     te_label = load(fullfile(tmpDir,sprintf('test_label_%i_split_cal101.mat', ...
%              split)));
%     te_label = te_label.test_label;

%     
   XTrain_gray = [c2gray{1}]; XTest_gray =[c2gray{2}];
  XTrain_do = [c2do{1}]; XTest_do =[c2do{2}];
    XTrain_so = [c2so{1}]; XTest_so =[c2so{2}];
   %  XTrain_shape = [c2do{1}; c2gray{1}]; XTest_shape =[c2do{2};c2gray{2}];

   % XTrain = [c2so{1};c2do{1}]; XTest =[c2so{2};c2do{2}];
    
    XTrain_s1 = zscore(XTrain_so',0,1);
    XTest_s1 = zscore(XTest_so',0,1);
   XTrain_c1 = zscore(XTrain_do',0,1);
    XTest_c1 = zscore(XTest_do',0,1);
   XTrain_gray = zscore(XTrain_gray',0,1);
   XTest_gray = zscore(XTest_gray',0,1);
%     XTrain_shape = zscore(XTrain_shape',0,1);
%     XTest_shape = zscore(XTest_shape',0,1);
%    XTrain = zscore(XTrain',0,1);
 %   XTest = zscore(XTest',0,1);
%       XTrain = XTrain_s1; XTest = XTest_s1;
%       XTrain = XTrain_c1; XTest = XTest_c1;
%      XTrain = XTrain_gray; XTest = XTest_gray;
%     
      XTrain_col = [XTrain_s1,XTrain_c1]; XTest_col =[XTest_s1,XTest_c1];
      XTrain = [XTrain_s1,XTrain_c1,XTrain_gray]; XTest =[XTest_s1,XTest_c1,XTest_gray];
      XTrain_shape = [XTrain_c1, XTrain_gray]; XTest_shape =[XTest_c1, XTest_gray];
 
    model_s1 = svmtrain_libsvm(yTrain', XTrain_s1, '-t 0 -c 1000 -b 1');%RBF kernel
    model_c1 = svmtrain_libsvm(yTrain', XTrain_c1, '-t 0 -c 1000 -b 1');%RBF kernel
    model_gray = svmtrain_libsvm(yTrain', XTrain_gray, '-t 0 -c 1000 -b 1');%RBF kerne
    model = svmtrain_libsvm(yTrain', XTrain, '-t 0 -c 1000 -b 1');%RBF kernel
    model_shape = svmtrain_libsvm(yTrain', XTrain_shape, '-t 0 -c 1000 -b 1');%RBF kernel
    model_col = svmtrain_libsvm(yTrain', XTrain_col, '-t 0 -c 1000 -b 1');%RBF kernel

    [predict_label_s1, accuracy_s1(:,split), dec_values_s1] = svmpredict(yTest', XTest_s1, model_s1,'-b 1');    
    [predict_label_c1, accuracy_c1(:,split), dec_values_c1] = svmpredict(yTest', XTest_c1, model_c1,'-b 1');    
    [predict_label_gray, accuracy_gray(:,split), dec_values_gray] = svmpredict(yTest', XTest_gray, model_gray,'-b 1');    
    [predict_label_sp, accuracy_sp(:,split), dec_values_sp] = svmpredict(yTest', XTest_shape, model_shape,'-b 1');    
    [predict_label, accuracy(:,split), dec_values] = svmpredict(yTest', XTest, model,'-b 1');    
    [predict_label_col, accuracy_col(:,split), dec_values_col] = svmpredict(yTest', XTest_col, model,'-b 1');    

%     for j = 1:Nclasses
%         for i = 1:Nclasses
%             Cm_s1(i,j,split) = 100*sum((yTest'==i) .* (predict_label_s1 ==j))/sum(yTest==i);
%             Cm_c1(i,j,split) = 100*sum((yTest'==i) .* (predict_label_c1 ==j))/sum(yTest==i);
%             Cm_gray(i,j,split) = 100*sum((yTest'==i) .* (predict_label_gray ==j))/sum(yTest==i);
%         end
%     end
%     perf_s1(split) = mean(diag(Cm_s1(:,:,split)));
%     perf_c1(split) = mean(diag(Cm_c1(:,:,split)));
%     perf_gray(split) = mean(diag(Cm_gray(:,:,split)));
    
    %% product
    dec_values =  dec_values_s1.* dec_values_c1.* dec_values_gray;
    dec_values_shape =  dec_values_c1.* dec_values_gray;
    dec_values_sodo =  dec_values_s1.* dec_values_c1;
% 
% %     dec_values =  dec_values_s1. dec_values_c1;
%     %% average
%     % dec_va    % dec_values =  (dec_values_s1 + dec_values_c1 + dec_values_gray)/3;
% 
 
%         [foo,predict_label(k)] = max(dec_values(k,:));
%       lues =  (dec_values_s1 + dec_values_c1 + dec_values_gray)/3;

    for k = 1:length(yTest)
        [foo,predict_label_all(k)] = max(dec_values(k,:));
        [foo_shape,predict_label_shape(k)] = max(dec_values_shape(k,:));
       [foo_sodo,predict_label_sodo(k)] = max(dec_values_sodo(k,:));
% 
    end
    
    for j = 1:Nclasses
        for i = 1:Nclasses
            Cm(i,j,split) = 100*sum((yTest'==i) .* (predict_label_all'==j))/sum(yTest'==i);
            Cm_shape(i,j,split) = 100*sum((yTest'==i) .* (predict_label_shape'==j))/sum(yTest'==i);
             Cm_sodo(i,j,split) = 100*sum((yTest'==i) .* (predict_label_sodo'==j))/sum(yTest'==i);
% 
        end
    end
    perf(split) = mean(diag(Cm(:,:,split)));%NP
    perf_shape(split) = mean(diag(Cm_shape(:,:,split)));%NP
     perf_sodo(split) = mean(diag(Cm_sodo(:,:,split)));%NP
%    

end

mean(accuracy(1,:))


%% errorbar for each category
for split = 1:3
    Cm_obj(split,:) = diag(Cm(:,:,split));
    Cm_obj_s1(split,:) =diag(Cm_s1(:,:,split));
    Cm_obj_c1(split,:) = diag(Cm_c1(:,:,split));
    Cm_obj_gray(split,:) = diag(Cm_gray(:,:,split));
end

for ii = 1:7
    my_barvalues(1,ii) = mean(Cm_obj_s1(:,ii));
    my_barvalues(2,ii) = mean(Cm_obj_c1(:,ii));
    my_barvalues(3,ii) = mean(Cm_obj_gray(:,ii));
    my_barvalues(4,ii) = mean(Cm_obj(:,ii));
    
    my_errors(1,ii) = std(Cm_obj_s1(:,ii));
    my_errors(2,ii) = std(Cm_obj_c1(:,ii));
    my_errors(3,ii) = std(Cm_obj_gray(:,ii));
    my_errors(4,ii) = std(Cm_obj(:,ii));    
end

figure;
% barweb(my_barvalues, my_errors, [], [], [], [], [], haxby, [], [], 2, 'axis');
% handles = barweb(my_barvalues,my_errors, 0.6, {'SO','DO','Gray','Full'}, 'Soccer', 'Model', 'Accuracy(%)', 'haxby', 'y', [], 2, 'plot');
% handles = barweb(my_barvalues, my_errors, [], [], [], [], [], bone, [], [], 1, 'axis')
barwitherr(my_errors, my_barvalues); 
colormap whed;
set(gca,'XTickLabel',{'SO','DO','Gray','Full'});
% errorbar_tick(h);
ylabel('Accuracy(%)');
xlabel('Models');
legend(catNames);
set(gca,'ygrid','on','GridLineStyle','--');
title('Soccer');


%% confusion matrix for each category
Cm = mean(Cm,3);
% [confus,numcorrect,precision,recall,F] = getcm (yTest,predict_label);
% figure;plotconfmat(Cm);%3d bar
% figure;plotconfmattext(Cm);
figure;prtScoreConfusionMatrix(predict_label',yTest,[ ],catNames);
colormap pink2;
title('Soccer');

% figure; plotConfusionMatrix(Cm);



% %% errorbar for each model
% % color
% my_barvalues(1,1) = mean(perf_s1);%so
% my_barvalues(1,2) = 75;%hue
% my_barvalues(1,3) = 75;%opponent
% my_barvalues(1,4) = 86;%color names
% my_barvalues(2,1) = mean(perf_c1);%do
% my_barvalues(2,2) = mean(perf_gray);%gray-hmax
% my_barvalues(2,3) = 58;%sift
% my_barvalues(2,4) = 0;
% my_barvalues(3,1) = mean(perf);%so
% my_barvalues(3,2) = 84;%hue
% my_barvalues(3,3) = 85;%opponent
% my_barvalues(3,4) = 89;%color names
% 
% % my_errors = [];
% my_errors(1,1) = std(perf_s1);
% my_errors(1,2) = 0;
% my_errors(1,3) = 0;
% my_errors(1,4) = 0;
% my_errors(2,1) = std(perf_c1);
% my_errors(2,2) = std(perf_gray);
% my_errors(2,3) = 0;
% my_errors(2,4) = 0;
% my_errors(3,1) = std(perf);
% my_errors(3,2) = 0;
% my_errors(3,3) = 0;
% my_errors(3,4) = 0;
% 
% 
% 
% figure;
% % handles = barweb(my_barvalues, my_errors, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend, error_sides, legend_type);
% handles = barweb(my_barvalues,my_errors, 0.6, {'Color','Shape','Color+shape'}, 'Soccer dataset', 'Model', 'Accuracy(%)', 'haxby', 'y', [], 2, 'plot');
% % barweb(my_barvalues, my_errors, [], 'Models', [], {'SO','DO','Gray','Full'}, [], summer, [], [], 2, 'axis');
% % barwitherr(my_errors, my_barvalues); 
% % set(gca,'XTickLabel',{'Color','Shape','Color+shape'});
% set(gca,'ygrid','on','GridLineStyle','--');
% % errorbar_tick(h);
% ylabel('Accuracy(%)');
% xlabel('Models');



%% overlapped errorbar
%% errorbar for each model
% bar1
clear my_barvalues my_errors
my_barvalues(1,1) = mean(perf);%so
my_barvalues(1,2) = 84;%hue
my_barvalues(1,3) = 85;%opponent
my_barvalues(1,4) = 89;%color names
my_barvalues(2,1) = mean(perf);%so
my_barvalues(2,2) = 84;%hue
my_barvalues(2,3) = 85;%opponent
my_barvalues(2,4) = 89;%color names

my_errors(1,1) = std(perf);
my_errors(1,2) = 0;
my_errors(1,3) = 0;
my_errors(1,4) = 0;
my_errors(2,1) = std(perf);
my_errors(2,2) = 0;
my_errors(2,3) = 0;
my_errors(2,4) = 0;

% bar2
my_barvalues2(1,1) = mean(perf_s1);%so
my_barvalues2(1,2) = 75;%hue
my_barvalues2(1,3) = 75;%opponent
my_barvalues2(1,4) = 86;%color names
my_barvalues2(2,1) = mean(perf_shape);%do+gray
my_barvalues2(2,2) = 58;%
my_barvalues2(2,3) = 58;%sift
my_barvalues2(2,4) = 58;

% my_errors = [];
my_errors2(1,1) = std(perf_s1);
my_errors2(1,2) = 0;
my_errors2(1,3) = 0;
my_errors2(1,4) = 0;
my_errors2(2,1) = std(perf_shape);
my_errors2(2,2) = 0;
my_errors2(2,3) = 0;
my_errors2(2,4) = 0;

figure;
k1 = 0.8;
bar1 = jun_barweb(my_barvalues,my_errors, 0.4, {'Color','Shape'}, 'Soccer', 'Model', 'Accuracy(%)', [], 'y', [], 2, 'plot','m');
set(bar1.bars,'BarWidth',k1); 

hold on;
bar2 = jun_barweb(my_barvalues2,my_errors2, 0.4, {'Color','Shape'}, [], [], [],[], 'y', [], 2, 'plot','c');
set(gca,'ygrid','on','GridLineStyle','--');
set(bar2.bars,'BarWidth',k1/2); 

hold off; 
legend('Color/Shape','Color+Shape')


% errorbar_tick(h);
ylabel('Accuracy(%)');
xlabel('Models');
bar1=bar(x, y1, 'FaceColor', 'b', 'EdgeColor', 'b'); 
set(bar1,'BarWidth',K); 



%% average for each category




return



