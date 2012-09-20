% Demo, scene recognition
function F = jun_gist_coNormCol_phase_nonmaxS1_indoor(mode)
% This script trains a SVM for classification.
%
% You can use any SVM code. 
%
% Here, the code assumes you have installed the next toolbox:
% Support Vector Machine toolbox
% Version 2.51, January 2002
% available at: http://ida.first.fraunhofer.de/~anton/software.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Parameters
% addpath('/users/jz7/data/ColorGist/gist/svm571');
addpath('/users/jz7/data/ColorGist/gist/gaussian_opponent_model_for_colorTuning');
addpath('/users/jz7/data/Tuning/The Steerable Pyramid/matlabPyrTools/matlabPyrTools');

datapath = '/users/jz7/data/database/indoorCVPR09';
train = importdata('/users/jz7/data/database/indoorCVPR09/TrainImages.txt');
test  = importdata('/users/jz7/data/database/indoorCVPR09/TestImages.txt');
cI = cell(2,1);
cI{1} = train;
cI{2} = test;
clear train test

% imageSize = 140; 
% orientationsPerScale = [8 8 8 8];
numberBlocks = 4;
fc_prefilt = 4;
NtrainingPerClass = 80;


numChannel = 8;
numPhases = 2;
% rot = [90 -45 0 45];
rot =  0:22.5:157.5;
% c1ScaleSS = [1:2:18];
c1ScaleSS = [1:2:8];
% % c1ScaleSS = [1 1.4 2 3.5 4 7];
% RF_siz    = [7:2:39];
RF_siz    = [7:6:39];
% c1SpaceSS = [8:2:22];
c1SpaceSS = [8:4:20];
minFS     = 7;
maxFS     = 39;
div = [4:-.05:3.2];
Div       = div(1:3:end);

% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute global features

% 
% % Precompute filter transfert functions (only need to do this once, unless
% % image size is changes):
% switch mode
%     case 0
%         rot = 0:22.5:157.5;
%     case 1
%         rot = [0 90];
% end

[fSiz,filter, filters,c1OL,numSimpleFilters] = init_color_gabor(rot, RF_siz, Div,numChannel,numPhases);

Nfeatures = length(rot)*(length(c1ScaleSS)-1)*numberBlocks^2 * numChannel * numPhases;
% % 
% Loop: Compute global features for all scenes


outDir = sprintf('/gpfs/data/jz7/ColorGist/results/indoor/%i_train/C1',NtrainingPerClass);
if ~exist(outDir)
    mkdir(outDir);
end

height  =200;


switch mode
    case 0
        ii = 1;
    case 1
        ii = 2;
end

F = zeros([numel(cI{ii}) Nfeatures]);
for n = 1:numel(cI{ii})

    outFile = fullfile(outDir,sprintf('F_phase_nonmaxS1_%i_orient_%i_scene_%i_dataset.mat', ...
            length(rot),n,ii));
        
    if exist(outFile,'file'),
        F = load(outFile);
        F = F.F;  
        
    else
    
    img = imread(fullfile(datapath, cI{ii}{n}));
    [mm,nn,unused] = size(img);
                ratio = height/min(mm,nn);
        img = imresize(img,ratio,'bilinear'); 
    if unused == 1
        im(:,:,1) = img;
        im(:,:,2) = img;
        im(:,:,3) = img;
        clear img
        img = im;
        clear im;
    end
    
    output = prefilt(double(img), fc_prefilt);
    g = C1_soNorm_phase_nonmaxS1(output, filters, fSiz, c1SpaceSS, ...
        c1ScaleSS, c1OL,numChannel,numPhases);

    F(n,:) = g;
    
    save(fullfile(outDir,sprintf('F_phase_nonmaxS1_%i_orient_%i_scene_%i_dataset.mat', ...
           length(rot), n,ii)) ,'F','-v7.3');
    end
end

    save(fullfile(outDir,sprintf('F_c1NormCol_phase_nonmaxS1_%i_orient_%i_dataset.mat', ...
           length(rot), ii)) ,'F','-v7.3');


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tmpDir = sprintf('../results/%i_train',NtrainingPerClass);
% train = load(fullfile(tmpDir,sprintf('train_%i_train.mat', ...
%              NtrainingPerClass)), 'train','-v7.3');
% train = train.train;         
% test = load(fullfile(tmpDir,sprintf('test_%i_train.mat', ...
%              NtrainingPerClass)), 'test','-v7.3');
% test = test.test;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Train SVM
% % train one versus all and then take maximal score
% 
% % fig = figure;
% for c = 1:Nclasses
%     netc = svm(Nfeatures, 'rbf', 0.003, 100);
%     netc = svmtrain(netc, F(train,:), 2*(C(train)==c)-1, [], 1);
% 
%     [Y, scores(c,:)] = svmfwd(netc, F(test,:)); %��c���Ӧ�Ĳ��Լ�Ԥ�����ͷ���
% 
% end
% for k = 1:length(test)
%     [foo, ctest_hat(k)] = max(scores(:,k));
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot performance
% % Confusion matrix:
% Nclasses = length(categories);
% Cm = zeros(Nclasses, Nclasses);
% for j = 1:Nclasses
%     for i = 1:Nclasses
%         % row i, col j is the percentage of images from class i that
%         % were missclassified as class j.
%         Cm(i,j) = 100*sum((C(test)==i) .* (ctest_hat'==j))/sum(C(test)==i);
%     end
% end
% 
% perf = mean(diag(Cm));%NP
% 
% save F_S1_6scales_opp_8col_4phase F
% save perf_S1_6scales_opp_8col_4phase perf


exit
% figure
% subplot(121)
% imagesc(Cm); axis('square'); colorbar
% subplot(122)
% bar(diag(Cm))%���Խ����ϵ�ֵ
% title(mean(diag(Cm)))%8������ʶ���ʾ�ֵ
% axis('square')
