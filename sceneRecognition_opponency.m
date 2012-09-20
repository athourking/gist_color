% Demo, scene recognition
%
% This script trains a SVM for classification.
%
% You can use any SVM code.
%
% Here, the code assumes you have installed the next toolbox:
% Support Vector Machine toolbox
% Version 2.51, January 2002
% available at: http://ida.first.fraunhofer.de/~anton/software.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath('E:\Projects@Brown\CNS\Y\hjpkg_v3_backup\hjpkg_v3');
addpath('/users/jz7/data/recent/ColorGist/gist/gaussian_opponent_model_for_colorTuning');
% % addpath('C:\Users\pencil\Documents\My Dropbox\recent\double_model_energy');

% Parameters
HOMEIMAGES = '/users/jz7/data/recent/ColorGist/images/spatial_envelope_256x256_static_8outdoorcategories';
HOMEANNOTATIONS = '/users/jz7/data/recent/ColorGist/annotations/spatial_envelope_256x256_static_8outdoorcategories';
% DB2 = LMdatabase(HOMEANNOTATIONS);
categories = {'tallbuilding','insidecity','street','highway','coast','opencountry','mountain','forest'};
% HOMEIMAGES = 'E:\paper and code\Modeling the shape of the scene a holistic representation of the spatial envelope\images2'
% categories = {'aisle','boston','car','face','kitchen','meeting','office','static_indoor'};
imageSize = 256;
% imageSize = 140;
orientationsPerScale = [8 8 8 8];
numberBlocks = 4;
fc_prefilt = 4;
NtrainingPerClass = 1;

Nclasses = length(categories);


READPATCHESFROMFILE = 1; %use patches that were already computed
%(e.g., from natural images)
patchSizes = [4 8 12 16]; %other sizes might be better, maybe not
%all sizes are required
numChannel = 8;

numPatchSizes = length(patchSizes);
numPatchesPerSize = 250; %more will give better results, but will
numFeatures = numPatchesPerSize * numPatchSizes;

%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute global features
scenes = dir(fullfile(HOMEIMAGES, '*.jpg'));
scenes = {scenes(:).name};
Nscenes = length(scenes);
%
% % Precompute filter transfert functions (only need to do this once, unless
% % image size is changes):

% %
% Loop: Compute global features for all scenes
C = zeros([Nscenes 1]);


c1Res = extractRandC1Patches(HOMEIMAGES, numPatchSizes, ...
    numPatchesPerSize, patchSizes,numChannel); %fix: extracting from positive only

save c1Res

exit

% Nfeatures = size(c1Res,1);
% 
% F = c1Res';
% clear c1Res;
% 
% 
% for n = 1:Nscenes
%     scenes{n}(1)= lower(scenes{n}(1)); %��ͼ��������ͼ��������ĸ��ΪСд
%     C(n) = strmatch(scenes{n}(1), categories);
% end
% 
% % load('E:\Projects@Brown\CNS\Y\hjpkg_v3_backup\hjpkg_v3\ColorGist\gist\0609\allgist.mat');
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Split training/test (train = index training samples, test = index test)
% train = [];
% for c = 1:Nclasses
%     j = find(C==c);%�ҳ���c�ೡ����ͼ����
%     t = randperm(length(j));%��j��ͼ���������
%     train = [train; j(t(1:NtrainingPerClass))];%ȡǰ100��ͼ����Ϊtraining set
% end %�õ�800��ͼ���training set��ÿ��100��
% test = setdiff(1:Nscenes, train)';%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Train SVM
% % train one versus all and then take maximal score
% 
% fig = figure;
% for c = 1:Nclasses
%     netc = svm(Nfeatures, 'rbf', 0.003, 100);
%     netc = svmtrain(netc, F(train,:), 2*(C(train)==c)-1, [], 1);
%     
%     [Y, scores(c,:)] = svmfwd(netc, F(test,:)); %��c���Ӧ�Ĳ��Լ�Ԥ�����ͷ���
%     
%     %     figure(fig)
%     %     clf
%     %     plot(Y==1,'.')
%     %     hold on
%     %     plot(C(test)==c);%��ʵ���������
%     %     title(100*sum((C(test)==c).*(Y==1))/sum(C(test)==c));%��ȷ�����
%     %     xlabel(c)
%     %     drawnow
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
% 
% save allSceneOpponencyC1
% 
% exit

% figure
% subplot(121)
% imagesc(Cm); axis('square'); colorbar
% subplot(122)
% bar(diag(Cm))%���Խ����ϵ�ֵ
% title(mean(diag(Cm)))%8������ʶ���ʾ�ֵ
% axis('square')
