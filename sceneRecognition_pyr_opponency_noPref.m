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
clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath('E:\Projects@Brown\CNS\Y\hjpkg_v3_backup\hjpkg_v3\ColorGist\Gist_NP\svm_v251\svm');
addpath('/users/jz7/data/recent/ColorGist/gist/gaussian_opponent_model_for_colorTuning');
% % addpath('C:\Users\pencil\Documents\My Dropbox\recent\double_model_energy');
addpath('/users/jz7/data/recent/Tuning/The Steerable Pyramid/matlabPyrTools/matlabPyrTools');


% Parameters
HOMEIMAGES = '/users/jz7/data/recent/ColorGist/images/spatial_envelope_256x256_static_8outdoorcategories';
HOMEANNOTATIONS = '/users/jz7/data/recent/ColorGist/annotations/spatial_envelope_256x256_static_8outdoorcategories';
% DB2 = LMdatabase(HOMEANNOTATIONS);
categories = {'tallbuilding','insidecity','street','highway','coast','opencountry','mountain','forest'};
% HOMEIMAGES = 'E:\paper and code\Modeling the shape of the scene a holistic representation of the spatial envelope\images2'
% categories = {'aisle','boston','car','face','kitchen','meeting','office','static_indoor'};
imageSize = 256;
% imageSize = 140;
% orientationsPerScale = [8 8 8 8];
numberBlocks = 4;
fc_prefilt = 4;
NtrainingPerClass = 100;
numChannel = 8;
Nclasses = length(categories);
w = 4;
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute global features
scenes = dir(fullfile(HOMEIMAGES, '*.jpg'));
scenes = {scenes(:).name};
Nscenes = length(scenes);
%
% % Precompute filter transfert functions (only need to do this once, unless
% % image size is changes):
% G = createGabor(orientationsPerScale, imageSize);
% [G, filters] = get_filters_gabor(orientationsPerScale, imageSize,numChannel);
% Nfeatures = size(G,3)*numberBlocks^2*numChannel;
Nfeatures = w*w*(2*6+2)-2*w*w;

% %
% Loop: Compute global features for all scenes
% F = zeros([Nscenes Nfeatures]);
C = zeros([Nscenes 1]);
numPrincomp = 80;


for n = 1:Nscenes
    % for n = 1:10
    disp([n Nscenes]);
    
    %% Load an image, and downsample to a size appropriate for the machine speed.
    img = imread(fullfile(HOMEIMAGES, scenes{n}));
    img = double(img);
    %     img = mean(img,3);%����ɫͼ���Ϊ�Ҷ�ͼ�񣬵���ֵ�ձ����255��������ʾ�󲿷�Ϊ��ɫ
    if size(img,1) ~= imageSize
        img = imresize(img, [imageSize imageSize], 'bilinear');
        % %         img = imresize(img,[imageSize NaN]);
    end
    
%     img = prefilt(img, fc_prefilt);
    tic; corrDn(img,[1 1; 1 1]/4,'reflect1',[2 2]); time = toc;
    imSubSample = min(max(floor(log2(time)/2+3),0),2);
    img = blurDn(img, imSubSample,'qmf9');
    
    g = gistPyr_opponency(img,imSubSample,w,numChannel);
    
    F(n,:) = g;
    
    scenes{n}(1)= lower(scenes{n}(1)); %��ͼ��������ͼ��������ĸ��ΪСд
    C(n) = strmatch(scenes{n}(1), categories);
    
end

save allPyr_opponency_rect_noPref


% load('E:\Projects@Brown\CNS\Y\hjpkg_v3_backup\hjpkg_v3\ColorGist\gist\pyr\allPyr_new.mat');
% 
% 
% %PCA
% [coeff,score] = princomp(F);
% g = F * coeff(:,1:numPrincomp);
% 
% [Vecs,Vals,Psi] = pc_evectors(F',numPrincomp);%Get top 3 PC evectors 
% g = F * Vecs(:,:); % Project onto evectors
% Nfeatures = numPrincomp;
% % % for n = 1:Nscenes
% % %     scenes{n}(1)= lower(scenes{n}(1)); %��ͼ��������ͼ��������ĸ��ΪСд
% % %     C(n) = strmatch(scenes{n}(1), categories);
% % % end
% 
% Nfeatures = size(F,2);
% 
% 
% load('E:\Projects@Brown\CNS\Y\hjpkg_v3_backup\hjpkg_v3\ColorGist\gist\0609\allgist_rgb.mat');
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
% end
% 
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
% figure
% subplot(121)
% imagesc(Cm); axis('square'); colorbar
% subplot(122)
% bar(diag(Cm))%���Խ����ϵ�ֵ
% title(mean(diag(Cm)))%8������ʶ���ʾ�ֵ
% axis('square')
