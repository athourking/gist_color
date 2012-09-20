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
% % Parameters
addpath('/users/jz7/data/recent/ColorGist/gist/svm571');
addpath('/users/jz7/data/recent/ColorGist/gist/gaussian_opponent_model_for_colorTuning');
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

Nclasses = length(categories);


numChannel = 8;
numPhases = 2;
% rot = [90 -45 0 45];
rot =  0:22.5:22.5*7;
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
scenes = dir(fullfile(HOMEIMAGES, '*.jpg'));
scenes = {scenes(:).name};
Nscenes = length(scenes);
% 
% % Precompute filter transfert functions (only need to do this once, unless
% % image size is changes):
% G = createGabor(orientationsPerScale, imageSize);
% [fSiz,filters,c1OL,numSimpleFilters] = init_gabor_phase(rot, RF_siz, Div,numPhases); %Thomas's gabor
[fSiz,filter, filters,c1OL,numSimpleFilters] = init_color_gabor(rot, RF_siz, Div,numChannel);

Nfeatures = length(rot)*length(RF_siz)*numberBlocks^2 * numChannel;
% % 
% Loop: Compute global features for all scenes
F = zeros([Nscenes Nfeatures]);
C = zeros([Nscenes 1]);

for n = 1:Nscenes
% for n = 1:100
    disp([n Nscenes]);
    img = imread(fullfile(HOMEIMAGES, scenes{n}));
%     img = mean(img,3);%����ɫͼ���Ϊ�Ҷ�ͼ�񣬵���ֵ�ձ����255��������ʾ�󲿷�Ϊ��ɫ
    if size(img,1) ~= imageSize
        img = imresize(img, [imageSize imageSize], 'bilinear');
%         img = imresize(img,[imageSize NaN]);
    end
    
%     output = prefilt(double(img), fc_prefilt);
%     g = gistGabor(output, numberBlocks, G);
%     g = gaborS1_phase(output, filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,numPhases);
%     output = prefilt(double(img), fc_prefilt);
    if max(img(:)) > 1
        img = double(img) / 255;
    end
    output = 2 * img - 1;
    
    g = gaborS1_opp(output, filter,filters, fSiz, c1SpaceSS, ...
        c1ScaleSS, c1OL,numChannel);

    F(n,:) = g;
    scenes{n}(1)= lower(scenes{n}(1)); %��ͼ��������ͼ��������ĸ��ΪСд
    C(n) = strmatch(scenes{n}(1), categories);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split training/test (train = index training samples, test = index test)
train = [];
for c = 1:Nclasses
    j = find(C==c);%�ҳ���c�ೡ����ͼ����
    t = randperm(length(j));%��j��ͼ���������
    train = [train; j(t(1:NtrainingPerClass))];%ȡǰ100��ͼ����Ϊtraining set
end %�õ�800��ͼ���training set��ÿ��100��
test = setdiff(1:Nscenes, train)';%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Train SVM
% train one versus all and then take maximal score

% fig = figure;
for c = 1:Nclasses
    netc = svm(Nfeatures, 'rbf', 0.003, 100);
    netc = svmtrain(netc, F(train,:), 2*(C(train)==c)-1, [], 1);

    [Y, scores(c,:)] = svmfwd(netc, F(test,:)); %��c���Ӧ�Ĳ��Լ�Ԥ�����ͷ���

end
for k = 1:length(test)
    [foo, ctest_hat(k)] = max(scores(:,k));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot performance
% Confusion matrix:
Nclasses = length(categories);
Cm = zeros(Nclasses, Nclasses);
for j = 1:Nclasses
    for i = 1:Nclasses
        % row i, col j is the percentage of images from class i that
        % were missclassified as class j.
        Cm(i,j) = 100*sum((C(test)==i) .* (ctest_hat'==j))/sum(C(test)==i);
    end
end

perf = mean(diag(Cm));%NP

save S1_6scales_opp_8col_rect F
save perf_6scales_opp_8col_rect perf
% figure
% subplot(121)
% imagesc(Cm); axis('square'); colorbar
% subplot(122)
% bar(diag(Cm))%���Խ����ϵ�ֵ
% title(mean(diag(Cm)))%8������ʶ���ʾ�ֵ
% axis('square')
