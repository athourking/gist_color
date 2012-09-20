% Demo, scene recognition
function jun_sceneRecognition(mode)
% This script trains a SVM for classification.
%
% You can use any SVM code. 
%
% Here, the code assumes you have installed the next toolbox:
% Support Vector Machine toolbox
% Version 2.51, January 2002
% available at: http://ida.first.fraunhofer.de/~anton/software.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
addpath('/users/jz7/data/ColorGist/gist/svm571');
addpath('/users/jz7/data/ColorGist/gist/gaussian_opponent_model_for_colorTuning');
addpath('/users/jz7/data/Tuning/The Steerable Pyramid/matlabPyrTools/matlabPyrTools');

% Parameters
HOMEIMAGES = '/users/jz7/data/ColorGist/images/spatial_envelope_256x256_static_8outdoorcategories';
HOMEANNOTATIONS = '/users/jz7/data/ColorGist/annotations/spatial_envelope_256x256_static_8outdoorcategories';
% DB2 = LMdatabase(HOMEANNOTATIONS);
categories = {'tallbuilding','insidecity','street','highway','coast','opencountry','mountain','forest'};
% HOMEIMAGES = 'E:\paper and code\Modeling the shape of the scene a holistic representation of the spatial envelope\images2'
% categories = {'aisle','boston','car','face','kitchen','meeting','office','static_indoor'};
imageSize = 256; 
% imageSize = 140; 
orientationsPerScale = [8 8 8 8];
numberBlocks = 4;
fc_prefilt = 4;
NtrainingPerClass = 100;

Nclasses = length(categories);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute global features
scenes = dir(fullfile(HOMEIMAGES, '*.jpg'));
scenes = {scenes(:).name};
Nscenes = length(scenes);
% 
% % Precompute filter transfert functions (only need to do this once, unless
% % image size is changes):
% % createGabor(orientationsPerScale, imageSize);
G = createGabor(orientationsPerScale, imageSize);
Nfeatures = size(G,3)*numberBlocks^2;
% % 
% Loop: Compute global features for all scenes
F = zeros([Nscenes Nfeatures]);
C = zeros([Nscenes 1]);


switch mode
    case 0
        split = 2;
    case 1
        split = 3;
end

outDir =  sprintf('../results/%i_train/%i_split',NtrainingPerClass,split);
        

for n = 1:Nscenes
% for n = 1:2
    disp([n Nscenes]);
    img = imread(fullfile(HOMEIMAGES, scenes{n}));
    img = mean(img,3);
    if size(img,1) ~= imageSize
        img = imresize(img, [imageSize imageSize], 'bilinear');
%         img = imresize(img,[imageSize NaN]);
    end
    
    output = prefilt(img, fc_prefilt);
    g = gistGabor(output, numberBlocks, G);
    
    F(n,:) = g;
%     scenes{n}(1)= lower(scenes{n}(1)); %��ͼ��������ͼ��������ĸ��ΪСд
%     C(n) = strmatch(scenes{n}(1), categories);
end

save(fullfile(outDir,sprintf('F_gray.mat')) ,'F','-v7.3');

% for n = 1:Nscenes
%     scenes{n}(1)= lower(scenes{n}(1)); %��ͼ��������ͼ��������ĸ��ΪСд
%     C(n) = strmatch(scenes{n}(1), categories);
% end



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
% perf = mean(diag(Cm))%NP
% 
% 
% 
% save gist_torralba
% 
% 
% exit
% % 
% % 
% % figure
% % subplot(121)
% % imagesc(Cm); axis('square'); colorbar
% % subplot(122)
% % bar(diag(Cm))%���Խ����ϵ�ֵ
% % title(mean(diag(Cm)))%8������ʶ���ʾ�ֵ
% % axis('square')
