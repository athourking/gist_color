% Demo, scene recognition
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
% orientationsPerScale = [8 8 8 8];
numberBlocks = 4;
fc_prefilt = 4;
NtrainingPerClass = 100;

Nclasses = length(categories);


% % 
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute global features
scenes = dir(fullfile(HOMEIMAGES, '*.jpg'));
scenes = {scenes(:).name};
Nscenes = length(scenes);
% 
% % Precompute filter transfert functions (only need to do this once, unless
% % image size is changes):

% % 
% Loop: Compute global features for all scenes
% F = zeros([Nscenes Nfeatures]);
% C = zeros([Nscenes 1]);
% 
% for n = 1:Nscenes
% % for n = 1:100
%     disp([n Nscenes]);
%     scenes{n}(1)= lower(scenes{n}(1)); %��ͼ��������ͼ��������ĸ��ΪСд
%     C(n) = strmatch(scenes{n}(1), categories);
% end
% Nfeatures = size(F,2);
% % for n = 1:Nscenes
% %     scenes{n}(1)= lower(scenes{n}(1)); %��ͼ��������ͼ��������ĸ��ΪСд
% %     C(n) = strmatch(scenes{n}(1), categories);
% % end

outDir = sprintf('../results/%i_train',NtrainingPerClass);
% if ~exist(outDir)
%     mkdir(outDir);
% end

load(fullfile(outDir,sprintf('C.mat')));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Split training/test (train = index training samples, test = index test)
train = [];
for c = 1:Nclasses
    j = find(C==c);%�ҳ���c�ೡ����ͼ����
    t = randperm(length(j));%��j��ͼ���������
    train = [train; j(t(1:NtrainingPerClass))];%ȡǰ100��ͼ����Ϊtraining set
end %�õ�800��ͼ���training set��ÿ��100��
test = setdiff(1:Nscenes, train)';%


save(fullfile(outDir,sprintf('train_%i_train_10_split.mat', ...
             NtrainingPerClass)), 'train','-v7.3');
save(fullfile(outDir,sprintf('test_%i_train_10_split.mat', ...
             NtrainingPerClass)), 'test','-v7.3');
         
         
         
         