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
% Parameters
addpath(genpath('/gpfs/data/jz7/Saliency/kLab/gbvs'));
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
G = createGabor(orientationsPerScale, imageSize);
Nfeatures = size(G,3)*numberBlocks^2*3;
% % 
% Loop: Compute global features for all scenes
F = zeros([Nscenes Nfeatures]);
C = zeros([Nscenes 1]);

outDir = sprintf('/gpfs/scratch/jz7/colorGist/rgb/%i_train',NtrainingPerClass);
if ~exist(outDir)
    mkdir(outDir);
end

for n = 1:Nscenes
    disp([n Nscenes]);
    outFile = fullfile(outDir,sprintf('F_rgb_%i_scene.mat', ...
            n));
        
    if exist(outFile,'file'),
        F = load(outFile);
        F = F.F;  
        
    else
        
    img = imread(fullfile(HOMEIMAGES, scenes{n}));
%     img = mean(img,3);%����ɫͼ���Ϊ�Ҷ�ͼ�񣬵���ֵ�ձ����255��������ʾ�󲿷�Ϊ��ɫ
    if size(img,1) ~= imageSize
        img = imresize(img, [imageSize imageSize], 'bilinear');
% %         img = imresize(img,[imageSize NaN]);
    end

    output = prefilt(double(img), fc_prefilt);
    g = gistGabor(output, numberBlocks, G);
    
    F(n,:) = g;
    
    save(fullfile(outDir,sprintf('F_rgb_%i_scene.mat', ...
            n)) ,'F','-v7.3');
    end
end

exit
