function  C2res_so = jun_recognition_hmax_8scenes(mode)

%demoRelease.m
%demonstrates how to use C2 standard model features in a pattern classification framework
% % clc
% % clear all


addpath(genpath('/gpfs/data/jz7/libsvm/libsvm-3.1/libsvm-3.1'));
% addpath(genpath('/gpfs/data/jz7/Object_recognition/VGG'));
addpath(genpath('/gpfs/data/jz7/Object_recognition/standardmodelrelease/standardmodelmatlabrelease'));
addpath('/gpfs/data/jz7/Object_recognition/caltech101');

HOMEIMAGES = '/users/jz7/data/ColorGist/images/spatial_envelope_256x256_static_8outdoorcategories';
HOMEANNOTATIONS = '/users/jz7/data/ColorGist/annotations/spatial_envelope_256x256_static_8outdoorcategories';
% DB2 = LMdatabase(HOMEANNOTATIONS);
categories = {'tallbuilding','insidecity','street','highway','coast','opencountry','mountain','forest'};
% HOMEIMAGES = 'E:\paper and code\Modeling the shape of the scene a holistic representation of the spatial envelope\images2'
% categories = {'aisle','boston','car','face','kitchen','meeting','office','static_indoor'};
imageSize = 256; 
% maxTest = 50;
numTrain = 100;

 
% mode = 2;
switch mode
    case 0
       C2res_so = jun_hmax_gist_8scenes(numTrain,mode);  
    case 1
        [C2res_do, C2res_gray] = jun_hmax_gist_8scenes_dogray(numTrain,mode);
    case 2
        C2res_so = jun_hmax_gist_8scenes(numTrain,mode);
                    
end

        
exit