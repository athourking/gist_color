% Demo, scene recognition
%
% This script illustrates Double-Opponent color Gist features
%



%% ---------------------------------------------------------------
%                                Parameters
% -------------------------------------------------------------------------
imgDir = '..\images\spatial_envelope_256x256_static_8outdoorcategories';
categories = {'tallbuilding','insidecity','street','highway','coast','opencountry','mountain','forest'};
imageSize = 256; 
numberBlocks = 4;
fc_prefilt = 4;
Nclasses = length(categories);

% params for Gabor filters
numPhases = 2;
numChannels = 8;
rot =  0:22.5:22.5*7;
c1ScaleSS = 1:2:8;
RF_siz    = 7:6:39;
c1SpaceSS = 8:6:20;
div = 4:-.05:3.2;
Div       = div(1:3:end);




%% ---------------------------------------------------------------
%                       Compute SO GIST features
%--------------------------------------------------------------------------
% Compute global features
scenes = dir(fullfile(imgDir, '*.jpg'));
scenes = {scenes(:).name};
Nscenes = length(scenes);
%

fprintf(1,'Initializing color gabor filters -- full set...');
%creates the gabor filters use to extract the S1 layer
[fSiz,gfilters,cfilters,c1OL,numOrients] = init_color_gabor(rot, RF_siz, Div,numChannels,numPhases);
fprintf(1,'done\n');

Nfeatures = length(rot)*length(RF_siz)*numberBlocks^2 * numChannels/2;



F = zeros([Nscenes Nfeatures]);
for n = 1:Nscenes
    disp([n Nscenes]);
    img = imread(fullfile(imgDir, scenes{n}));
    if size(img,1) ~= imageSize
        img = imresize(img, [imageSize imageSize], 'bilinear');
    end
    
    output = prefilt(double(img), fc_prefilt);
    F(n,:) = computeDoGist(output,gfilters, cfilters, fSiz, ...
        c1ScaleSS,numPhases,numChannels);
end

outDir = sprintf('../results');
if ~exist(outDir,'dir')
    mkdir(outDir);
end


save(fullfile(outDir,sprintf('F.mat')) ,'F','-v7.3');

%