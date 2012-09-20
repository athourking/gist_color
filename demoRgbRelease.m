% Demo, scene recognition
%
% This script illustrates the RGB-Gist feature extraction which is modified
% % from Antonio Torralba's MATLAB implementation while using Thomas Serre's
% Gabor in the space domain.

%% ---------------------------------------------------------------
%                                Parameters
% -------------------------------------------------------------------------
imgDir = 'C:\cognitive science\scene perception\gist\Modeling the shape of the scene a holistic representation of the spatial envelope\images\spatial_envelope_256x256_static_8outdoorcategories';
categories = {'tallbuilding','insidecity','street','highway','coast','opencountry','mountain','forest'};
imageSize = 256; 
numberBlocks = 4;
fc_prefilt = 4;
Nclasses = length(categories);


% params for Gabor filters
numPhases = 2;
rot =  0:22.5:22.5*7;
c1ScaleSS = 1:2:8;
RF_siz    = 7:6:39;
c1SpaceSS = 8:6:20;
div = 4:-.05:3.2;
Div       = div(1:3:end);




%% ---------------------------------------------------------------
%                       Compute RGB GIST features
%--------------------------------------------------------------------------
scenes = dir(fullfile(imgDir, '*.jpg'));
scenes = {scenes(:).name};
Nscenes = length(scenes);
% 
[fSiz,filters,c1OL,numSimpleFilters] = init_gabor_phase(rot, RF_siz, Div,numPhases);
Nfeatures = length(rot)*(length(c1ScaleSS)-1)*numberBlocks^2*3;



F = zeros([Nscenes Nfeatures]);
for n = 1:Nscenes
    disp([n Nscenes]);
    img = imread(fullfile(imgDir, scenes{n}));
    if size(img,1) ~= imageSize
        img = imresize(img, [imageSize imageSize], 'bilinear');
    end
    
    output = prefilt(double(img), fc_prefilt);
    
    g = zeros(Nfeatures/3,3);
    for kk = 1:size(img,3)
        g(:,kk) = computeGist(output(:,:,kk), filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,numPhases);
    end
    F(n,:) = reshape(g, [size(g,1)*3 size(g,2)/3]);

end


outDir = sprintf('../results');
if ~exist(outDir,'dir')
    mkdir(outDir);
end

save(fullfile(outDir,sprintf('F.mat')) ,'F','-v7.3');
