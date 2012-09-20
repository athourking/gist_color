function c1Res = extractRandC1Patches_255(HOMEIMAGES, numPatchSizes, numPatchesPerSize, patchSizes,numChannel)
%extracts random prototypes as part of the training of the C2 classification
%system.
%Note: we extract only from BAND 2. Extracting from all bands might help
%cPatches the returned prototypes
%cItrainingOnly the training images
%numPatchesPerSize is the number of sizes in which the prototypes come
%numPatchesPerSize is the number of prototypes extracted for each size
%patchSizes is the vector of the patche sizes

if nargin<2
    numPatchSizes = 4;
    numPatchesPerSize = 250;
    patchSizes = 4:4:16;
end

images = dir(fullfile(HOMEIMAGES, '*.jpg'));
images = {images(:).name};
nImages = length(images);


%----Settings for Training the random patches--------%
% % 	rot = [90 -45 0 45];
% 	c1ScaleSS = [1 3];
% 	RF_siz    = [11 13];
% 	c1SpaceSS = [10];
% % 	minFS     = 11;
% % 	maxFS     = 13;
% % 	div = [4:-.05:3.2];
% % 	Div       = div(3:4);
rot = [90 -45 0 45];
c1ScaleSS = [1:2:18];
RF_siz    = [7:2:39];
c1SpaceSS = [8:2:22];
minFS     = 7;
maxFS     = 39;
div = [4:-.05:3.2];
Div       = div;
%--- END Settings for Training the random patches--------%

fprintf(1,'Initializing gabor filters -- partial set...');
% 	[fSiz, filters,c1OL,numSimpleFilters] = init_color_gaussian_from_gabor(rot, RF_siz, Div);
[fSiz,filter, filters,c1OL,numSimpleFilters] = init_color_gabor(rot, RF_siz, Div,numChannel);
fprintf(1,'done\n');

cPatches = cell(numPatchSizes,1);
bsize = [0 0];

pind = zeros(numPatchSizes,1);


for j = 1:numPatchSizes
    
    cPatches{j} = zeros(patchSizes(j)^2*numChannel/2,numPatchesPerSize);
    %         cPatches{j} = zeros(patchSizes(j)^2,numPatchesPerSize);
end

imageSize = 256;
numScaleBands=length(c1ScaleSS)-1; 
c1Res = [];

% 	for i = 1:numPatchesPerSize,

% 		ii = floor(rand*nImages) + 1;
for ii = 1:nImages
%     for ii = 1:2
    img = double(imread(fullfile(HOMEIMAGES, images{ii})));
    %     img = mean(img,3);%����ɫͼ���Ϊ�Ҷ�ͼ�񣬵���ֵ�ձ����255��������ʾ�󲿷�Ϊ��ɫ
    if size(img,1) ~= imageSize
        img = imresize(img, [imageSize imageSize], 'bilinear');
        % %         img = imresize(img,[imageSize NaN]);
    end
    
    if max(img(:))>1
        stim = img./255;
    end
    img = 2*img -1;
    
    [c1source,s1source] = C1(img, filter,filters, fSiz, c1SpaceSS, ...
        c1ScaleSS, c1OL,numChannel);%double opponent simple cell and complex cell
    % 		b = c1source{1}; %new C1 interface;
    % 		bsize(1) = size(b,1);
    % 		bsize(2) = size(b,2);
    % 		for j = 1:numPatchSizes,
    % 			xy = floor(rand(1,2).*(bsize-patchSizes(j)))+1;
    % 			tmp = b(xy(1):xy(1)+patchSizes(j)-1,xy(2):xy(2)+patchSizes(j)-1,:,:);
    % 			pind(j) = pind(j) + 1; % counting how many patches per size
    % 			cPatches{j}(:,pind(j)) = tmp(:);
    % 		end
    c1_band = [];
    for iBand = 1:numScaleBands
        tmpC1 = reshape(c1source{iBand}, [size(c1source{iBand},1)*size(c1source{iBand},2)*numChannel/2*4,1]);
        c1_band = [c1_band ; tmpC1];
    end
  
    c1Res = [c1Res,c1_band];
    
%     c1Res{ii} = c1source;
    
end


fprintf('\n');
