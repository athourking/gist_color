function cPatches = extractRandC1Patches_normPhase_mode(cItrainingOnly, numPatchSizes, numPatchesPerSize, patchSizes, numChannel,numTrain,mode)
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

     imageSize = 256; 
	nImages = length(cItrainingOnly);
%     nImages = size(cItrainingOnly,4);

	%----Settings for Training the random patches--------%
switch mode
    case 0
    numPhases = 2;
    rot =  0:22.5:157.5;
    c1ScaleSS = [1 3];
%     RF_siz    = [7:6:39];
     RF_siz    = [7 13];
%     c1SpaceSS = [8:4:20];
    c1SpaceSS = [8];
    minFS     = 7;
    maxFS     = 13;
    div = [4:-.05:3.2];        
    Div       = [div(1) div(4)];
    case 1
    numPhases = 1;    
    rot = 0:22.5:157.5;
	c1ScaleSS = [1 3];
	RF_siz    = [11 13];
	c1SpaceSS = [10];
	minFS     = 11;
	maxFS     = 13;
	div = [4:-.05:3.2];
	Div       = div(3:4);   
    case 2
    numPhases = 4;    
    rot = 0:22.5:157.5;
	c1ScaleSS = [1 3];
	RF_siz    = [11 13];
	c1SpaceSS = [10];
	minFS     = 11;
	maxFS     = 13;
	div = [4:-.05:3.2];
	Div       = div(3:4);   
end

imageSize = 256;
outDir = sprintf('/gpfs/data/jz7/ColorGist/results/%i_train/C2/%i_phases',numTrain,numPhases);
if ~exist(outDir)
    mkdir(outDir);
end
	%--- END Settings for Training the random patches--------%

	fprintf(1,'Initializing gabor filters -- partial set...');
	[fSiz,filter, filters,c1OL,numSimpleFilters] = init_color_gabor_phases(rot, RF_siz, Div,numChannel,numPhases);
	fprintf(1,'done\n');

	cPatches = cell(numPatchSizes,1);
	bsize = [0 0];

	pind = zeros(numPatchSizes,1);


	for j = 1:numPatchSizes

% 		cPatches{j} = zeros(patchSizes(j)^2*length(rot)*numChannel,numPatchesPerSize);%4������8����ɫ�ŵ�
        		cPatches{j} = zeros(patchSizes(j)^2*length(rot)*numChannel/2,numPatchesPerSize);%4������8����ɫ�ŵ�
% cPatches{j} = zeros(patchSizes(j)^2*4,numPatchesPerSize);%4������8����ɫ�ŵ�ȡ������?
	end

% 	textprogressbar('Generating dictionary: ');

	for i = 1:numPatchesPerSize
%         i = 20;
		ii = floor(rand*nImages) + 1;

		stim = cItrainingOnly{ii};
       [m,n,unused] = size(stim);
       if m ~= imageSize || n ~= imageSize 
        stim = imresize(stim, [imageSize imageSize], 'bilinear');
       end        
     
         if max(stim(:))>1
                stim = double(stim)./255;
         end
        stim = 2*stim - 1;
		img_siz = size(stim);
        
        if unused == 1;
            im = zeros(size(stim));
            im(:,:,1) = stim;
            im(:,:,2) = stim;
            im(:,:,3) = stim;
            
             clear stim;
        stim = im;
        clear im;
        end
       

		[c1source,d1source,s1source] = C1_dopp_maxCol_normPhase(stim, filter,filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,numChannel,numPhases);
		b = c1source{1}; %new C1 interface;
		bsize(1) = size(b,1);
		bsize(2) = size(b,2);
		for j = 1:numPatchSizes,
			xy = floor(rand(1,2).*(bsize-patchSizes(j)))+1;
			tmp = b(xy(1):xy(1)+patchSizes(j)-1,xy(2):xy(2)+patchSizes(j)-1,:,:,:);
			pind(j) = pind(j) + 1; % counting how many patches per size
			cPatches{j}(:,pind(j)) = tmp(:);
		end
	end

% 	textprogressbar('done')
    save(fullfile(outDir,sprintf('dict_c1do_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases, numPatchesPerSize, length(patchSizes))) ,'cPatches','-v7.3');

	fprintf('\n');
