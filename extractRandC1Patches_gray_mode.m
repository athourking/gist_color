function cPatches = extractRandC1Patches_gray_mode(cItrainingOnly, numPatchSizes, numPatchesPerSize, patchSizes,numTrain,mode);
%extracts random prototypes as part of the training of the C2 classification 
%system. 
%Note: we extract only from BAND 2. Extracting from all bands might help
%cPatches the returned prototypes
%cItrainingOnly the training images
%numPatchesSize is the number of sizes in which the prototypes come
%numPatchesPerSize is the number of prototypes extracted for each size
%patchSizes is the vector of the patche sizes

if nargin<2
  numPatchSizes = 4;
%   numPatchSizes = 3;
  numPatchesPerSize = 250;
  patchSizes = 4:4:16;
%   patchSizes = 4:4:12;
end
outDir = sprintf('/gpfs/data/jz7/ColorGist/results/%i_train/C2',numTrain);

nImages = length(cItrainingOnly);
%%%%ע�⣺������ȡ����band 2�е�������С���˲���
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
%--- END Settings for Training the random patches--------%

fprintf(1,'Initializing gabor filters -- partial set...');%�������ʾInitializing gabor filters -- partial set...K
% [fSiz,filters,c1OL,numSimpleFilters] = init_gabor(rot, RF_siz, Div);
[fSiz,filters,c1OL,numSimpleFilters] = init_gabor_phase(rot, RF_siz, Div,numPhases); %Thomas's gabor

fprintf(1,'done\n');

% %%%%%%%%%%%%gabor��ά%%%%%%%%%%%%%
%  A = filters;
%        n =30;
%        [Vecs1,Vals1,Psi1] = pc_evectors(A,n);
%        filters = Vecs1;
       
cPatches = cell(numPatchSizes,1);
bsize = [0 0];

pind = zeros(numPatchSizes,1);
for j = 1:numPatchSizes
  cPatches{j} = zeros(patchSizes(j)^2*length(rot),numPatchesPerSize);
end

for i = 1:numPatchesPerSize,
  ii = floor(rand*nImages) + 1;
  fprintf(1,'.');
  stim = cItrainingOnly{ii}; %һ��ѵ��ͼ��
  [m,n,unused] = size(stim);
         if m ~= imageSize || n ~= imageSize 
        stim = imresize(stim, [imageSize imageSize], 'bilinear');
       end

         if max(stim(:))>1
                stim = double(stim)./255;
         end
  if unused ~= 1;
      stim = rgb2gray(stim);
  end
  
 
  [c1source,s1source] = C1_phase(stim, filters, fSiz, c1SpaceSS, ...
      c1ScaleSS, c1OL,numPhases);
  b = c1source{1}; %new C1 interface;
  bsize(1) = size(b,1);
  bsize(2) = size(b,2);
		for j = 1:numPatchSizes,
			xy = floor(rand(1,2).*(bsize-patchSizes(j)))+1;
			tmp = b(xy(1):xy(1)+patchSizes(j)-1,xy(2):xy(2)+patchSizes(j)-1,:);
			pind(j) = pind(j) + 1; % counting how many patches per size
			cPatches{j}(:,pind(j)) = tmp(:);
		end

end

save(fullfile(outDir,sprintf('dict_c1gray_%i_train_%i_phase_%i_patches_%i_sizes.mat', ...
             numTrain,numPhases, numPatchesPerSize, length(patchSizes))) ,'cPatches','-v7.3');


% fprintf('\n');
