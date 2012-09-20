% function mC2 = extractC2forcell(filter, filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,cPatches,cImages,numPatchSizes, nprocess, ext, show_progress);
% 	function mC2 = extractC2forcell(filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,cPatches,cImages,numPatchSizes)
function mC2_so = extractC2forcell_so_mode(cPatches,cImages,numPatchSizes,numChannel,numTrain,outDir,dataset,mode)
	%
	%this function is a wrapper of C2. For each image in the cell cImages, 
	%it extracts all the values of the C2 layer 
	%for all the prototypes in the cell cPatches.
	%The result mC2 is a matrix of size total_number_of_patches \times number_of_images where
	%total_number_of_patches is the sum over i = 1:numPatchSizes of length(cPatches{i})
	%and number_of_images is length(cImages)
	%The C1 parameters used are given as the variables filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL
	%for more detail regarding these parameters see the help entry for C1
	%
	%See also C1

	%% a bug was fixed on Jul 01 2005
switch mode
    case 0
    numPhases = 2;    
    rot =  [0 90];
    c1ScaleSS = [1:2:8];
    RF_siz    = [7:6:39];
    c1SpaceSS = [8:4:20];
    minFS     = 7;
    maxFS     = 39;
    div = [4:-.05:3.2];        
    Div       = div(1:3:end);
    case 1
    numPhases = 1;    
    rot = [0 90];
    c1ScaleSS = [1:2:18];
    RF_siz    = [7:2:39];
    c1SpaceSS = [8:2:22];
    minFS     = 7;
    maxFS     = 39;
    div = [4:-.05:3.2];
    Div       = div;
    case 2
    numPhases = 4;    
    rot = [0 90];  
    c1ScaleSS = [1:2:18];
    RF_siz    = [7:2:39];
    c1SpaceSS = [8:2:22];
    minFS     = 7;
    maxFS     = 39;
    div = [4:-.05:3.2];
    Div       = div;  
end

% outDir = sprintf('/gpfs/data/jz7/ColorGist/results/%i_train/C2/%i_phases',numTrain,numPhases);


fprintf(1,'Initializing gabor filters -- full set...');
%creates the gabor filters use to extract the S1 layer
[fSiz,filter, filters,c1OL,numSimpleFilters] = init_color_gabor_phases(rot, RF_siz, Div,numChannel,numPhases);
fprintf(1,'done\n');



	%all the patches are being flipped. This is becuase in matlab conv2 is much faster than filter2
	for i = 1:numPatchSizes,
        [siz_so,numpatch_so] = size(cPatches{i});
        siz_so = sqrt(siz_so/(numChannel*numPhases*length(rot)));
        
		for j = 1:numpatch_so,
            tmp_so = reshape(cPatches{i}(:,j),[siz_so,siz_so,numChannel,length(rot),numPhases]);
            tmp_so = tmp_so(end:-1:1,end:-1:1,:,:,:);
            cPatches{i}(:,j) = tmp_so(:);
		end
	end


 mC2_so = [];
 mC1_so = [];

imageSize = 256;

nImages = length(cImages);
% 	for i = 1:length(cImages), %for every input image
for i = 1:nImages
outFile = fullfile(outDir,sprintf('c2so_%i_dataset_%i_images_%i_phase.mat', ...
            dataset, i,numPhases));

if exist(outFile,'file') 

mC2_so = load(outFile);
mC2_so = mC2_so.mC2_so;        

else
    
 fprintf(1,'%d:',i);
  stim = cImages{i};
  [m,n,unused] = size(stim);
         if m ~= imageSize || n ~= imageSize 
        stim = imresize(stim, [imageSize imageSize], 'bilinear');
         end      
            if max(stim(:))>1
                stim = double(stim)./255;
            end

            stim = 2*stim -1;
	
        if unused == 1;
            im(:,:,1) = stim;
            im(:,:,2) = stim;
            im(:,:,3) = stim;
            
            clear stim;
        stim = im;
        clear im;
        
        end
       
             c1_so  = [];
			ic2_so = [];            

			for j = 1:numPatchSizes, %for every unique patch size        
				if isempty(c1_so),  %compute C2
					tmpC2_so = C2_so_mode(stim,filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,cPatches{j},numChannel,numPhases);
				else
					tmpC2_so = C2_so_mode(stim, filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,cPatches{j},c1,numChannel,numPhases);
                end               
                ic2_so = [ic2_so;tmpC2_so];
            end
            mC2_so = [mC2_so, ic2_so];
            

save(fullfile(outDir,sprintf('./c2so_%i_dataset_%i_images_%i_phase.mat', ...
            dataset,i,numPhases)), 'mC2_so','-v7.3');
end
end
fprintf('\n');
	%matlabpool close
% 	textprogressbar('Done\n');
