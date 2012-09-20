% function mC2 = extractC2forcell(filter, filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,cPatches,cImages,numPatchSizes, nprocess, ext, show_progress);
% 	function mC2 = extractC2forcell(filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,cPatches,cImages,numPatchSizes)
function [mC2,mC2_gray] = extractC2forcell_dogray_mode(cPatches_gray,cPatches,cImages,numPatchSizes,numChannel,outDir,dataset,mode)
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
    rot =  0:22.5:157.5;
    c1ScaleSS = [1:2:8];
    RF_siz    = [7:6:39];
    c1SpaceSS = [8:4:20];
    minFS     = 7;
    maxFS     = 39;
    div = [4:-.05:3.2];        
    Div       = div(1:3:end); 
    case 1
    numPhases = 1;    
    rot = 0:22.5:157.5;
c1ScaleSS = [1:2:18];
RF_siz    = [7:2:39];
c1SpaceSS = [8:2:22];
minFS     = 7;
maxFS     = 39;
div = [4:-.05:3.2];
Div       = div; 
    case 2
    numPhases = 4;    
    rot = 0:22.5:157.5;
c1ScaleSS = [1:2:18];
RF_siz    = [7:2:39];
c1SpaceSS = [8:2:22];
minFS     = 7;
maxFS     = 39;
div = [4:-.05:3.2];
Div       = div;  
end

imageSize = 256;   

 fprintf(1,'Initializing gabor filters -- full set...');
%creates the gabor filters use to extract the S1 layer
[fSiz,filter, filters,c1OL,numSimpleFilters] = init_color_gabor_phases(rot, RF_siz, Div,numChannel,numPhases);
[fSiz_gray,filters_gray,c1OL_gray,numSimpleFilters_gray] =init_gabor_phase(rot, RF_siz, Div,numPhases);  
fprintf(1,'done\n');


	%all the patches are being flipped. This is becuase in matlab conv2 is much faster than filter2
	for i = 1:numPatchSizes,
		[siz,numpatch] = size(cPatches{i});
        [siz_gray,numpatch_gray] = size(cPatches_gray{i});
		siz = sqrt(siz/(numChannel/2*length(rot)));%8���ŵ���4������
        siz_gray = sqrt(siz_gray/length(rot));%���ŵ�ȡ���ֵ��keep4������
        
		for j = 1:numpatch,
			tmp = reshape(cPatches{i}(:,j),[siz,siz,numChannel/2,length(rot)]);
            tmp_gray = reshape(cPatches_gray{i}(:,j),[siz_gray,siz_gray,length(rot)]);
			tmp = tmp(end:-1:1,end:-1:1,:,:);
            tmp_gray = tmp_gray(end:-1:1,end:-1:1,:);
			cPatches{i}(:,j) = tmp(:);
            cPatches_gray{i}(:,j) = tmp_gray(:);
		end
	end

mC2 = [];
 mC1 = [];
 mC2_gray = [];
 mC1_gray = [];

nImages = length(cImages);
% 	for i = 1:length(cImages), %for every input image
for i = 1:nImages
outFile = fullfile(outDir,sprintf('c2do_%i_dataset_%i_images_%i_phase.mat', ...
            dataset, i,numPhases));
outFile_gray = fullfile(outDir,sprintf('c2gray_%i_dataset_%i_images_%i_phase.mat', ...
            dataset, i,numPhases));
        
if exist(outFile,'file') && exist(outFile_gray,'file')

mC2 = load(outFile);
mC2 = mC2.mC2;      

mC2_gray = load(outFile_gray);
mC2_gray = mC2_gray.mC2_gray; 

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
            
              if unused ~= 1;
                  stim_gray = rgb2gray(stim);
              else
                  stim_gray = stim;
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
        
            
			c1  = [];
			ic2 = []; %bug fix
            c1_gray  = [];
			ic2_gray = [];
            
%             cc1 = [];
            
			for j = 1:numPatchSizes, %for every unique patch size        
				if isempty(c1),  %compute C2
					[tmpC2,tmpC2_gray] = C2_dogray_mode(stim_gray,filters_gray,fSiz_gray,cPatches_gray{j},stim,filter,filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,cPatches{j},numChannel,numPhases);
				else
					[tmpC2,tmpC2_gray] = C2_dogray_mode(stim_gray,filters_gray,fSiz_gray,cPatches_gray{j},stim, filter, filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,cPatches{j},c1,numChannel,numPhases);
				end

			ic2 = [ic2;tmpC2];
            ic2_gray = [ic2_gray;tmpC2_gray];

            end

			mC2 = [mC2, ic2];
            mC2_gray = [mC2_gray, ic2_gray];  
            
save(fullfile(outDir,sprintf('./c2do_%i_dataset_%i_images_%i_phase.mat', ...
            dataset,i,numPhases)), 'mC2','-v7.3');
        
save(fullfile(outDir,sprintf('./c2gray_%i_dataset_%i_images_%i_phase.mat', ...
            dataset,i,numPhases)), 'mC2_gray','-v7.3');

end
end
fprintf('\n');
	%matlabpool close
% 	textprogressbar('Done\n');
