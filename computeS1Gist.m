function g = computeS1Gist(stim, filters, fSiz, c1ScaleSS, numPhases,INCLUDEBORDERS)

% Modified from the matlab implementation of the C1 code of Thomas Serre.
%
%  Returns the C1 and the S1 units' activation given the
%  input image, stim.
%  filters, fSiz, c1ScaleSS, c1ScaleSS, c1OL, INCLUDEBORDERS are the
%  parameters of the C1 system
%
%  stim   - the input image must be grayscale (single channel) and
%   type ''double''
%
%%% For S1 unit computation %%%
%
% filters -  Matrix of Gabor filters of size max_fSiz x num_filters,
% where max_fSiz is the length of the largest filter and num_filters the
% total number of filters. Column j of filters matrix contains a n_jxn_j
% filter (reshaped as a column vector and padded with zeros).
%
% fSiz - Vector of size num_filters containing the various filter
% sizes. fSiz(j) = n_j if filters j is n_j x n_j (see variable filters
% above).
%
%%% For C1 unit computation %%%
%
% c1ScaleSS  - Vector defining the scale bands, i.e. a group of filter
% sizes over which a local max is taken to get the C1 unit responses,
% e.g. c1ScaleSS = [1 k num_filters+1] means 2 scale bands, the first
% one contains filters(:,1:k-1) and the second one contains
% filters(:,k:num_filters). If N pooling bands, c1ScaleSS should be of
% length N+1.
%
% c1SpaceSS - Vector defining the spatial pooling range for each scale
% band, i.e. c1SpaceSS(i) = m_i means that each C1 unit response in band
% i is obtained by taking a max over a local neighborhood of m_ixm_i S1
% units. If N bands then c1SpaceSS should be of size N.
%
% c1OL - Scalar value defining the overlap between C1 units. In scale
% band i, the C1 unit responses are computed every c1Space(i)/c1OL.
%
% INCLUDEBORDERS - the type of treatment for the image borders.

USECONV2 = 1; %should be faster if 1

USE_NORMXCORR_INSTEAD = 0;
if(nargin < 7)
    INCLUDEBORDERS = 0;
end
numScaleBands=length(c1ScaleSS)-1;  % convention: last element in c1ScaleSS is max index + 1
numScales=c1ScaleSS(end)-1;
%   last index in scaleSS contains scale index where next band would start, i.e., 1 after highest scale!!
numSimpleFilters=floor(length(fSiz)/numScales);
for iBand = 1:numScaleBands
    ScalesInThisBand{iBand} = c1ScaleSS(iBand):(c1ScaleSS(iBand+1) -1);
end

% Rebuild all filters (of all sizes)
%%%%%%%%
nFilts = length(fSiz); 
for iPhase = 1:numPhases
    for i = 1:nFilts
        sqfilter{i}{iPhase} = reshape(filters(1:(fSiz(i)^2),i,iPhase),fSiz(i),fSiz(i));
        if USECONV2
            sqfilter{i}{iPhase} = sqfilter{i}{iPhase}(end:-1:1,end:-1:1); %flip in order to use conv2 instead of imfilter
        end
    end
end

% Calculate all filter responses (s1)
%%%%%%%%
sqim = stim.^2;
% iUFilterIndex = 0;
% precalculate the normalizations for the usable filter sizes
uFiltSizes = unique(fSiz); 
for i = 1:length(uFiltSizes) 
    s1Norm{uFiltSizes(i)} = (sumfilter(sqim,(uFiltSizes(i)-1)/2)).^0.5;
    %avoid divide by zero
    s1Norm{uFiltSizes(i)} = s1Norm{uFiltSizes(i)} + ~s1Norm{uFiltSizes(i)};
end



%% --------------------------------------------------------------
%                           Calculate S1 
% ------------------------------------------------------------------------
for iPhase = 1:numPhases
    iUFilterIndex = 0;
    
    for iBand = 1:numScaleBands 
        for iScale = 1:length(ScalesInThisBand{iBand})
            for iFilt = 1:numSimpleFilters 
                iUFilterIndex = iUFilterIndex+1;
                if ~USECONV2
                    s1(:,:,iBand,iScale,iFilt,iPhase) = abs(imfilter(stim,sqfilter{iUFilterIndex}{iPhase},'symmetric','same','corr'));
                    
                    if(~INCLUDEBORDERS)
                        s1(:,:,iBand,iScale,iFilt,iPhase) = removeborders(s1(:,:,iBand,iScale,iFilt,iPhase),fSiz(iUFilterIndex));
                    end
                    s1{iBand}{iScale}{iFilt} = im2double(s1{iBand}{iScale}{iFilt}) ./ s1Norm{fSiz(iUFilterIndex)};
                else %not 100% compatible but 20% faster at least
 
                    s1(:,:,iBand,iScale,iFilt,iPhase) = abs(conv2padded(stim,sqfilter{iUFilterIndex}{iPhase}));
                    
                    if(~INCLUDEBORDERS)
                        s1(:,:,iBand,iScale,iFilt,iPhase) = removeborders(s1(:,:,iBand,iScale,iFilt,iPhase),fSiz(iUFilterIndex));
                    end
                    s1(:,:,iBand,iScale,iFilt,iPhase) = im2double(s1(:,:,iBand,iScale,iFilt,iPhase)) ./ s1Norm{fSiz(iUFilterIndex)};
                end
                
            end
        end
    end
    
%     ss1 = ss1 + s1(:,:,:,:,:,iPhase) ./ numPhases;
end

s1_energy = sqrt(s1(:,:,:,:,:,1).^2 + s1(:,:,:,:,:,2).^2);


g = gistS1Block(s1_energy,numScaleBands,2,numSimpleFilters);






function sout = removeborders(sin,siz)
sin = unpadimage(sin, [(siz+1)/2,(siz+1)/2,(siz-1)/2,(siz-1)/2]);
sin = padarray(sin, [(siz+1)/2,(siz+1)/2],0,'pre');
sout = padarray(sin, [(siz-1)/2,(siz-1)/2],0,'post');


