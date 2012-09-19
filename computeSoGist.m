function g = computeSoGist(stim, filters, fSiz, c1SpaceSS, c1ScaleSS, numPhases, numChannels)

% USECONV2 = 1; %should be faster if 1

k = 1;
sigma = 0.225;



if(nargin < 8)
    INCLUDEBORDERS = 0;
end
numScaleBands=length(c1ScaleSS)-1;  % convention: last element in c1ScaleSS is max index + 1
numScales=c1ScaleSS(end)-1;
%   last index in scaleSS contains scale index where next band would start, i.e., 1 after highest scale!!
numSimpleFilters = floor(length(fSiz)/numScales/numChannels);

ScalesInThisBand = cell(1,numScaleBands);
for iBand = 1:numScaleBands
    ScalesInThisBand{iBand} = c1ScaleSS(iBand):(c1ScaleSS(iBand+1) -1);
end



%% ---------------------------------------------------------------
%                Calculate all filter responses (s1)
% -------------------------------------------------------------------------
s1 = {};
for iBand = 1:numScaleBands
    for iScale = 1:length(ScalesInThisBand{iBand})
        sc = (iBand-1)*length(ScalesInThisBand{iBand}) + iScale;
        
        for iPhase = 1:numPhases
            iUFilterIndex = 0;
            
            % -----------Compute Single-Opponency-----------%
            s1{iBand}{iScale}(:,:,:,:,iPhase) = computeSOhmax(stim,filters{sc}{iPhase},numChannels,numSimpleFilters);%SOS1 unit
            
            if(~INCLUDEBORDERS)
                for jj=1:numChannels
                    for ii = 1:numSimpleFilters
                        iUFilterIndex = iUFilterIndex+1;
                        s1{iBand}{iScale}(:,:,jj,ii,iPhase) = removeborders(s1{iBand}{iScale}(:,:,jj,ii,iPhase),fSiz(iUFilterIndex+numChannels*numSimpleFilters*(sc-1)));
                    end
                end
            end
            s1{iBand}{iScale} = im2double(s1{iBand}{iScale});
        end
        
        
        % ------Divisive normalization over opponent color channels-----%
        s1{iBand}{iScale} = divNorm_so(s1{iBand}{iScale},k,sigma,numChannels);
        
    end
end



%% ---------------------------------------------------------------
%               Calculate local pooling (assuming:SOC1 units)
% -------------------------------------------------------------------------

%   (1) pool over scales within band
c1 = {};
for iBand = 1:numScaleBands
    for jj=1:numChannels
        for iFilt = 1:numSimpleFilters
            for iPhase = 1:numPhases
                
                c1{iBand}(:,:,jj,iFilt,iPhase) = zeros(size(s1{iBand}{1}(:,:,jj,iFilt,iPhase)));
                for iScale = 1:length(ScalesInThisBand{iBand});
                    c1{iBand}(:,:,jj,iFilt,iPhase) = max(c1{iBand}(:,:,jj,iFilt,iPhase),s1{iBand}{iScale}(:,:,jj,iFilt,iPhase));
                end
                
            end
        end
    end
end



%   (2) pool over local neighborhood
for iBand = 1:numScaleBands
    poolRange = (c1SpaceSS(iBand));
    for jj=1:numChannels
        for iFilt = 1:numSimpleFilters
            for iPhase = 1:numPhases
                c1{iBand}(:,:,jj,iFilt,iPhase) = maxfilter(c1{iBand}(:,:,jj,iFilt,iPhase),[0 0 poolRange-1 poolRange-1]);
            end
        end
    end
end

%
%   (3) subsample
for iBand = 1:numScaleBands
    sSS=ceil(c1SpaceSS(iBand)/c1OL);
    clear T;
    
    for jj=1:numChannels
        for iFilt = 1:numSimpleFilters
            for iPhase = 1:numPhases
                T(:,:,jj,iFilt,iPhase) = c1{iBand}(1:sSS:end,1:sSS:end,jj,iFilt,iPhase);
            end
        end
    end
    c1{iBand} = T;
    
end



g = gistSoBlock(c1,numScaleBands,numPhases,numChannels,numSimpleFilters);




function sout = removeborders(sin,siz)
sin = unpadimage(sin, [(siz+1)/2,(siz+1)/2,(siz-1)/2,(siz-1)/2]);
sin = padarray(sin, [(siz+1)/2,(siz+1)/2],0,'pre');
sout = padarray(sin, [(siz-1)/2,(siz-1)/2],0,'post');
return
