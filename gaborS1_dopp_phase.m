function g = gaborS1_dopp_phase(stim, filter,filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,numChannel)

USECONV2 = 1; %should be faster if 1

USE_NORMXCORR_INSTEAD = 0;
if(nargin < 9)
    INCLUDEBORDERS = 0;
end
numScaleBands=length(c1ScaleSS)-1;  % convention: last element in c1ScaleSS is max index + 1
numScales=c1ScaleSS(end)-1;
%   last index in scaleSS contains scale index where next band would start, i.e., 1 after highest scale!!
numSimpleFilters = floor(length(fSiz)/numScales/numChannel);

for iBand = 1:numScaleBands
    ScalesInThisBand{iBand} = c1ScaleSS(iBand):(c1ScaleSS(iBand+1) -1);
end

% Rebuild all filters (of all sizes)
%%%%%%%%
numPhases = 2;

% Calculate all filter responses (s1)
%%%%%%%%
sqim = stim.^2;
iUFilterIndex = 0;
% precalculate the normalizations for the usable filter sizes
uFiltSizes = unique(fSiz);


for iBand = 1:numScaleBands
    for iScale = 1:length(ScalesInThisBand{iBand})
        
        sc = (iBand-1)*length(ScalesInThisBand{iBand}) + iScale;
        
        for iPhase = 1:numPhases
            s1{iPhase}{iBand}{iScale} = double_opponents_S1(stim,filter{sc}{1},filters{sc}{iPhase},numChannel,numSimpleFilters);%double opponent simple cell

            s1{iPhase}{iBand}{iScale} = im2double(s1{iPhase}{iBand}{iScale});
            
        end
        
%         ss1{iBand}{iScale} = sqrt(s1{1}{iBand}{iScale}.^2 + s1{2}{iBand}{iScale}.^2);
        ss1{iBand}{iScale} = (s1{1}{iBand}{iScale} + s1{2}{iBand}{iScale}) / 2;   
       
         for iFilt = 1:numSimpleFilters %
            for jj = 1:numChannel/2
                 sss1{iBand}{iScale}(:,:,jj,iFilt) = max( ss1{iBand}{iScale}(:,:,jj,iFilt), ss1{iBand}{iScale}(:,:,jj+numChannel/2,iFilt));
            end
         end
        
      
    end
end


% ss1 = sqrt(s1{1}{iBand}{iScale}(:,:,iFilt,channel) .^2 + s1{2}{iBand}{iScale} .^2);
% ss1 = sqrt(s1(:,:,:,:,1,:,:) .^2 + s1(:,:,:,:,2,:,:) .^2);


% [xx,yy,m,s,f] = size(ss1);
g = gistGabor_S1_opp(sss1,numScaleBands,2,numChannel/2,numSimpleFilters);

return
