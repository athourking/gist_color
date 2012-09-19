function g = computeDoGist(stim, gfilters, cfilters,  fSiz, c1ScaleSS, numPhases, numChannels)

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




%% --------------------------------------------------------------
%                Calculate all filter responses (s1)
% -------------------------------------------------------------------------
s1 = {};
for iBand = 1:numScaleBands
    for iScale = 1:length(ScalesInThisBand{iBand})
        sc = (iBand-1)*length(ScalesInThisBand{iBand}) + iScale;
        
         % -----------Compute Single-Opponency----------%
        for iPhase = 1:numPhases
             iUFilterIndex = 0;
             
             s1{iBand}{iScale}(:,:,:,:,iPhase) = computeSOhmax(stim,cfilters{sc}{iPhase},numChannels,numSimpleFilters);%SOS1 unit
          
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
      
        % --------Divisive normalization over orientations---------%
        % Note: this is different from normalization used for computing SO
        % descriptors
        s1{iBand}{iScale} = divNorm_do(s1{iBand}{iScale},k,sigma,numSimpleFilters);
        
        
        
        % -----------Compute Double-Opponency----------------%
        % ds: Double-Opponent simple cell(DOS1)
        % dc: Double-Opponnet complex cell(DOC1)
        % gfilters is used at the DO stage is the same as the one used at the SO
        % stage but in the general case any filter with excitatory and inhibitory 
        % components could be used.
        
        tmpdc = zeros(size(stim,1),size(stim,2),numChannels,numSimpleFilters);
        for iPhase = 1:numPhases
            % DOS1 unit
            ds = computeDOS1hmax(s1{iBand}{iScale}(:,:,:,:,iPhase),gfilters{sc}{1},numChannels,numSimpleFilters);
            % DOC1 unit
            tmpdc =  tmpdc + ds ./ numPhases;
        end

        % yield invariance to figure-ground reversal
        for jj = 1:numChannels/2
            dc{iBand}{iScale}(:,:,jj,:) = sqrt(tmpdc(:,:,jj,:).^2 + tmpdc(:,:,jj+numChannels/2,:).^2);
%             dc{iBand}{iScale}(:,:,jj,:) = max( tmpdc(:,:,jj,:), tmpdc(:,:,jj+numChannels/2,:));
        end
 
    end
end



g = gistDoBlock(dc,numScaleBands,length(ScalesInThisBand{1}),numChannels/2,numSimpleFilters);




function sout = removeborders(sin,siz)
sin = unpadimage(sin, [(siz+1)/2,(siz+1)/2,(siz-1)/2,(siz-1)/2]);
sin = padarray(sin, [(siz+1)/2,(siz+1)/2],0,'pre');
sout = padarray(sin, [(siz-1)/2,(siz-1)/2],0,'post');
return
