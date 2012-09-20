function g = gaborS1_gaussopp_norm(stim,filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,numChannel)

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
% numPhases = 2;
% nFilts = length(fSiz);
% for i = 1:nFilts
%     sqfilter{i} = reshape(filters(1:(fSiz(i)^2),i),fSiz(i),fSiz(i));
%     if USECONV2
%         sqfilter{i} = sqfilter{i}(end:-1:1,end:-1:1); %flip in order to use conv2 instead of imfilter (%bug_fix 6/28/2007);
%     end
% end

% Calculate all filter responses (s1)
%%%%%%%%
sqim = stim.^2;
iUFilterIndex = 0;
% precalculate the normalizations for the usable filter sizes
uFiltSizes = unique(fSiz);
% for i = 1:length(uFiltSizes)
%     s1Norm{uFiltSizes(i)} = (sumfilter(sqim,(uFiltSizes(i)-1)/2)).^0.5;
%     %avoid divide by zero
%     s1Norm{uFiltSizes(i)} = s1Norm{uFiltSizes(i)} + ~s1Norm{uFiltSizes(i)};
% end
% ss1 = {};
% s1 = [];
%
% for iPhase = 1:numPhases
%     iUFilterIndex = 0;

for iBand = 1:numScaleBands
    for iScale = 1:length(ScalesInThisBand{iBand})
        
        sc = (iBand-1)*length(ScalesInThisBand{iBand}) + iScale;
        

            s1{iBand}{iScale} = double_opponents_gaussian_norm(stim,filters{sc},numChannel);%double opponent simple cell
            s1{iBand}{iScale} = im2double(s1{iBand}{iScale});
    end
        
%         ss1{iBand}{iScale} = sqrt(s1{1}{iBand}{iScale}.^2 + s1{2}{iBand}{iScale}.^2);
        
%          for iFilt = 1:numSimpleFilters %
%             for jj = 1:numChannel/2
%                  sss1{iBand}{iScale}(:,:,jj,iFilt) = max( ss1{iBand}{iScale}(:,:,jj,iFilt), ss1{iBand}{iScale}(:,:,jj+numChannel/2,iFilt));
%             end
%          end
end


% ss1 = sqrt(s1{1}{iBand}{iScale}(:,:,iFilt,channel) .^2 + s1{2}{iBand}{iScale} .^2);
% ss1 = sqrt(s1(:,:,:,:,1,:,:) .^2 + s1(:,:,:,:,2,:,:) .^2);


% [xx,yy,m,s,f] = size(ss1);
g = gistGauss_S1_opp(s1,numScaleBands,2,numChannel);


% function sout = removeborders(sin,siz)
% sin = unpadimage(sin, [(siz+1)/2,(siz+1)/2,(siz-1)/2,(siz-1)/2]);
% sin = padarray(sin, [(siz+1)/2,(siz+1)/2],0,'pre');
% sout = padarray(sin, [(siz-1)/2,(siz-1)/2],0,'post');
% return
