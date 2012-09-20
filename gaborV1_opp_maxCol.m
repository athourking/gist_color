% function g = gaborV1_opp_maxCol(stim, filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,INCLUDEBORDERS)
%function [c1,s1] = C1(stim, filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,INCLUDEBORDERS)
function g = gaborV1_opp_maxCol(stim, filter,filters, fSiz, c1SpaceSS,  c1ScaleSS, c1OL,numChannel)

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
        
        ss1{iBand}{iScale} = sqrt(s1{1}{iBand}{iScale}.^2 + s1{2}{iBand}{iScale}.^2);
        
         for iFilt = 1:numSimpleFilters %
            for jj = 1:numChannel/2
                 sss1{iBand}{iScale}(:,:,jj,iFilt) = max( ss1{iBand}{iScale}(:,:,jj,iFilt), ss1{iBand}{iScale}(:,:,jj+numChannel/2,iFilt));
            end
         end
        
      
    end
end


% Calculate local pooling (c1)
%%%%%%%%

%   (1) pool over scales within band
% for iBand = 1:numScaleBands %�ӵ�һ���߶ȴ�ʼѭ��
%     for iFilt = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
%          for channel=1:numChannel
%         c1{iBand}(:,:,iFilt) = zeros(size(sss1(:,:,iBand,1,iFilt)));
%         for iScale = 1:length(ScalesInThisBand{iBand}); %�ӵ�һ���߶ȿ�ʼѭ��
%             c1{iBand}(:,:,iFilt) = max(c1{iBand}(:,:,iFilt),ss1(:,:,iBand,iScale,iFilt)); %ȡ��ÿ������ϴ�Ĳ�ͬ�߶��µ�s1��Ӧ
%         end
%     end
%     end
% end
c1 = {};
for iBand = 1:numScaleBands
    for jj=1:numChannel/2
        for iFilt = 1:numSimpleFilters
            
            c1{iBand}(:,:, jj,iFilt) = zeros(size(sss1{iBand}{1}(:,:,jj,iFilt)));
            
            for iScale = 1:length(ScalesInThisBand{iBand});
                c1{iBand}(:,:,jj,iFilt) = max(c1{iBand}(:,:, jj,iFilt),sss1{iBand}{iScale}(:,:,jj,iFilt));
            end
            
        end
    end
end


    
    
%   (2) pool over local neighborhood
for iBand = 1:numScaleBands %��ÿ���߶ȴ�ʼѭ��
    poolRange = (c1SpaceSS(iBand)); %����������С10
    for jj=1:numChannel/2
        for iFilt = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
        c1{iBand}(:,:,jj,iFilt) = maxfilter(c1{iBand}(:,:,jj,iFilt),[0 0 poolRange-1 poolRange-1]); %��ÿ������ȡ���Ľϴ��s1������̬ѧ���ʹ���
        end
    end
end

% 
% %   (3) subsample
% for iBand = 1:numScaleBands %�ӵ�һ���߶ȴ�ʼѭ��
%     sSS=ceil(c1SpaceSS(iBand)/c1OL);
%     clear T;
%     for iFilt = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
%         T(:,:,iFilt) = c1{iBand}(1:sSS:end,1:sSS:end,iFilt); %�������õ�ÿ���������Ӧ
%         %     T(:,:,iFilt) = c1{iBand}(1:2:end,1:2:end,iFilt); %�������õ�ÿ���������Ӧ
%     end
%     c1{iBand} = T;
% end


g = gistGabor_V1_maxCol(c1,numScaleBands,numChannel/2,numSimpleFilters);







function sout = removeborders(sin,siz)
sin = unpadimage(sin, [(siz+1)/2,(siz+1)/2,(siz-1)/2,(siz-1)/2]);
sin = padarray(sin, [(siz+1)/2,(siz+1)/2],0,'pre');
sout = padarray(sin, [(siz-1)/2,(siz-1)/2],0,'post');


