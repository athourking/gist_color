% function [c1,c1_DO,c1_SO_energy,s1_SO] = C1_SO(stim, filter, filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,INCLUDEBORDERS) %energy model
%function [c1,s1] = C1(stim, filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,INCLUDEBORDERS)
function [c1,s1] = C1_heeger(stim, filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,INCLUDEBORDERS) %heeger

USECONV2 = 0; %should be faster if 1

USE_NORMXCORR_INSTEAD = 0;
if(nargin < 8)
    INCLUDEBORDERS = 0;
end
numScaleBands=length(c1ScaleSS)-1;  % convention: last element in c1ScaleSS is max index + 1
numScales=c1ScaleSS(end)-1;
%   last index in scaleSS contains scale index where next band would start, i.e., 1 after highest scale!!
numSimpleFilters=floor(length(fSiz)/numScales);
for iBand = 1:numScaleBands
    ScalesInThisBand{iBand} = c1ScaleSS(iBand):(c1ScaleSS(iBand+1) -1);
end

numPhases = 4; %4个相位

% Rebuild all filters (of all sizes)
%%%%%%%%
nFilts = length(fSiz);
% for i = 1:nFilts
%     sqfilter{i} = reshape(filters(1:(fSiz(i)^2),i),fSiz(i),fSiz(i));
%     if USECONV2
%         sqfilter{i} = sqfilter{i}(end:-1:1,end:-1:1); %flip in order to use conv2 instead of imfilter (%bug_fix 6/28/2007);
%     end
% end

% Calculate all filter responses (s1)
%%%%%%%%
sqim = stim.^2;

% precalculate the normalizations for the usable filter sizes
uFiltSizes = unique(fSiz);
% for i = 1:length(uFiltSizes)
%     s1Norm{uFiltSizes(i)} = (sumfilter(sqim,(uFiltSizes(i)-1)/2)).^0.5;
%     %avoid divide by zero
%     s1Norm{uFiltSizes(i)} = s1Norm{uFiltSizes(i)} + ~s1Norm{uFiltSizes(i)};
% end


%%
% single oppponency stage
for iPhase = 1:numPhases
    iUFilterIndex = 0;
    
    for iBand = 1:numScaleBands
        for iScale = 1:length(ScalesInThisBand{iBand})

            sc = (iBand-1)*length(ScalesInThisBand{iBand}) + iScale;
            s1_SO{iPhase}{iBand}{iScale} = single_opponents(stim, filters{sc}{iPhase});
        
%             for channel=1:8
%                 for iFilt = 1:4
%                     iUFilterIndex = iUFilterIndex+1;
%                     
%                     if(~INCLUDEBORDERS)
%                         s1_SO{iPhase}{iBand}{iScale}{channel}{iFilt} = removeborders(s1_SO{iPhase}{iBand}{iScale}{channel}{iFilt},fSiz(iUFilterIndex));
%                     end
%                     s1_SO{iPhase}{iBand}{iScale}{channel}{iFilt} = im2double(s1_SO{iPhase}{iBand}{iScale}{channel}{iFilt});
%                     
%                 end
%             end
            
        end
    end
end


%%
% sum and square single opponent cells out of 90 degree phases
 for iBand = 1:numScaleBands
        for iScale = 1:length(ScalesInThisBand{iBand})
            
            for channel=1:8
                for iFilt = 1:4

                    c1_SO_energy{iBand}{iScale}{channel}{iFilt} = sqrt(s1_SO{1}{iBand}{iScale}{channel}{iFilt}.^2 + s1_SO{2}{iBand}{iScale}{channel}{iFilt}.^2);
                end
            end
        end
 end
 clear  s1_SO
 
%%
% double oppponency stage   
iUFilterIndex = 0;
for iBand = 1:numScaleBands
        for iScale = 1:length(ScalesInThisBand{iBand})
            
            sc = (iBand-1)*length(ScalesInThisBand{iBand}) + iScale;
            c1_DO{iBand}{iScale} = double_opponents_seperate(c1_SO_energy{iBand}{iScale}, filter{sc}{1});
      
            
%             for channel=1:8
%                 for iFilt = 1:4
%                     iUFilterIndex = iUFilterIndex+1;
%             
%                     if(~INCLUDEBORDERS)
%                         c1_DO{iBand}{iScale}{channel}{iFilt} = removeborders(c1_DO{iBand}{iScale}{channel}{iFilt},fSiz(iUFilterIndex));
%                     end
%             %s1{iBand}{iScale}{iFilt} = im2double(s1{iBand}{iScale}{iFilt}) ./ s1Norm{fSiz(iUFilterIndex)};
%                     c1_DO{iBand}{iScale}{channel}{iFilt} = im2double(c1_DO{iBand}{iScale}{channel}{iFilt});
%                 end
%             end
        end
end  

clear c1_SO_energy 
% Calculate local pooling (c1)
%%%%%%%%

% sum and square s1 units out of 90 degree phases
%  for iBand = 1:numScaleBands
%         for iScale = 1:length(ScalesInThisBand{iBand})
%             
%             for channel=1:8
%                 for iFilt = 1:4
% 
%                     c1_energy{iBand}{iScale}{channel}{iFilt} = sqrt(s1{1}{iBand}{iScale}{channel}{iFilt}.^2 + s1{2}{iBand}{iScale}{channel}{iFilt}.^2);
%                 end
%             end
%         end
%  end
 
 
%   (1) pool over scales within band
for iBand = 1:numScaleBands
    for channel=1:8
        for iFilt = 1:4
            c1{iBand}(:,:, channel, iFilt) = zeros(size(c1_DO{iBand}{1}{channel}{iFilt}));
%             c1{iBand}(:,:, channel, iFilt) = zeros(size(c1_energy{iBand}{1}{channel}{iFilt}));

            for iScale = 1:length(ScalesInThisBand{iBand});
                c1{iBand}(:,:,channel, iFilt) = max(c1{iBand}(:,:, channel, iFilt),c1_DO{iBand}{iScale}{channel}{iFilt});%double opponency不需要rectification，因为max over scales时已经相当于将进行了rectification
            end
        end
    end
end


% % taking the max over orientations
% for iBand=1:numScaleBands
%     cc1{iBand}(:,:,:) = max(c1{iBand}, [], 4);
%     
% %     cc2{iBand} = max(cc1{iBand}(:,:,:) ,[],3);
% end
% c1 = cc1;
% % c1 = cc2;
% 
% % taking the max over channels
% for iBand=1:numScaleBands
%     cc1{iBand}(:,:,:) = max(c1{iBand}, [], 3);
% end
% c1 = cc1;




%   (2) pool over local neighborhood
for iBand = 1:numScaleBands
    poolRange = (c1SpaceSS(iBand));
    for channel =1:8
        for iFilt = 1:4
            c1{iBand}(:,:,channel, iFilt) = maxfilter(c1{iBand}(:,:,channel, iFilt),[0 0 poolRange-1 poolRange-1]);
            
%             c1{iBand}(:,:, channel) = maxfilter(c1{iBand}(:,:, channel),[0 0 poolRange-1 poolRange-1]);

%             c1{iBand}(:,:, iFilt) = maxfilter(c1{iBand}(:,:, iFilt),[0 0 poolRange-1 poolRange-1]);
            
%             c1{iBand} = maxfilter(c1{iBand},[0 0 poolRange-1 poolRange-1]);

        end
    end
end


% % %   (3) subsample
% % for iBand = 1:numScaleBands
% %     sSS=ceil(c1SpaceSS(iBand)/c1OL);
% %     clear T;
% %     for channel =1:8
% %         for iFilt = 1:4
% %             T(:,:, channel, iFilt) = c1{iBand}(1:sSS:end,1:sSS:end,channel,iFilt);
% %  
% %             T(:,:, channel) = c1{iBand}(1:sSS:end,1:sSS:end,channel);
% % 
% %             T(:,:, iFilt) = c1{iBand}(1:sSS:end,1:sSS:end,iFilt);
% %             
% %             T = c1{iBand}(1:sSS:end,1:sSS:end);
% % 
% %         end
% %     end
% %     c1{iBand} = T;
% % end



function sout = removeborders(sin,siz)

sin = unpadimage(sin, [(siz+1)/2,(siz+1)/2,(siz-1)/2,(siz-1)/2]);
sin = padarray(sin, [(siz+1)/2,(siz+1)/2],0,'pre');
sout = padarray(sin, [(siz-1)/2,(siz-1)/2],0,'post');

