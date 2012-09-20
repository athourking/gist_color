function g = S1_soNorm_phase(stim, filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,numChannel,numPhases)

USECONV2 = 1; %should be faster if 1
k = 1;
sigma = 0.225;

USE_NORMXCORR_INSTEAD = 0;
if(nargin < 10)
    INCLUDEBORDERS = 0;
end
numScaleBands=length(c1ScaleSS)-1;  % convention: last element in c1ScaleSS is max index + 1
numScales=c1ScaleSS(end)-1;
%   last index in scaleSS contains scale index where next band would start, i.e., 1 after highest scale!!
numSimpleFilters = floor(length(fSiz)/numScales/numChannel);

for iBand = 1:numScaleBands
    ScalesInThisBand{iBand} = c1ScaleSS(iBand):(c1ScaleSS(iBand+1) -1);
end


% Calculate all filter responses (s1)
%%%%%%%%
sqim = stim.^2;
iUFilterIndex = 0;
% precalculate the normalizations for the usable filter sizes
uFiltSizes = unique(fSiz);

%
% for iPhase = 1:numPhases
%     iUFilterIndex = 0;
% iUFilterIndex = 0;
for iBand = 1:numScaleBands
    for iScale = 1:length(ScalesInThisBand{iBand})
        sc = (iBand-1)*length(ScalesInThisBand{iBand}) + iScale;
        ss1{iBand}{iScale} = zeros(size(stim, 1), size(stim,2),numChannel,numSimpleFilters);
        
        for iPhase = 1:numPhases
                    iUFilterIndex = 0;
            s1{iBand}{iScale}(:,:,:,:,iPhase) = single_opponents_S1(stim,filters{sc}{iPhase},numChannel,numSimpleFilters);%double opponent simple cell
             if(~INCLUDEBORDERS)
                 
                  for jj=1:numChannel
                      for ii = 1:numSimpleFilters
                          iUFilterIndex = iUFilterIndex+1;
                          s1{iBand}{iScale}(:,:,jj,ii,iPhase) = removeborders(s1{iBand}{iScale}(:,:,jj,ii,iPhase),fSiz(iUFilterIndex+numChannel*numSimpleFilters*(sc-1)));      
                      end
                  end
             end
            
            s1{iBand}{iScale} = im2double(s1{iBand}{iScale});
        end
     
         %normalization
        E_single = zeros(size(stim, 1), size(stim,2),numSimpleFilters);
        
   for ii = 1:numSimpleFilters
      for jj=1:numChannel
        for iPhase = 1:numPhases
          E_single(:,:,ii) = E_single(:,:,ii) + 0.25 *  s1{iBand}{iScale} (:,:,jj,ii,iPhase).^2;
        end     
      end
   end
% 
   
   for iPhase = 1:numPhases
for ii = 1:numSimpleFilters
    for jj=1:numChannel
           s1{iBand}{iScale}(:,:,jj,ii,iPhase) = sqrt(k * s1{iBand}{iScale}(:,:,jj,ii,iPhase).^2 ./ (sigma.^2 + E_single(:,:,ii)));
    end
end
   end
%    
%    
%    for iPhase = 1:numPhases
% % d{iPhase}{iBand{iScale} = double_opponents_d1(stim, s1{iPhase}{iBand}{iScale},filter{sc}{1},numChannel,numSimpleFilters);
%     ss1{iBand}{iScale} =  ss1{iBand}{iScale} +  s1{iBand}{iScale}(:,:,:,:,iPhase) / numPhases;
%    end
% %    

      
    end
end

g = gistGabor_S1_opp_phase(s1,numScaleBands,numel(s1{iBand}),numChannel,numSimpleFilters,numPhases);




% ss1 = sqrt(s1{1}{iBand}{iScale}(:,:,iFilt,channel) .^2 + s1{2}{iBand}{iScale} .^2);
% ss1 = sqrt(s1(:,:,:,:,1,:,:) .^2 + s1(:,:,:,:,2,:,:) .^2);

% c1 = {};
% for iBand = 1:numScaleBands
%         for jj=1:numChannel
%             for iFilt = 1:numSimpleFilters
%                 for iPhase = 1:numPhases
%                 c1{iBand}(:,:,jj,iFilt,iPhase) = zeros(size(s1{iBand}{1}(:,:,jj,iFilt,iPhase)));
%                 for iScale = 1:length(ScalesInThisBand{iBand});
%                     c1{iBand}(:,:,jj,iFilt,iPhase) = max(c1{iBand}(:,:,jj,iFilt,iPhase),s1{iBand}{iScale}(:,:,jj,iFilt,iPhase));
%                 end
%             end
%             
%         end
%     end
% end
% 
%     
%     
% %   (2) pool over local neighborhood
% for iBand = 1:numScaleBands %��ÿ���߶ȴ�ʼѭ��
%     poolRange = (c1SpaceSS(iBand)); %����������С10
%     for jj=1:numChannel
%         for iFilt = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
%             for iPhase = 1:numPhases
%                 c1{iBand}(:,:,jj,iFilt,iPhase) = maxfilter(c1{iBand}(:,:,jj,iFilt,iPhase),[0 0 poolRange-1 poolRange-1]); %��ÿ������ȡ���Ľϴ��s1������̬ѧ���ʹ���
%         end
%     end
% end
% end
% 
% % 
% %   (3) subsample
% for iBand = 1:numScaleBands %�ӵ�һ���߶ȴ�ʼѭ��
%     sSS=ceil(c1SpaceSS(iBand)/c1OL);
%     clear T;
%     
%     for jj=1:numChannel
%     for iFilt = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
%     for iPhase = 1:numPhases
%         T(:,:,jj,iFilt,iPhase) = c1{iBand}(1:sSS:end,1:sSS:end,jj,iFilt,iPhase); %�������õ�ÿ���������Ӧ
%         %     T(:,:,iFilt) = c1{iBand}(1:2:end,1:2:end,iFilt); %�������õ�ÿ���������Ӧ
%     end
%     end
%     end
%     c1{iBand} = T;
% 
% end

function sout = removeborders(sin,siz)
sin = unpadimage(sin, [(siz+1)/2,(siz+1)/2,(siz-1)/2,(siz-1)/2]);
sin = padarray(sin, [(siz+1)/2,(siz+1)/2],0,'pre');
sout = padarray(sin, [(siz-1)/2,(siz-1)/2],0,'post');
return
