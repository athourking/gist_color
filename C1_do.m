function g = C1_do(stim, filter,filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,numChannel,numPhases)

USECONV2 = 1; %should be faster if 1
k = 1;
sigma = 0.125;

USE_NORMXCORR_INSTEAD = 0;
if(nargin < 13)
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
% numPhases = 4;
% nFilts_gray = length(fSiz_gray);
% for i = 1:nFilts_gray
%     sqfilter{i} = reshape(filters_gray(1:(fSiz_gray(i)^2),i),fSiz_gray(i),fSiz_gray(i));
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
% for i = 1:length(uFiltSizes_gray)
%     s1Norm{uFiltSizes_gray(i)} = (sumfilter(sqim_gray,(uFiltSizes_gray(i)-1)/2)).^0.5;
%     %avoid divide by zero
%     s1Norm{uFiltSizes_gray(i)} = s1Norm{uFiltSizes_gray(i)} + ~s1Norm{uFiltSizes_gray(i)};
% end

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
            s1{iPhase}{iBand}{iScale} = single_opponents_S1(stim,filters{sc}{iPhase},numChannel,numSimpleFilters);%double opponent simple cell
             if(~INCLUDEBORDERS)
                 
                  for jj=1:numChannel
                      for ii = 1:numSimpleFilters
                          iUFilterIndex = iUFilterIndex+1;
                 s1{iPhase}{iBand}{iScale}(:,:,jj,ii) = removeborders(s1{iPhase}{iBand}{iScale}(:,:,jj,ii),fSiz(iUFilterIndex+numChannel*numSimpleFilters*(sc-1)));
            
                      end
                  end
             end
            
            s1{iPhase}{iBand}{iScale} = im2double(s1{iPhase}{iBand}{iScale});            
        end
      %normalization
        E_single = zeros(size(stim, 1), size(stim,2),numChannel);
   
      for jj=1:numChannel
       for ii = 1:numSimpleFilters
        
        for iPhase = 1:numPhases
          E_single(:,:,jj) = E_single(:,:,jj) + 0.25 *  s1{iPhase}{iBand}{iScale} (:,:,jj,ii).^2;
        end
        
       end
       
      end

   
   for iPhase = 1:numPhases
for ii = 1:numSimpleFilters
    for jj=1:numChannel
           s1{iPhase}{iBand}{iScale}(:,:,jj,ii) = sqrt(k * s1{iPhase}{iBand}{iScale}(:,:,jj,ii).^2 ./ (sigma.^2 + E_single(:,:,jj)));
    end
end
   end
   
  for iPhase = 1:numPhases
d{iPhase}{iBand}{iScale} = double_opponents_d1(stim, s1{iPhase}{iBand}{iScale},filter{sc}{1},numChannel,numSimpleFilters);
   
    ss1{iBand}{iScale} =  ss1{iBand}{iScale} + d{iPhase}{iBand}{iScale} / numPhases;
   end
   
        
%         ss1{iBand}{iScale} = sqrt(s1{1}{iBand}{iScale}.^2 + s1{2}{iBand}{iScale}.^2);
        
         for iFilt = 1:numSimpleFilters %
            for jj = 1:numChannel/2
                 sss1{iBand}{iScale}(:,:,jj,iFilt) = sqrt( ss1{iBand}{iScale}(:,:,jj,iFilt).^2 + ss1{iBand}{iScale}(:,:,jj+numChannel/2,iFilt).^2);
            end
         end
 
    end
end

clear  s1 ss1


% ss1 = sqrt(s1{1}{iBand}{iScale}(:,:,iFilt,channel) .^2 + s1{2}{iBand}{iScale} .^2);
% ss1 = sqrt(s1(:,:,:,:,1,:,:) .^2 + s1(:,:,:,:,2,:,:) .^2);

c1 = {};
for iBand = 1:numScaleBands
        for jj=1:numChannel/2
            for iFilt = 1:numSimpleFilters
%                 for iPhase = 1:numPhases
                c1{iBand}(:,:,jj,iFilt) = zeros(size(sss1{iBand}{1}(:,:,jj,iFilt)));
                for iScale = 1:length(ScalesInThisBand{iBand});
                    c1{iBand}(:,:,jj,iFilt) = max(c1{iBand}(:,:,jj,iFilt),sss1{iBand}{iScale}(:,:,jj,iFilt));
                end
%             end
            
        end
    end
end

    
    
%   (2) pool over local neighborhood
for iBand = 1:numScaleBands %��ÿ���߶ȴ�ʼѭ��
    poolRange = (c1SpaceSS(iBand)); %����������С10
    for jj=1:numChannel/2
        for iFilt = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
%             for iPhase = 1:numPhases
                c1{iBand}(:,:,jj,iFilt) = maxfilter(c1{iBand}(:,:,jj,iFilt),[0 0 poolRange-1 poolRange-1]); %��ÿ������ȡ���Ľϴ��s1������̬ѧ���ʹ���
%         end
    end
end
end

% 
%   (3) subsample
for iBand = 1:numScaleBands %�ӵ�һ���߶ȴ�ʼѭ��
    sSS=ceil(c1SpaceSS(iBand)/c1OL);
    clear T;
    
    for jj=1:numChannel/2
    for iFilt = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
%     for iPhase = 1:numPhases
        T(:,:,jj,iFilt) = c1{iBand}(1:sSS:end,1:sSS:end,jj,iFilt); %�������õ�ÿ���������Ӧ
        %     T(:,:,iFilt) = c1{iBand}(1:2:end,1:2:end,iFilt); %�������õ�ÿ���������Ӧ
%     end
    end
    end
    c1{iBand} = T;

end

g = gistGabor_C1_opp(c1,numScaleBands,numChannel/2,numSimpleFilters);


function sout = removeborders(sin,siz)
sin = unpadimage(sin, [(siz+1)/2,(siz+1)/2,(siz-1)/2,(siz-1)/2]);
sin = padarray(sin, [(siz+1)/2,(siz+1)/2],0,'pre');
sout = padarray(sin, [(siz-1)/2,(siz-1)/2],0,'post');
return
