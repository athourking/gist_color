function [c1,c1_gray] = C1_dogray_mode(stim_gray,filters_gray,fSiz_gray,stim, filter,filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,numChannel,numPhases)

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
nFilts_gray = length(fSiz_gray);

for iPhase = 1:numPhases
    for i = 1:nFilts_gray
        sqfilter{i}{iPhase} = reshape(filters_gray(1:(fSiz_gray(i)^2),i,iPhase),fSiz_gray(i),fSiz_gray(i)); %ÿ���߶���ÿ�������ϵ��˲���
        if USECONV2
            sqfilter{i}{iPhase} = sqfilter{i}{iPhase}(end:-1:1,end:-1:1); %flip in order to use conv2 instead of imfilter (%bug_fix 6/28/2007);
        end
    end
end
% Calculate all filter responses (s1)
%%%%%%%%
sqim = stim.^2;
sqim_gray = stim_gray.^2;
iUFilterIndex = 0;
iUFilterIndex_gray = 0;
% precalculate the normalizations for the usable filter sizes
uFiltSizes = unique(fSiz);
uFiltSizes_gray = unique(fSiz_gray);
for i = 1:length(uFiltSizes_gray)
    s1Norm{uFiltSizes_gray(i)} = (sumfilter(sqim_gray,(uFiltSizes_gray(i)-1)/2)).^0.5;
    %avoid divide by zero
    s1Norm{uFiltSizes_gray(i)} = s1Norm{uFiltSizes_gray(i)} + ~s1Norm{uFiltSizes_gray(i)};
end
% ss1 = {};
% s1 = [];
%
% for iPhase = 1:numPhases
%     iUFilterIndex = 0;
% iUFilterIndex = 0;
for iBand = 1:numScaleBands
    for iScale = 1:length(ScalesInThisBand{iBand})
        
        ss1{iBand}{iScale} = zeros(size(stim,1),size(stim,2),numChannel,numSimpleFilters);
        ss1_gray{iBand}{iScale} = zeros(size(stim,1),size(stim,2),numSimpleFilters);

        sc = (iBand-1)*length(ScalesInThisBand{iBand}) + iScale;

        %% do
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
         
         
         %% gray
        for iPhase = 1:numPhases
             iUFilterIndex_gray = 0;
      
                    %% gray
           for iFilt = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
                iUFilterIndex_gray = iUFilterIndex_gray+1;
                
                if ~USECONV2
                    s1{iBand}{iScale}{iFilt} = abs(imfilter(stim,sqfilter{iUFilterIndex_gray},'symmetric','same','corr'));
                    
                    if(~INCLUDEBORDERS)
                        s1{iBand}{iScale}{iFilt} = removeborders(s1{iBand}{iScale}{iFilt},fSiz(iUFilterIndex)); %��face��ݼ�ʱ����ȡC2ʱ������?Ҳ���ǵ��߶ȴ�СΪ21ʱ���?��Ҫԭ���ǣ����ʱ���s1��21*21��unpagimage����ִ�к󷵻�[]
                    end
                    s1{iBand}{iScale}{iFilt} = im2double(s1{iBand}{iScale}{iFilt}) ./ s1Norm{fSiz(iUFilterIndex)}; %�õ�s1��ÿ���߶���ÿ�������˲���ͼ��
                else %not 100% compatible but 20% faster at least
                    s1_gray{iPhase}{iBand}{iScale}(:,:,iFilt) = abs(conv2(stim_gray,sqfilter{iUFilterIndex_gray}{iPhase},'same'));
                    s1_gray{iPhase}{iBand}{iScale}(:,:,iFilt) = im2double(s1_gray{iPhase}{iBand}{iScale}(:,:,iFilt)) ./ s1Norm{fSiz_gray(iUFilterIndex_gray)};
                end
                
           end
           
            ss1_gray{iBand}{iScale} = ss1_gray{iBand}{iScale} + s1_gray{iPhase}{iBand}{iScale}./numPhases;
        end
 

    end
end


% ss1 = sqrt(s1{1}{iBand}{iScale}(:,:,iFilt,channel) .^2 + s1{2}{iBand}{iScale} .^2);
% ss1 = sqrt(s1(:,:,:,:,1,:,:) .^2 + s1(:,:,:,:,2,:,:) .^2);

c1 = {};

% c1_gray = {};

for iBand = 1:numScaleBands
        for iFilt = 1:numSimpleFilters
            
            c1_gray{iBand}(:,:,iFilt) = zeros(size(ss1_gray{iBand}{1}(:,:,iFilt)));
            for iScale = 1:length(ScalesInThisBand{iBand});
                c1_gray{iBand}(:,:,iFilt) = max(c1_gray{iBand}(:,:,iFilt),ss1_gray{iBand}{iScale}(:,:,iFilt)); %ȡ��ÿ������ϴ�Ĳ�ͬ�߶��µ�s1��Ӧ
            end

            for jj=1:numChannel/2
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
    for iFilt = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
        c1_gray{iBand}(:,:,iFilt) = maxfilter(c1_gray{iBand}(:,:,iFilt),[0 0 poolRange-1 poolRange-1]); %��ÿ������ȡ���Ľϴ��s1������̬ѧ���ʹ���
        
        for jj=1:numChannel/2
            c1{iBand}(:,:,jj,iFilt) = maxfilter(c1{iBand}(:,:,jj,iFilt),[0 0 poolRange-1 poolRange-1]); %��ÿ������ȡ���Ľϴ��s1������̬ѧ���ʹ���
        end
      
    end
end

% 
%   (3) subsample
for iBand = 1:numScaleBands %�ӵ�һ���߶ȴ�ʼѭ��
    sSS=ceil(c1SpaceSS(iBand)/c1OL);
    clear T;
    clear T_gray;
    for iFilt = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
            T_gray(:,:,iFilt) = c1_gray{iBand}(1:sSS:end,1:sSS:end,iFilt); %�������õ�ÿ���������Ӧ
            
            for jj=1:numChannel/2
                T(:,:,jj,iFilt) = c1{iBand}(1:sSS:end,1:sSS:end,jj,iFilt); %�������õ�ÿ���������Ӧ
            end
            
            
    end
    c1{iBand} = T;
    c1_gray{iBand} = T_gray;
end

function sout = removeborders(sin,siz)
sin = unpadimage(sin, [(siz+1)/2,(siz+1)/2,(siz-1)/2,(siz-1)/2]);
sin = padarray(sin, [(siz+1)/2,(siz+1)/2],0,'pre');
sout = padarray(sin, [(siz-1)/2,(siz-1)/2],0,'post');
return
