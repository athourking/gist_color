function g = gaborV1(stim, filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,INCLUDEBORDERS)
%function [c1,s1] = C1(stim, filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,INCLUDEBORDERS)
%
%  A matlab implementation of the C1 code originally by Max Riesenhuber
%  and Thomas Serre.
%  Adapted by Stanley Bileschi
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
numScaleBands=length(c1ScaleSS)-1;  % convention: last element in c1ScaleSS is max index + 1 尺度带的个数
numScales=c1ScaleSS(end)-1; %尺度的个数
%   last index in scaleSS contains scale index where next band would start, i.e., 1 after highest scale!!
numSimpleFilters=floor(length(fSiz)/numScales); %4个方向
% numSimpleFilters=4;
for iBand = 1:numScaleBands
    ScalesInThisBand{iBand} = c1ScaleSS(iBand):(c1ScaleSS(iBand+1) -1); %在band 2中有2种尺度
end

% Rebuild all filters (of all sizes)
%%%%%%%%
nFilts = length(fSiz); %2种尺度，每个尺度有4个方向，总共8个
% nFilts = size(filters,2); %2种尺度，每个尺度有4个方向，总共8个
for iPhase = 1:2
    for i = 1:nFilts
        sqfilter{i}{iPhase} = reshape(filters(1:(fSiz(i)^2),i,iPhase),fSiz(i),fSiz(i)); %每个尺度在每个方向上的滤波器
        if USECONV2
            sqfilter{i}{iPhase} = sqfilter{i}{iPhase}(end:-1:1,end:-1:1); %flip in order to use conv2 instead of imfilter (%bug_fix 6/28/2007);
        end
    end
end

% Calculate all filter responses (s1)
%%%%%%%%
sqim = stim.^2;
% iUFilterIndex = 0;
% precalculate the normalizations for the usable filter sizes
uFiltSizes = unique(fSiz); %两个尺度
for i = 1:length(uFiltSizes) %从第一个尺度开始循环
    s1Norm{uFiltSizes(i)} = (sumfilter(sqim,(uFiltSizes(i)-1)/2)).^0.5;
    %avoid divide by zero
    s1Norm{uFiltSizes(i)} = s1Norm{uFiltSizes(i)} + ~s1Norm{uFiltSizes(i)};
end

for iPhase = 1:2
    iUFilterIndex = 0;
    
    for iBand = 1:numScaleBands %从第一个尺度带开始循环
        for iScale = 1:length(ScalesInThisBand{iBand}) %从第一个尺度开始循环
            for iFilt = 1:numSimpleFilters %从第一个方向开始循环
                iUFilterIndex = iUFilterIndex+1;
                if ~USECONV2
                    s1{iBand}{iScale}{iFilt} = abs(imfilter(stim,sqfilter{iUFilterIndex},'symmetric','same','corr'));
                    
                    if(~INCLUDEBORDERS)
                        s1{iBand}{iScale}{iFilt} = removeborders(s1{iBand}{iScale}{iFilt},fSiz(iUFilterIndex)); %用face数据集时，提取C2时这里出错，也就是当尺度大小为21时出错，主要原因是：这个时候的s1是21*21，unpagimage函数执行后返回[]
                    end
                    s1{iBand}{iScale}{iFilt} = im2double(s1{iBand}{iScale}{iFilt}) ./ s1Norm{fSiz(iUFilterIndex)}; %得到s1，每个尺度下每个方向滤波后图像
                else %not 100% compatible but 20% faster at least
                    %                 s1{iBand}{iScale}{iFilt} = abs(conv2(stim,sqfilter{iUFilterIndex},'same'));
                    %                 s1(:,:,iBand,iScale,iFilt) = abs(conv2(stim,sqfilter{iUFilterIndex},'same'));
                    s1(:,:,iBand,iScale,iFilt,iPhase) = abs(conv2padded(stim,sqfilter{iUFilterIndex}{iPhase}));
                    
%                     if(~INCLUDEBORDERS)
%                         %                     s1{iBand}{iScale}{iFilt} = removeborders(s1{iBand}{iScale}{iFilt},fSiz(iUFilterIndex));
%                         s1(:,:,iBand,iScale,iFilt,iPhase) = removeborders(s1(:,:,iBand,iScale,iFilt,iPhase),fSiz(iUFilterIndex));
%                     end
%                     s1{iBand}{iScale}{iFilt} = im2double(s1{iBand}{iScale}{iFilt}) ./ s1Norm{fSiz(iUFilterIndex)};
                    s1(:,:,iBand,iScale,iFilt,iPhase) = im2double(s1(:,:,iBand,iScale,iFilt,iPhase)) ./ s1Norm{fSiz(iUFilterIndex)};
                end
                
            end
        end
    end
    
end
ss1 = sqrt(s1(:,:,:,:,:,1).^2 + s1(:,:,:,:,:,2).^2);
% ss1 = s1(:,:,:,:,:,1) + s1(:,:,:,:,:,2);
% clear s1
% s1 = ss1;

% Calculate local pooling (c1)
%%%%%%%%

%   (1) pool over scales within band
for iBand = 1:numScaleBands %从第一个尺度带开始循环
    for iFilt = 1:numSimpleFilters %从第一个方向开始循环
        c1{iBand}(:,:,iFilt) = zeros(size(ss1(:,:,iBand,1,iFilt)));
        for iScale = 1:length(ScalesInThisBand{iBand}); %从第一个尺度开始循环
            c1{iBand}(:,:,iFilt) = max(c1{iBand}(:,:,iFilt),ss1(:,:,iBand,iScale,iFilt)); %取出每个方向较大的不同尺度下的s1响应
        end
    end
end

%   (2) pool over local neighborhood
for iBand = 1:numScaleBands %从每个尺度带开始循环
    poolRange = (c1SpaceSS(iBand)); %将采样网格大小10
    for iFilt = 1:numSimpleFilters %从第一个方向开始循环
        c1{iBand}(:,:,iFilt) = maxfilter(c1{iBand}(:,:,iFilt),[0 0 poolRange-1 poolRange-1]); %对每个方向取出的较大的s1进行形态学膨胀处理
    end
end

% 
% %   (3) subsample
% for iBand = 1:numScaleBands %从第一个尺度带开始循环
%     sSS=ceil(c1SpaceSS(iBand)/c1OL);
%     clear T;
%     for iFilt = 1:numSimpleFilters %从第一个方向开始循环
%         T(:,:,iFilt) = c1{iBand}(1:sSS:end,1:sSS:end,iFilt); %将采样后得到每个方向的响应
%         %     T(:,:,iFilt) = c1{iBand}(1:2:end,1:2:end,iFilt); %将采样后得到每个方向的响应
%     end
%     c1{iBand} = T;
% end


g = gistGabor_V1(c1,numScaleBands,numSimpleFilters);






function sout = removeborders(sin,siz)
sin = unpadimage(sin, [(siz+1)/2,(siz+1)/2,(siz-1)/2,(siz-1)/2]);
sin = padarray(sin, [(siz+1)/2,(siz+1)/2],0,'pre');
sout = padarray(sin, [(siz-1)/2,(siz-1)/2],0,'post');


