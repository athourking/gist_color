function doLocalPb_orient = doLocalPb_normPhase_sum_smooth_norm_gray(stim, filter,filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,numChannel,numPhases)


%% mPb
% [tg1, tg2, tg3, textons] = multiscalePb_tg(stim, 1);


stim = 2 * stim - 1;
USECONV2 = 1; %should be faster if 1
k = 1;
sigma = 0.125;


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
% numPhases = 4;
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
        ss1{iBand}{iScale} = zeros(size(stim,1),size(stim,2),numChannel,numSimpleFilters);
        
        for iPhase = 1:numPhases
            s1{iPhase}{iBand}{iScale} = single_opponents_S1(stim,filters{sc}{iPhase},numChannel,numSimpleFilters);%double opponent simple cell
            s1{iPhase}{iBand}{iScale} = im2double(s1{iPhase}{iBand}{iScale});
        end
        
        %normalization
        E_single = zeros(size(stim, 1), size(stim,2),numChannel);
        
        for ii = 1:numSimpleFilters
            for jj=1:numChannel
                
                for iPhase = 1:numPhases
                    %            E_single = E_single + s1{jj}{ii}.^2 * 0.25; % modified heeger's model
                    %           E_single(:,:,jj) = E_single(:,:,jj) + 0.25 * ( s1(:,:,jj,ii).^n(2) +  s1(:,:,jj+numChannel,ii).^2 + s1(:,:,jj+numChannle*2,ii).^2 +  s1(:,:,jj+numChannel*3,ii).^2);
                    E_single(:,:,jj) = E_single(:,:,jj) +s1{iPhase}{iBand}{iScale} (:,:,jj,ii).^2./(numSimpleFilters+numPhases);
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
        
        for iFilt = 1:numSimpleFilters %
            for jj = 1:numChannel/2
                sss1{iBand}{iScale}(:,:,jj,iFilt) = max( ss1{iBand}{iScale}(:,:,jj,iFilt), ss1{iBand}{iScale}(:,:,jj+numChannel/2,iFilt));
            end
        end
%        sss1{iBand}{iScale} = shiftdim(shiftdim(max(ss1{iBand}{iScale},[],3),3),1);%max color channels
      
    end
end


c1 = {};
doLocalPb_orient =  zeros(size(stim,1),size(stim,2),numSimpleFilters);

for iBand = 1:numScaleBands
    for jj=1:numChannel/2
        for iFilt = 1:numSimpleFilters
            
            c1{iBand}(:,:, jj,iFilt) = zeros(size(sss1{iBand}{1}(:,:,jj,iFilt)));
%             c1{iBand}(:,:,iFilt) = zeros(size(sss1{iBand}{1}(:,:,iFilt)));
%             doLocalPb_orient(:,:,iFilt) = zeros(size(sss1{iBand}{1}(:,:,iFilt)));
            
            for iScale = 1:length(ScalesInThisBand{iBand});
                c1{iBand}(:,:,jj,iFilt) = max(c1{iBand}(:,:, jj,iFilt),sss1{iBand}{iScale}(:,:,jj,iFilt));
%                 c1{iBand}(:,:,iFilt) = max(c1{iBand}(:,:,iFilt),sss1{iBand}{iScale}(:,:,iFilt));%max scales
            end
            
%             doLocalPb_orient(:,:,iFilt) = max(doLocalPb_orient(:,:,iFilt),c1{iBand}(:,:,iFilt));%max bands
%             doLocalPb_orient(:,:,iFilt) = doLocalPb_orient(:,:,iFilt) + c1{iBand}(:,:,jj,iFilt);%max bands
            
        end
    end
end

% weights = [0  0  0.0039  0.0050  0.0058  0.0069  0.0040  0.0044  0.0049  0.0024  0.0027  0.0170  0.0074];


% [bg1, bg2, bg3, cga1, cga2, cga3, cgb1, cgb2, cgb3, tg1, tg2, tg3, textons] = det_mPb(stim);
% smooth cues

% load('/users/jz7/data/recent/Contour_Detection_and_Image_Segmentation/BSR_full/BSR/grouping/data/tg1.mat');
% load('/users/jz7/data/recent/Contour_Detection_and_Image_Segmentation/BSR_full/BSR/grouping/data/tg2.mat');
% load('/users/jz7/data/recent/Contour_Detection_and_Image_Segmentation/BSR_full/BSR/grouping/data/tg3.mat');


gtheta = [1.5708    1.1781    0.7854    0.3927   0    2.7489    2.3562    1.9635];
tic;
filters = make_filters([9 5 10 20], gtheta);

for iFilt = 1 : numSimpleFilters
    for iBand = 1:numScaleBands
        for jj=1:numChannel/2
            c1{iBand}(:,:,jj,iFilt) = fitparab(c1{iBand}(:,:,jj,iFilt),9,9/4,gtheta(iFilt),filters{1,iFilt});
            
%     bg1(:,:,o) = fitparab(bg1(:,:,o),3,3/4,gtheta(o),filters{1,o});
%     bg2(:,:,o) = fitparab(bg2(:,:,o),5,5/4,gtheta(o),filters{2,o});
%     bg3(:,:,o) = fitparab(bg3(:,:,o),10,10/4,gtheta(o),filters{3,o});

%     cga1(:,:,o) = fitparab(cga1(:,:,o),5,5/4,gtheta(o),filters{2,o});
%     cga2(:,:,o) = fitparab(cga2(:,:,o),10,10/4,gtheta(o),filters{3,o});
%     cga3(:,:,o) = fitparab(cga3(:,:,o),20,20/4,gtheta(o),filters{4,o});
% 
%     cgb1(:,:,o) = fitparab(cgb1(:,:,o),5,5/4,gtheta(o),filters{2,o});
%     cgb2(:,:,o) = fitparab(cgb2(:,:,o),10,10/4,gtheta(o),filters{3,o});
%     cgb3(:,:,o) = fitparab(cgb3(:,:,o),20,20/4,gtheta(o),filters{4,o});
        end
    end

%     tg1(:,:,iFilt) = fitparab(tg1(:,:,iFilt),5,5/4,gtheta(iFilt),filters{2,iFilt});
%     tg2(:,:,iFilt) = fitparab(tg2(:,:,iFilt),10,10/4,gtheta(iFilt),filters{3,iFilt});
%     tg3(:,:,iFilt) = fitparab(tg3(:,:,iFilt),20,20/4,gtheta(iFilt),filters{4,iFilt});

end
fprintf('Cues smoothing: %g\n', toc);


% compute mPb at full scale
% mPb_all = zeros(size(tg1));
for iFilt = 1 : numSimpleFilters
     for iBand = 1:numScaleBands
        for jj=1:numChannel/2
            
            doLocalPb_orient(:,:,iFilt) = doLocalPb_orient(:,:,iFilt) + normaliseminmax(c1{iBand}(:,:,jj,iFilt));%max bands
%     l1 = weights(1)*bg1(:, :, o);
%     l2 = weights(2)*bg2(:, :, o);
%     l3 = weights(3)*bg3(:, :, o);
% 
%     a1 = weights(4)*cga1(:, :, o);
%     a2 = weights(5)*cga2(:, :, o);
%     a3 = weights(6)*cga3(:, :, o);
% 
%     b1 = weights(7)*cgb1(:, :, o);
%     b2 = weights(8)*cgb2(:, :, o);
%     b3 = weights(9)*cgb3(:, :, o);

%     t1 = weights(10)*tg1(:, :, iFilt);
%     t2 = weights(11)*tg2(:, :, iFilt);
%     t3 = weights(12)*tg3(:, :, iFilt);
        end
     end
%     
    
%     doLocalPb_orient(:,:,iFilt) = normaliseminmax(doLocalPb_orient(:,:,iFilt));
    
%     t1 = normaliseminmax(tg1(:, :, iFilt));
%     t2 = normaliseminmax(tg2(:, :, iFilt));
%     t3 = normaliseminmax(tg3(:, :, iFilt));

%     doLocalPb_orient(:, :, iFilt) = doLocalPb_orient(:,:,iFilt) + t1 + t2 + t3;

end


% % non-maximum suppression
% doPb_nmax = nonmax_channels(doLocalPb_orient);
% doPb_nmax = max(0, min(1, 1.2*doPb_nmax));
% %An optional non-maximum suppression step [22] produces thinned, real-valued contours.   
% 
% 
% %% outputs
% doLocalPb = max(doLocalPb_orient, [], 3);
% 
% doLocalPb_thin = doLocalPb .* (doPb_nmax>0.05);
% doLocalPb_thin = doLocalPb_thin .* bwmorph(doLocalPb_thin, 'skel', inf);
% 
% % if ~strcmp(outFile,''), save(outFile,'doLocalPb_thin', 'doLocalPb_orient'); end
% 
% 
% 
% % %   (2) pool over local neighborhood
% % for iBand = 1:numScaleBands %��ÿ���߶ȴ�ʼѭ��
% %     poolRange = (c1SpaceSS(iBand)); %����������С10
% %     for jj=1:numChannel/2
% %         for iFilt = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
% %         c1{iBand}(:,:,jj,iFilt) = maxfilter(c1{iBand}(:,:,jj,iFilt),[0 0 poolRange-1 poolRange-1]); %��ÿ������ȡ���Ľϴ��s1������̬ѧ���ʹ���
% %         end
% %     end
% % end
% 
% % 
% % %   (3) subsample
% % for iBand = 1:numScaleBands %�ӵ�һ���߶ȴ�ʼѭ��
% %     sSS=ceil(c1SpaceSS(iBand)/c1OL);
% %     clear T;
% %     for iFilt = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
% %         T(:,:,iFilt) = c1{iBand}(1:sSS:end,1:sSS:end,iFilt); %�������õ�ÿ����������?%         %     T(:,:,iFilt) = c1{iBand}(1:2:end,1:2:end,iFilt); %�������õ�ÿ����������?%     end
% %     c1{iBand} = T;
% % end



%%
function filters = make_filters(radii, gtheta)

d = 2; 

filters = cell(numel(radii), numel(gtheta));
for r = 1:numel(radii),
    for t = 1:numel(gtheta),
        
        ra = radii(r);
        rb = ra / 4;
        theta = gtheta(t);
        
        ra = max(1.5, ra);
        rb = max(1.5, rb);
        ira2 = 1 / ra^2;
        irb2 = 1 / rb^2;
        wr = floor(max(ra, rb));
        wd = 2*wr+1;
        sint = sin(theta);
        cost = cos(theta);
        
        % 1. compute linear filters for coefficients
        % (a) compute inverse of least-squares problem matrix
        filt = zeros(wd,wd,d+1);
        xx = zeros(2*d+1,1);
        for u = -wr:wr,
            for v = -wr:wr,
                ai = -u*sint + v*cost; % distance along major axis
                bi = u*cost + v*sint; % distance along minor axis
                if ai*ai*ira2 + bi*bi*irb2 > 1, continue; end % outside support
                xx = xx + cumprod([1;ai+zeros(2*d,1)]);
            end
        end
        A = zeros(d+1,d+1);
        for i = 1:d+1,
            A(:,i) = xx(i:i+d);
        end
        
        % (b) solve least-squares problem for delta function at each pixel
        for u = -wr:wr,
            for v = -wr:wr,
                ai = -u*sint + v*cost; % distance along major axis
                bi = u*cost + v*sint; % distance along minor axis
                if (ai*ai*ira2 + bi*bi*irb2) > 1, continue; end % outside support
                yy = cumprod([1;ai+zeros(d,1)]);
                filt(v+wr+1,u+wr+1,:) = A\yy;
            end
        end
        
        filters{r,t}=filt;
    end
end






function sout = removeborders(sin,siz)
sin = unpadimage(sin, [(siz+1)/2,(siz+1)/2,(siz-1)/2,(siz-1)/2]);
sin = padarray(sin, [(siz+1)/2,(siz+1)/2],0,'pre');
sout = padarray(sin, [(siz-1)/2,(siz-1)/2],0,'post');

return