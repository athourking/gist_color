function [fSiz,filters,c1OL,numSimpleFilters] = init_gabor_opp(rot, RF_siz, Div,numPhases)
% function init_gabor(rot, RF_siz, Div)
% Thomas R. Serre
% Feb. 2003

c1OL             = 2;
numFilterSizes   = length(RF_siz);%不同大小滤波器个数，scale
numSimpleFilters = length(rot);%方向个数
numFilters       = numFilterSizes*numSimpleFilters;%不同尺度和方向的滤波器个数
fSiz             = zeros(numFilters,1);	% vector with filter sizes
filters          = zeros(max(RF_siz)^2,numFilters,numPhases);

lambda = RF_siz*2./Div; %gabor滤波器的波长
sigma  = lambda.*0.8;%gabor滤波器的有效带宽
G      = 0.3;   % spatial aspect ratio: 0.23 < gamma < 0.92



filter1 = GenerateGabor(size_f, numSimpleFilters, G, lambda, sigma, 'positive', pCount);
filter2 = GenerateGabor(size_f, numSimpleFilters, G, lambda, sigma, 'negative', pCount);
% filter1 = GenerateGabor(size_f,4, 0.3, 5.6410, 4.5128, 'positive', pCount);
% filter2 = GenerateGabor(size_f,4, 0.3, 5.6410, 4.5128, 'negative', pCount);

%seperate to pos and neg parts
for p = 1:pCount
    for j=1:4
        filter{p}(:,:,1,j) = abs(filter1{p}(:,:,j));
        filter{p}(:,:,2,j) = abs(filter2{p}(:,:,j));
    end
end

%---------------------------gaussian----------------------
% size_kernel = 2;
filter3 = fspecial('gauss', size_f, size_kernel);
% filter3 = normalise(filter3);
filter3(filter3<0) = 0;



%%
if numChannel == 8
    weights = [1/sqrt(2),-1/sqrt(2),0; 2/sqrt(6),-1/sqrt(6),-1/sqrt(6); 1/sqrt(6),1/sqrt(6),-2/sqrt(6); 1/sqrt(3),1/sqrt(3),1/sqrt(3)]';
    %       weights = [1/sqrt(2),-1/sqrt(2),0; 2/sqrt(6),-1/sqrt(6),-1/sqrt(6); 1/sqrt(6),1/sqrt(6),-2/sqrt(6); 1/sqrt(6),-2/sqrt(6),1/sqrt(6)]';
    weights = [weights , -weights];
    
else if numChannel == 7
        weights = [1/sqrt(2),-1/sqrt(2),0; 2/sqrt(6),-1/sqrt(6),-1/sqrt(6); 1/sqrt(6),1/sqrt(6),-2/sqrt(6)]';
        w_bw = [1/sqrt(3),1/sqrt(3),1/sqrt(3)]';
        weights = [weights , -weights, w_bw];
        
    else if numChannel == 6
            weights = [1/sqrt(2),-1/sqrt(2),0; 2/sqrt(6),-1/sqrt(6),-1/sqrt(6); 1/sqrt(6),1/sqrt(6),-2/sqrt(6)]';
        else if numChannel == 10
                weights = [1/sqrt(2),-1/sqrt(2),0; 2/sqrt(6),-1/sqrt(6),-1/sqrt(6); 1/sqrt(6),1/sqrt(6),-2/sqrt(6); 1/sqrt(6),-2/sqrt(6),1/sqrt(6); 1/sqrt(3),1/sqrt(3),1/sqrt(3)]';
            end
        end
    end
end
% weights = [weights , -weights];


filter_all = {};

% for p = 1:pCount
for i=1:3
    
    for j=1:numChannel
        
        if weights(i,j) > 0
            hh = 1;
            
        elseif weights(i,j) < 0
            hh = 2;
            
        elseif weights(i,j) == 0
            hh = 1;
        end
        
        for p = 1:pCount
            filter_all{p}(:,:,i,j ,:) =  weights(i,j) * filter{p}(:,:,hh,:);
        end
        filters_gaussian(:,:,i,j) =  weights(i,j) * filter3;
    end
end

% filters_gabor = reshape(filter_all{p}, size_f, size_f, 3, numChannel*4);

for p = 1:pCount
    filters_gabor{p} = reshape(filter_all{p}, size_f, size_f, 3, numChannel*4);
    
    for i=1:3
        for j=1:numChannel*4
            nn = norm(filters_gabor{p}(:,:,i,j),2);
            if nn ~= 0
                filters_gabor{p}(:,:,i,j) = filters_gabor{p}(:,:,i,j)/nn;
            end
        end
    end
    
end




function fVals = GenerateGabor(rfCount, fCount, aspectRatio, lambda, sigma, gabor_sign, numPhases)


for pp = 1:numPhases
    phase =  (pp - 1) / 2 * pi;%0,90,180,270
    
    for k = 1:numFilterSizes  %从第一个大小的滤波器开始循环
        for r = 1:numSimpleFilters %从第一个方向开始循环
            theta     = rot(r)*pi/180; %转换为极坐标
            filtSize  = RF_siz(k); %第k个滤波器的尺度
            center    = ceil(filtSize/2); %取最近的整数
            filtSizeL = center-1;
            filtSizeR = filtSize-filtSizeL-1;
            sigmaq    = sigma(k)^2;
            
            for i = -filtSizeL:filtSizeR
                for j = -filtSizeL:filtSizeR
                    
                    if ( sqrt(i^2+j^2)>filtSize/2 )
                        E = 0;
                    else
                        x = i*cos(theta) - j*sin(theta);
                        y = i*sin(theta) + j*cos(theta);
                        E = exp(-(x^2+G^2*y^2)/(2*sigmaq))*cos(2*pi*x/lambda(k) + phase); %cvpr05中第2节第一个公式,gabor滤波器
                    end
                    f(j+center,i+center) = E;
                end
            end
            
            f = f - mean(mean(f));
            f = f ./ sqrt(sum(sum(f.^2)));
            p = numSimpleFilters*(k-1) + r;
            filters(1:filtSize^2,p,pp)=reshape(f,filtSize^2,1); %将11×11的滤波器变为121×1
            fSiz(p)=filtSize;
            
            %             %         name = RF_siz(k)_rot(r);
            %      fname = ['E:\cognitive science\serre\5、2005 Object recognition with features inspired by visual cortex\NewResults','/','gabor','/',strcat(num2str(RF_siz(k)),'-',num2str(rot(r))),'.tif'];
            %         imwrite(f,fname,'resolution',300);
        end
    end
    
    clear f
    
end

return