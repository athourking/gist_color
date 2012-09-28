function [fSiz,filters,c1OL,numSimpleFilters] = init_gabor_phase(rot, RF_siz, Div,numPhases)
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