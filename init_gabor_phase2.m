function [fSiz,filters,c1OL,numSimpleFilters] = init_gabor_phase(rot, RF_siz, Div,numPhases)
% function init_gabor(rot, RF_siz, Div)
% Thomas R. Serre
% Feb. 2003

c1OL             = 2;
numFilterSizes   = length(RF_siz);%��ͬ��С�˲���������scale
numSimpleFilters = length(rot);%�������
numFilters       = numFilterSizes*numSimpleFilters;%��ͬ�߶Ⱥͷ�����˲�������
fSiz             = zeros(numFilters,1);	% vector with filter sizes
filters          = zeros(max(RF_siz)^2,numFilters,numPhases);

lambda = RF_siz*2./Div; %gabor�˲����Ĳ���
sigma  = lambda.*0.8;%gabor�˲�������Ч����
G      = 0.3;   % spatial aspect ratio: 0.23 < gamma < 0.92


for pp = 1:numPhases
    phase =  (pp - 1) / 2 * pi;%0,90,180,270
    
    for k = 1:numFilterSizes  %�ӵ�һ����С���˲�����ʼѭ��
        for r = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
            theta     = rot(r)*pi/180; %ת��Ϊ������
            filtSize  = RF_siz(k); %��k���˲����ĳ߶�
            center    = ceil(filtSize/2); %ȡ���������
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
                        E = exp(-(x^2+G^2*y^2)/(2*sigmaq))*cos(2*pi*x/lambda(k) + phase); %cvpr05�е�2�ڵ�һ����ʽ,gabor�˲���
                    end
                    f(j+center,i+center) = E;
                end
            end
            
            f = f - mean(mean(f));
            f = f ./ sqrt(sum(sum(f.^2)));
            p = numSimpleFilters*(k-1) + r;
            filters(1:filtSize^2,p,pp)=reshape(f,filtSize^2,1); %��11��11���˲�����Ϊ121��1
            fSiz(p)=filtSize;
            
            %             %         name = RF_siz(k)_rot(r);
            %      fname = ['E:\cognitive science\serre\5��2005 Object recognition with features inspired by visual cortex\NewResults','/','gabor','/',strcat(num2str(RF_siz(k)),'-',num2str(rot(r))),'.tif'];
            %         imwrite(f,fname,'resolution',300);
        end
    end
    
    clear f
    
end