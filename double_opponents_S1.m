function d = double_opponents_S1(im, filter,filters_gabor,numChannel,numSimpleFilters)
sigma = 0.125;
k = 1;

s1 = zeros(size(im, 1), size(im,2),numChannel,numSimpleFilters);

for ii = 1:numSimpleFilters
    for jj=1:numChannel
                
        for kk=1:3
            % 			tmp = conv2(im(:,:,kk), squeeze(filters_gabor(:,:,kk, jj, ii)), 'same');
            tmp = conv2padded(im(:,:,kk), squeeze(filters_gabor(:,:,kk, jj, ii)));
            s1(:,:,jj,ii)  = s1(:,:,jj,ii)  + tmp;%+R-center,-G-surround...
        end
        
%         s1(:,:,ii,jj) = ss1;
        
%         s1(ii,jj)(s1(ii,jj)<0) = 0; %rectification
        
    end
end
s1(s1<0) = 0;


%     tmp = 0;
% 	for ii = 1:4
% 		for jj=1:8
%             tmp = s1{jj}{ii} + alpha.s;
% 			tmp(tmp<0) = 0;
%             s1{jj}{ii} = tmp;
%             
% %             s1{jj}{ii} = sgolayfilt(s1{jj}{ii},2,9); %second-order Savitzky-Golay filtering to enhance local maxima
% 
% 		end
% 	end


	%normalise single opponent
% 	for jj = 1:8
% 		norm_s(:,:,jj) = sqrt(s1{jj}{1}.^2 +s1{jj}{2}.^2 + s1{jj}{3}.^2 + s1{jj}{4}.^2);
% 		norm_s(:,:,jj) = norm_s(:,:,jj) + ~norm_s(:,:,jj);
% 		%     norm_s(:,:,jj) = norm_s(:,:,jj) + 0.1;%avoid divide by zero 
% 	end
% 
% 
% 	for jj = 1:8
% 		for ii =1:4
% 			s1{jj}{ii} = s1{jj}{ii}./norm_s(:,:,jj);
% 		end
% 	end

% 	clear norm_s;  
    
%     % normaliztion
%    E_single = zeros(size(im, 1), size(im,2));
%    
%    for ii = 1:4
%        for jj = 1:8
%            E_single = E_single + s1{jj}{ii}.^2 * 0.25; % modified heeger's model
%        end
%    end
%    
%    for ii = 1:4
%        for jj = 1:8
%            s1{jj}{ii} = sqrt(k * s1{jj}{ii}.^2 ./ (sigma.^2 + E_single));
%        end
%    end



% 	w_sc = size(s1{1}{1},2);
% 	h_sc = size(s1{1}{1},1);
	d = zeros(size(im, 1), size(im,2),numChannel,numSimpleFilters);
    
for ii = 1:numSimpleFilters
    for jj=1:numChannel
        
			d_pos = conv2padded(s1(:,:,jj,ii),filter(:,:,1,ii));
			d_neg = conv2padded(s1(:,:,jj,ii),-1*filter(:,:,2,ii));
			d(:,:,jj,ii) = d_pos + d_neg;%+R-G-center,+G-R-surround...
            
%             d{jj}{ii}(d{jj}{ii}<0) = 0;
%                 smoothing
%             d{jj}{ii} = sgolayfilt(d{jj}{ii},2,9); %second-order Savitzky-Golay filtering to enhance local maxima
            
		end
	end
d(d<0) = 0;
	    
% for ii = 1:numSimpleFilters
%     for jj=1:numChannel
%             tmp  = d{ll}{jj} + alpha.d;
% 			tmp(tmp<0) = 0;%��У��Ļ���energy model��رյ��෴���ԵĽ�� %half-wave rectification
%             d{ll}{jj} = tmp; %squred
% 		end
% 	end

