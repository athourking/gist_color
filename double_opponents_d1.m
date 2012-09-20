function d = double_opponents_d1(im, s1,filter,numChannel,numSimpleFilters)


	d = zeros(size(im, 1), size(im,2),numChannel,numSimpleFilters);
    
for ii = 1:numSimpleFilters
    for jj=1:numChannel
        
			d_pos = conv2(s1(:,:,jj,ii),filter(:,:,1,ii),'same');
			d_neg = conv2(s1(:,:,jj,ii),-1*filter(:,:,2,ii),'same');
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

