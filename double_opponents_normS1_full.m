function d = double_opponents_normS1(im, filter,filters_gabor,numChannel,numSimpleFilters)
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
    
    % normaliztion
   E_single = zeros(size(im, 1), size(im,2),numChannel);
   
for ii = 1:numSimpleFilters
    for jj=1:numChannel
%            E_single = E_single + s1{jj}{ii}.^2 * 0.25; % modified heeger's model
%           E_single(:,:,jj) = E_single(:,:,jj) + 0.25 * ( s1(:,:,jj,ii).^n(2) +  s1(:,:,jj+numChannel,ii).^2 + s1(:,:,jj+numChannle*2,ii).^2 +  s1(:,:,jj+numChannel*3,ii).^2);
          E_single(:,:,jj) = E_single(:,:,jj) + 0.25 *  s1(:,:,jj,ii).^2;
    end
end
   
   
for ii = 1:numSimpleFilters
    for jj=1:numChannel
           s1(:,:,jj,ii) = sqrt(k * s1(:,:,jj,ii).^2 ./ (sigma.^2 + E_single(:,:,jj)));
    end
end

% 
%         for pp = 1:pCount
%                 for jj = 1:numChannel
%                     E_lum(:,:,jj,pp) = E_lum(:,:,jj,pp) + 0.25 * (soRes_lum{ii,mm,nn}(:,:,jj,pp).^n(2) + soRes_lum{ii,mm,nn}(:,:,jj+numChannel,pp).^n(2) + soRes_lum{ii,mm,nn}(:,:,jj+numChannel*2,pp).^n(2) + soRes_lum{ii,mm,nn}(:,:,jj+numChannel*3,pp).^n(2));
%                     E_equil(:,:,jj,pp) = E_equil(:,:,jj,pp) + 0.25 * (soRes_equil{ii,mm,nn}(:,:,jj,pp).^n(2) + soRes_equil{ii,mm,nn}(:,:,jj+numChannel,pp).^n(2) + soRes_equil{ii,mm,nn}(:,:,jj+numChannel*2,pp).^n(2) + soRes_equil{ii,mm,nn}(:,:,jj+numChannel*3,pp).^n(2));
%                 end
%             end
%             for pp = 1:pCount
%                 for ss = 1:4
%                     for jj = 1+numChannel*(ss-1):numChannel+numChannel*(ss-1)
%                         soRes_lum{ii,mm,nn}(:,:,jj,pp) = sqrt(k .* soRes_lum{ii,mm,nn}(:,:,jj,pp).^n(2) ./ (sigma(ii).^n(2) + E_lum(:,:,jj-numChannel*(ss-1),pp)));
%                         soRes_equil{ii,mm,nn}(:,:,jj,pp) = sqrt(k .* soRes_equil{ii,mm,nn}(:,:,jj,pp).^n(2) ./ (sigma(ii).^n(2) + E_equil(:,:,jj-numChannel*(ss-1),pp)));
%                     end
%                 end
%             end
   
   
   
   

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
d = abs(d);
	    
% for ii = 1:numSimpleFilters
%     for jj=1:numChannel
%             tmp  = d{ll}{jj} + alpha.d;
% 			tmp(tmp<0) = 0;%��У��Ļ���energy model��رյ��෴���ԵĽ�� %half-wave rectification
%             d{ll}{jj} = tmp; %squred
% 		end
% 	end

