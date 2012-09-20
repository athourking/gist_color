function s1 = single_opponents(im, filters_gabor,numChannel)

sigma = 0.125;
k = 1;

	s1 = {};
    
for ii = 1:4
	for jj=1:numChannel
		s1{jj}{ii} = zeros(size(im, 1), size(im,2));
            
		for kk=1:3 
			tmp = conv2(im(:,:,kk), squeeze(filters_gabor(:,:,kk, jj, ii)), 'same');
			s1{jj}{ii} = s1{jj}{ii} + tmp;%+R-center,-G-surround...
        end
            
        s1{jj}{ii}(s1{jj}{ii}<0) = 0; %rectification
            
	end
end

    
   % normaliztion

E_single = zeros(size(im, 1), size(im,2),4);
for ii = 1:4
    for jj = 1:numChannel
        E_single(:,:,ii) = E_single(:,:,ii) + s1{jj}{ii}.^2 / numChannel;
    end
end

for ii = 1:4
    for jj = 1:numChannel
        s1{jj}{ii} = sqrt(k * s1{jj}{ii}.^2 ./ (sigma.^2 + E_single(:,:,ii)));
    end
end 