function s1 = single_opponents_S1(im, filters_gabor,numChannel,numSimpleFilters)

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

% normaliztion

E_single = zeros(size(im, 1), size(im,2),numSimpleFilters);
for ii = 1:numSimpleFilters
    for jj = 1:numChannel
        E_single(:,:,ii) = E_single(:,:,ii) + s1(:,:,jj,ii) .^2 / numChannel;
    end
end

for ii = 1:numSimpleFilters
    for jj = 1:numChannel
        s1(:,:,jj,ii)  = sqrt(k * s1(:,:,jj,ii) .^2 ./ (sigma.^2 + E_single(:,:,ii)));
    end
end

