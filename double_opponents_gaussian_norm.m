function s1 = double_opponents_gaussian_norm(im, filters_gaussian,numChannel)


E_type2 = 0;
sigma = 0.125;
k = 1;



s1 = zeros(size(im,1), size(im,2),numChannel);

% for ii = 1:4
for jj = 1:numChannel
        
        for kk=1:3
            tmp = conv2padded(im(:,:,kk), squeeze(filters_gaussian(:,:,kk, jj)));
            s1(:,:,jj) = s1(:,:,jj) + tmp;%+R-G,no center-surround structure
        end
%     end
%
     s1(s1<0) = 0;
%     s1{jj}(s1{jj}<0) = 0;%rectify
%     
    E_type2 = E_type2 + s1(:,:,jj).^2 /numChannel;%normalization
    
end




for jj = 1:numChannel
    s1(:,:,jj) = sqrt(k * s1(:,:,jj).^2 ./ (sigma.^2 + E_type2));
end


return