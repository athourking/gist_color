function filters = get_filters_gabor(or, n, numChannel)


Nscales = length(or);
Nfilters = sum(or);

l=0;
for i=1:Nscales
    for j=1:or(i)
        l=l+1;
        param(l,:)=[.35 .3/(1.85^(i-1)) 16*or(i)^2/32^2 pi/(or(i))*(j-1)];
    end
end

% Frequencies:
[fx, fy] = meshgrid(-n/2:n/2-1);
fr = fftshift(sqrt(fx.^2+fy.^2));
t = fftshift(angle(fx+sqrt(-1)*fy));


filter1 = createGabor_color(Nfilters,param, fr, t, n, 'positive');
filter2 = createGabor_color(Nfilters,param, fr, t, n, 'negative');


filter = zeros(n,n,2,Nfilters);
%seperate to pos and neg parts
for j=1:Nfilters
    filter(:,:,1,j) = abs(filter1(:,:,j));
    filter(:,:,2,j) = abs(filter2(:,:,j));
end


if numChannel == 8
    weights = [1/sqrt(2),-1/sqrt(2),0; 2/sqrt(6),-1/sqrt(6),-1/sqrt(6); 1/sqrt(6),1/sqrt(6),-2/sqrt(6); 1/sqrt(3),1/sqrt(3),1/sqrt(3)]';
else if numChannel == 6;
        weights = [1/sqrt(2),-1/sqrt(2),0; 2/sqrt(6),-1/sqrt(6),-1/sqrt(6); 1/sqrt(6),1/sqrt(6),-2/sqrt(6)]';
    end
end

% weights = [1,-1,0;1,-1,-1;1,1,-1;1,1,1]';
weights = [weights , -weights];

filters = zeros(n,n, 3, 8, Nfilters);

for i=1:3
    for j=1:numChannel
        if weights(i,j) > 0
            %take the positive gabor
            hh = 1;
        elseif weights(i,j) < 0
            %take the negative gabor
            hh = 2;
        elseif weights(i,j) == 0
            hh = 1;
        end
        
        filters(:,:,i,j ,:) =  weights(i,j) * filter(:,:,hh,:);
        
    end
end

% %filters_gabor = reshape(filters_all, size_f, size_f, 3, 32);
% for i=1:3
%     for j=1:numChannel
%         for k=1:4
%             nn = norm(filters(:,:,i,j,k),2);
%             if nn ~= 0
%                 filters(:,:,i,j,k) = filters(:,:,i,j,k)/nn;
%             end
%         end
%     end
% end


function fVals = createGabor_color(Nfilters,param, fr, t, n, gabor_sign)                                              

% Transfer functions:
fVals = zeros([n n Nfilters]);%ʱ��
% G = zeros([n n Nfilters]);%Ƶ��

for i = 1:Nfilters
    par = param(i,:);
    tr = t+param(i,4); 
    tr = tr+2*pi*(tr<-pi)-2*pi*(tr>pi);

    G = exp(-10*param(i,1)*(fr/n/param(i,2)-1).^2-2*param(i,3)*pi*tr.^2);
    
%     a = real(ifft2(G(:,:,i)));
    a = fftshift(real(ifft2(fftshift(G))));
    
    
    if strcmp(gabor_sign, 'positive')
        a(a<0) = 0;
    elseif  strcmp(gabor_sign, 'negative')
        a(a>0) = 0;
    end
    
    fVals(:,:,i) = a;
    
end


% if nargout == 0
%     figure
%     for i=1:Nfilters
%         max(max(G(:,:,i)))
%         contour(fftshift(G(:,:,i)),[1 .7 .6],'r');
%         hold on
%         drawnow
%     end
%     axis('on')
%     axis('square')
%     axis('ij')
% end
    
 

return;

