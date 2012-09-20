function [filter, filters] = get_filters_gabor_phases(size_f, rot, div,numChannel,pCount)
%---------------------------gabor------------------------

lambda = size_f*2/div;%wavelength
sigma  = lambda.*0.8;%effective bandwidth
G      = 0.3; %spatial aspect ratio
% pCount = 2; %4 phases out of 90 degree

filter1 = GenerateGabor(size_f, rot, G, lambda, sigma, 'positive',pCount);
filter2 = GenerateGabor(size_f, rot, G, lambda, sigma, 'negative',pCount);

filter = {};
%seperate to pos and neg parts
for p = 1:pCount
    for j=1:length(rot)
        filter{p}(:,:,1,j) = abs(filter1{p}(:,:,j));
        filter{p}(:,:,2,j) = abs(filter2{p}(:,:,j));
    end
end

if numChannel == 8
    weights = [1/sqrt(2),-1/sqrt(2),0; 2/sqrt(6),-1/sqrt(6),-1/sqrt(6); 1/sqrt(6),1/sqrt(6),-2/sqrt(6); 1/sqrt(3),1/sqrt(3),1/sqrt(3)]';
else if numChannel == 6;
        weights = [2/sqrt(6),-1/sqrt(6),-1/sqrt(6);1/sqrt(6),1/sqrt(6),-2/sqrt(6); 1/sqrt(3),1/sqrt(3),1/sqrt(3)]';
    end
end

% weights = [1,-1,0;1,-1,-1;1,1,-1;1,1,1]';
weights = [weights , -weights];

% filters = zeros(size_f,size_f, 3, 8, length(rot));
filters = {};

for p = 1:pCount
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
            
            filters{p}(:,:,i,j ,:) =  weights(i,j) * filter{p}(:,:,hh,:);
        
        end
    end
end

%filters_gabor = reshape(filters_all, size_f, size_f, 3, 32);
for p = 1:pCount
    for i=1:3
        for j=1:numChannel
            for k=1:length(rot)
                nn = norm(filters{p}(:,:,i,j,k),2);
                if nn ~= 0
                    filters{p}(:,:,i,j,k) = filters{p}(:,:,i,j,k)/nn;
                end
            end
        end
    end    
end


function fVals = GenerateGabor(rfCount, rot, aspectRatio, lambda, sigma, gabor_sign,pCount)                                              

fVals = {};

points = (1 : rfCount) - ((1 + rfCount) / 2);


for p = 1:pCount
%     phase = (p - 1) / pCount * pi * 2;
    phase =  (p - 1) / 2 * pi;
    
    for f = 1 : length(rot)
    
%     theta = (f - 1) / fCount * pi;%0,45,90,135
    theta = rot(f) / 180 * pi;%0,45,90,135

    

    for j = 1 : rfCount
        for i = 1 : rfCount

            x = points(j) * cos(theta) - points(i) * sin(theta);
            y = points(j) * sin(theta) + points(i) * cos(theta);
       
        
                
                if sqrt(x * x + y * y) <= rfCount / 2
                    e = exp(-(x * x + aspectRatio * aspectRatio * y * y) / (2 * sigma * sigma));
                    e = e * cos(2 * pi * x / lambda + phase);
                else
                    e = 0;
                end
             fVals{p}(i, j, f) = e;   
        end    
    end
        
         %to make gabor mean 0, std 1
         a = fVals{p}(:,:,f);
         a = a - mean(a(:)); 
%          a = a / sqrt(sum(a(:) .* a(:)));
         a = a / std(a(:));
         fVals{p}(:,:,f) = a;      
         
    end
    
    if strcmp(gabor_sign, 'positive')
        fVals{p}(fVals{p}<0) = 0;
    elseif  strcmp(gabor_sign, 'negative')
        fVals{p}(fVals{p}>0) = 0;
    end   
    
end
    
 

return;

