function filters = init_color_realgabor(or, n, numChannel)

Nscales = length(or);
Nfilters = sum(or);

% l=0;
for i=1:Nscales
    for j=1:or(i)
%         l=l+1;
        param{i}(j,:) = [.35 .3/(1.85^(i-1)) 16*or(i)^2/32^2 pi/(or(i))*(j-1)];
    end
end




filters = {};
filter = {};

% fSiz_temp = repmat(RF_siz, length(rot)*numChannel,1);%17*24
% % fSiz_temp = repmat(RF_siz, length(rot),1);%17*24
% 
% fSiz = fSiz_temp(:);

% rr = unique(RF_siz);

for i =1:Nscales
    
%     [f, ff]   = get_filters_gabor(rr(i), rot, Div(i)); %filter:pos and neg gabor;
    [f, ff]   = get_filters_realgabor(param{i}, or,n,numChannel); %filter:pos and neg gabor;
    
    filters{i} = ff;
    filter{i} = f;
    
end

c1OL = 2;
numSimpleFilters = or(1);


