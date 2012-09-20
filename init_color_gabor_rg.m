function [fSiz,filter, filters,c1OL,numSimpleFilters] = init_color_gabor(rot, RF_siz, Div,numChannel,numPhases)

filters = {};
filter = {};

fSiz_temp = repmat(RF_siz, length(rot)*numChannel,1);%17*24
% fSiz_temp = repmat(RF_siz, 1,length(rot)*numChannel)';%17*24
% fSiz_temp = repmat(RF_siz, length(rot),1);%17*24

fSiz = fSiz_temp(:);

rr = unique(RF_siz);

for i =1:length(rr)
    
%     [f, ff]   = get_filters_gabor(rr(i), rot, Div(i)); %filter:pos and neg gabor;
    [f, ff]   = get_filters_gabor_phases_rg(rr(i), rot, Div(i),numChannel,numPhases); %filter:pos and neg gabor;
    
    filters{i} = ff;
    filter{i} = f;
    
end

c1OL = 2;
numSimpleFilters = length(rot);




