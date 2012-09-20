function cI = jun_readAllImages_8scenes(train,test,datapath)
%Reads all training and testing images into a cell of length 4
%cI{1} = train_set.pos,
%cI{2} = train_set.neg,
%cI{3} = test_set.pos,
%cI{4} = test_set.neg,
if nargin<2
  maximagesperdir = inf;
end

files = dir(fullfile(datapath, '*jpg'));
train_files = files(train);
test_files = files(test);

for i = 1 : numel(train)
    train_set{i} = fullfile(datapath,train_files(i).name);
end
for i = 1 : numel(test)
    test_set{i} = fullfile(datapath,test_files(i).name);
end

% fprintf('Reading images...');    
% cI = cell(length(trainName),1);
% for i = 1:length(trainName),
%     cI{i} = double(imread(fullfile('/users/jz7/data/ColorGist/images/spatial_envelope_256x256_static_8outdoorcategories',sprintf(trainName{i}))));  
% end
dnames = {train_set,test_set};

fprintf('Reading images...');    
cI = cell(2,1);
for i = 1:2,
  cI{i} = cell(length(dnames{i}),1);
  for j = 1:length(dnames{i}),
    cI{i}{j} = double(imread(dnames{i}{j}))./255;  
  end
end

fprintf('done.\n');
