function train_set = generate_data_8scenes(datapath,numTrain)

files = dir(fullfile(datapath, '*jpg'));
% files = setdiff(files, {'.', '..'});

train_files = [];
% train_label = [];
test_files  = [];
% test_label  = [];

test_set = {};
% dirs = dirs;


train_files = files(inds(1:numTrain));


for i = 1 : numel(train_files)
    train_set{i} = fullfile(datapath,train_files(i).name);
end

return