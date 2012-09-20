% Demo, scene recognition
%
% This script trains a SVM for classification.
%
% You can use any SVM code. 
%
% Here, the code assumes you have installed the next toolbox:
% Support Vector Machine toolbox
% Version 2.51, January 2002
% available at: http://ida.first.fraunhofer.de/~anton/software.html

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
datapath = '/users/jz7/data/database/indoorCVPR09';
train = importdata('/users/jz7/data/database/indoorCVPR09/TrainImages.txt');
test  = importdata('/users/jz7/data/database/indoorCVPR09/TestImages.txt');
cI = cell(2,1);
cI{1} = train;
cI{2} = test;
clear train test

numberBlocks = 4;
fc_prefilt = 4;


orientationsPerScale = [8 8 8 8];
numberBlocks = 4;
fc_prefilt = 4;
NtrainingPerClass = 80;
imageSize = 200;
% Nclasses = length(categories);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute global features
% scenes = dir(fullfile(HOMEIMAGES, '*.jpg'));
% scenes = {scenes(:).name};
% Nscenes = length(scenes);
% 
% % Precompute filter transfert functions (only need to do this once, unless
% % image size is changes):
% % createGabor(orientationsPerScale, imageSize);
G = createGabor(orientationsPerScale, imageSize);
Nfeatures = size(G,3)*numberBlocks^2*3;
% % 
% Loop: Compute global features for all scenes
% F = zeros([Nscenes Nfeatures]);
% C = zeros([Nscenes 1]);
outDir = sprintf('/gpfs/data/jz7/ColorGist/results/indoor/%i_train/C1',NtrainingPerClass);
if ~exist(outDir)
    mkdir(outDir);
end
% %  height = 200;
F = cell(2,1);

for ii = 1:2
    F{ii} = zeros([numel(cI{ii}) Nfeatures]);
    for n = 1:numel(cI{ii})
n

    img = imread(fullfile(datapath, cI{ii}{n}));
    
[mm,nn,unused] = size(img);
     if mm ~= imageSize || nn ~= imageSize
         img = imresize(img, [imageSize imageSize], 'bilinear');
%         img = imresize(img,[imageSize NaN]);
    end
    if unused == 1
        im(:,:,1) = img;
        im(:,:,2) = img;
        im(:,:,3) = img;
        clear img
        img = im;
        clear im;
    end
   
    output = prefilt(double(img), fc_prefilt);
    g = gistGabor(output, numberBlocks, G);
    
    F{ii}(n,:) = g;
    end
end

save(fullfile(outDir,sprintf('F_rgb_indoor.mat')) ,'F','-v7.3');




% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Split training/test (train = index training samples, test = index test)
% train = [];
% for c = 1:Nclasses
%     j = find(C==c);%�ҳ���c�ೡ����ͼ����
%     t = randperm(length(j));%��j��ͼ���������
%     train = [train; j(t(1:NtrainingPerClass))];%ȡǰ100��ͼ����Ϊtraining set
% end %�õ�800��ͼ���training set��ÿ��100��
% test = setdiff(1:Nscenes, train)';%
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Train SVM
% % train one versus all and then take maximal score
% 
% % fig = figure;
% for c = 1:Nclasses
%     netc = svm(Nfeatures, 'rbf', 0.003, 100);
%     netc = svmtrain(netc, F(train,:), 2*(C(train)==c)-1, [], 1);
% 
%     [Y, scores(c,:)] = svmfwd(netc, F(test,:)); %��c���Ӧ�Ĳ��Լ�Ԥ�����ͷ���
%     
% end
% for k = 1:length(test)
%     [foo, ctest_hat(k)] = max(scores(:,k));
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Plot performance
% % Confusion matrix:
% Nclasses = length(categories);
% Cm = zeros(Nclasses, Nclasses);
% for j = 1:Nclasses
%     for i = 1:Nclasses
%         % row i, col j is the percentage of images from class i that
%         % were missclassified as class j.
%         Cm(i,j) = 100*sum((C(test)==i) .* (ctest_hat'==j))/sum(C(test)==i);
%     end
% end
% 
% perf = mean(diag(Cm))%NP
% 
% 
% 
% save gist_torralba
% 
% 
exit
% % 
% % 
% % figure
% % subplot(121)
% % imagesc(Cm); axis('square'); colorbar
% % subplot(122)
% % bar(diag(Cm))%���Խ����ϵ�ֵ
% % title(mean(diag(Cm)))%8������ʶ���ʾ�ֵ
% % axis('square')
