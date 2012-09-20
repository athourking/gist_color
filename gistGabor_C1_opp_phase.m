function g = gistGabor_C1_opp_phase(s1,numScaleBands,numChannel,numSimpleFilters,numPhases)
%
% Input:
%   img = input image (it can be a block: [nrows, ncols, c, Nimages])
%   w = number of windows (w*w)
%   G = precomputed transfer functions
%
% Output:
%   g: are the global features = [Nfeatures Nimages],
%                    Nfeatures = w*w*Nfilters*c

% if ndims(img)==2
%     c = 1;
%     N = 1;
% end
% if ndims(img)==3
%     [nrows ncols c] = size(img);
%     N = c;
% end
% if ndims(img)==4
%     [nrows ncols c N] = size(img);
%     img = reshape(img, [nrows ncols c*N]);
%     N = c*N;
% end
% c = 1;
N = numChannel;
w = 4;
Nfilters = numScaleBands * numSimpleFilters * numPhases;

% [n n Nfilters] = size(G);
W = w*w;
g = zeros([W*Nfilters N]);

% img = single(fft2(img));
k=0;


for iBand = 1:numScaleBands
%     for iScale = 1:numScales %�ӵ�һ���߶ȿ�ʼѭ��
        for iPhase = 1:numPhases
        for iFilt = 1:numSimpleFilters %�ӵ�һ������ʼѭ��
            
            ig = s1{iBand}(:,:,:,iFilt,iPhase);
            v = downN(ig, w);
            g(k+1:k+W,:) = reshape(v, [W N]);
            k = k + W;
            drawnow
        end
        end
%     end
end




% if c == 3
%     % If the input was a color image, then reshape 'g' so that one column
%     % is one images output:
%     g = reshape(g, [size(g,1)*3 size(g,2)/3]);
% end
g = reshape(g, [size(g,1)*N size(g,2)/N]);


function y=downN(x, N)
%
% averaging over non-overlapping square image blocks
%
% Input
%   x = [nrows ncols nchanels]
% Output
%   y = [N N nchanels]

nx = fix(linspace(0,size(x,1),N+1));
ny = fix(linspace(0,size(x,2),N+1));
y  = zeros(N, N, size(x,3));
for xx=1:N
    for yy=1:N
        v=mean(mean(x(nx(xx)+1:nx(xx+1), ny(yy)+1:ny(yy+1),:),1),2);
        y(xx,yy,:)=v(:);
    end
end
