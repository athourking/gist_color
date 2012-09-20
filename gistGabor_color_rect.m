function g = gistGabor_color_rect(img, w, filters)
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
% % 
% % [nrows ncols c N] = size(img);
% % img = reshape(img, [nrows ncols c*N]);
% % N = c*N;

[n n c N Nfilters] = size(filters);
W = w*w;
g = zeros([W*Nfilters N]);

% img = single(fft2(img)); 
k=0;
% s = zeros(size(img,1),size(img,2),N,Nfilters);
% ig = 0;

for n = 1:Nfilters
%     n
    s = zeros(size(img,1),size(img,2),N);
    
    for jj = 1:N
        
        for kk = 1:3
            %             ig = abs(ifft2(img.*repmat(G(:,:,n), [1 1 N])));
            tmp = conv2padded(img(:,:,kk), squeeze(filters(:,:,kk,jj,n)));
            s(:,:,jj) = s(:,:,jj) + tmp;%e.g. +R-center,-G-surround...
        end
        
    end
%     s = abs(s);%8*32
    s(s<0) = 0;
%     ig = s(:,:,:,n);%256*256*8
    v = downN(s, w);
    g(k+1:k+W,:) = reshape(v, [W N]);
    k = k + W;
    drawnow
end

g = reshape(g, [size(g,1)*N size(g,2)/N]);


% if c == 3
%     % If the input was a color image, then reshape 'g' so that one column
%     % is one images output:
%     g = reshape(g, [size(g,1)*3 size(g,2)/3]);
% end


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
