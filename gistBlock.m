function g = gistBlock(c1,numScaleBands,numSimpleFilters)
% 


N = 1;
w = 4;
Nfilters = numScaleBands*numSimpleFilters;


W = w*w;
g = zeros([W*Nfilters N]);


k=0;
for iBand = 1:numScaleBands
    for iFilt = 1:numSimpleFilters 
        ig = c1{iBand}(:,:,iFilt);
        v = downN(ig, w);
        g(k+1:k+W,:) = reshape(v, [W N]);
        k = k + W;
        drawnow
    end
end




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
