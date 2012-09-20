function output = gaussWh(img)

D = 16;
[x y] = meshgrid(-D/2:D/2-1);
G = exp(-0.5 * ((x.^2 + y.^2) / (D/2).^2));
G = G / sum(G(:));
imv = conv2padded (img.^2,G);
output = img ./ sqrt(imv + ~imv);

return