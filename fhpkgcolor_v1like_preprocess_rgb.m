function im = fhpkgcolor_v1like_preprocess_rgb(imRaw, q)

	imrg = double(imRaw);
% 	% image resize bicubic, the biggest edge should 150
% 	if size(imr,3) == 1
% 		% input image is gray
% 		imrg = imr;
% 	else
% 		%input image is rgb
% 		imrg = rgb2gray(imr);
%     end

    for kk = 1:3
        tmp = imrg(:,:,kk);
	    tmp = double(tmp);
	    tmp = tmp - min(tmp(:));
	    imrg(:,:,kk) = tmp / max(tmp(:));
    end
    

	k = 1/3*ones(1,3);

	imga0 = convn(convn(imrg, k, 'same'), k', 'same');

	imga0 = imga0 - mean(imga0(:));
	if std(imga0(:)) ~=0
		imga0 = imga0/std(imga0(:));
	end
	kshape = [3,3];
	imga1 = v1like_norm(imga0, 'same', kshape, 1, q);
        im = imga1;
end



