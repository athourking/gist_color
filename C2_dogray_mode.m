% function [c2,s2,c1,s1] = C2(stim, filter, filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,s2Target,c1)
%function [c2,s2,c1,s1] = C2(stim,filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,s2Target,c1)
function [c2,c2_gray] = C2_dogray_mode(stim_gray,filters_gray,fSiz_gray,s2Target_gray,stim,filter,filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,s2Target,numChannel,numPhases)
% given an image extracts layers s1 c1 s2 and finally c2
% for inputs stim, filters, fSiz, c1SpaceSS,c1ScaleeSS, and c1OL
% see the documentation for C1 (C1.m)
%
% briefly, 
% stim is the input image. 
% filters fSiz, c1SpaceSS, c1ScaleSS, c1OL are the parameters of
% the c1 process
%
% s2Target are the prototype (patches) to be used in the extraction
% of s2.  Each patch of size [n,n,d] is stored as a column in s2Target,
% which has itself a size of [n*n*d, n_patches];
%
% if available, a precomputed c1 layer can be used to save computation
% time.  The proper format is the output of C1.m
%
% See also C1


if nargin<15
%   [c1,s1] = C1_dopp_maxCol(stim, filter,filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,numChannel);
[c1,c1_gray] = C1_dogray_mode(stim_gray,filters_gray,fSiz_gray,stim, filter,filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,numChannel,numPhases);
end

nbands = length(c1);
c1BandImage = c1;
c1BandImage_gray = c1_gray;

nfilts = size(c1{1},3);
nfilts_gray = size(c1_gray{1},3);
n_rbf_centers = size(s2Target,2);
L = size(s2Target,1) / (nfilts*size(c1{1},4));
L_gray = size(s2Target_gray,1) / nfilts_gray;
PatchSize = [L^.5,L^.5,nfilts,size(c1{1},4)];
PatchSize_gray = [L_gray^.5,L_gray^.5,nfilts_gray];


s2 = cell(n_rbf_centers,1);
s2_gray = cell(n_rbf_centers,1);

%Build s2:
%  for all prototypes in s2Target (RBF centers)
%   for all bands
%    calculate the image response
for iCenter = 1:n_rbf_centers
  Patch = reshape(s2Target(:,iCenter),PatchSize);
  s2{iCenter} = cell(nbands,1);
  Patch_gray = reshape(s2Target_gray(:,iCenter),PatchSize_gray);
  s2_gray{iCenter} = cell(nbands,1);
  for iBand = 1:nbands
     s2{iCenter}{iBand} = jun_WindowedPatchDistance(c1BandImage{iBand},Patch);  
     s2_gray{iCenter}{iBand} = WindowedPatchDistance(c1BandImage_gray{iBand},Patch_gray);  %����ԣ�ģ��ƥ�䷽��
  end
end

%Build c2:
% calculate minimum distance (maximum stimulation) across position and scales
c2 = inf(n_rbf_centers,1);
c2_gray = inf(n_rbf_centers,1);
for iCenter = 1:n_rbf_centers
  for iBand = 1:nbands
     c2(iCenter) = min(c2(iCenter),min(min(s2{iCenter}{iBand})));
     c2_gray(iCenter) = min(c2_gray(iCenter),min(min(s2_gray{iCenter}{iBand})));    
  end
end
