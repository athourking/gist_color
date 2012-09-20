% function [c2,s2,c1,s1] = C2(stim, filter, filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,s2Target,c1)
%function [c2,s2,c1,s1] = C2(stim,filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,s2Target,c1)
function c2_so = C2_so_mode(stim,filters,fSiz,c1SpaceSS,c1ScaleSS,c1OL,s2Target,numChannel,numPhases)
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


if nargin<10
%   [c1,s1] = C1_dopp_maxCol(stim, filter,filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,numChannel);
c1_so =  C1_soppNorm_maxCol(stim, filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,numChannel,numPhases);
end

nbands = length(c1_so);

c1BandImage_so = c1_so;

nfilts_so1 = size(c1_so{1},3);
nfilts_so2 = size(c1_so{1},4);
nfilts_so3 = size(c1_so{1},5);
n_rbf_centers = size(s2Target,2);

L_so = size(s2Target,1) / (nfilts_so1*nfilts_so2*nfilts_so3);

PatchSize_so = [L_so^.5,L_so^.5,nfilts_so1,nfilts_so2,nfilts_so3];

s2_so = cell(n_rbf_centers,1);
%Build s2:
%  for all prototypes in s2Target (RBF centers)
%   for all bands
%    calculate the image response
for iCenter = 1:n_rbf_centers

  Patch_so = reshape(s2Target(:,iCenter),PatchSize_so);
  s2_so{iCenter} = cell(nbands,1);
  for iBand = 1:nbands
     s2_so{iCenter}{iBand} = jun_WindowedPatchDistance_so(c1BandImage_so{iBand},Patch_so);
  end
end

%Build c2:
% calculate minimum distance (maximum stimulation) across position and scales

c2_so = inf(n_rbf_centers,1);
for iCenter = 1:n_rbf_centers
  for iBand = 1:nbands
     c2_so(iCenter) = min(c2_so(iCenter),min(min(s2_so{iCenter}{iBand}))); 
   
  end
end
