 a biologically inspired grayscale/SO/DO GIST model for scene categorization
===============

GIST model was original proposed in:
Oliva, A., Torralba, A.: Modeling the shape of the scene : A holistic representation
of the spatial envelope.IJCV, 2001.

demoRelease  : grayscale Gist
demoSoRelease: SOGist (For most cases, two orientations of SO is sufficient due to the weekly oriented property)
demoDoRelease: DOGist



Filtering in the original algorithm was done in the frequency domain. We thus had to design
filters in the spatial domain that would approximate the original system as best as possible. 
In practice, we found that the Gabor lter parameters used in the
HMAX model with six scales (7x7 to 39x39 in steps of 6 pixels), two phases
(0 and 90 degrees), and eight orientations (0-180 degrees in steps of 22.5 degrees) led to results
comparable to those obtained with the grayscale gist model.


If you use it, please cite:
Zhang J., Barhomi Y., and Serre T. A new biologically inspired color image descriptor.In: ECCV, Florence, Italy, October 2012. 



For comments or questions, please contact Jun Zhang (zhangjun1126@gmail.com)

