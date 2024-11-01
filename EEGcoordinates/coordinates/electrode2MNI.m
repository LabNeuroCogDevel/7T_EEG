function [MNI] = electrode2MNI(x,y,z,Affine)


V=[x,y,z,1];
MNI=inv(Affine)\V';
