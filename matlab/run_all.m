% run all
clear all;
close all;
load('../matlab_support/bounds.mat')
load('../matlab_support/Wpts.mat')
I1 = imread('../images/target_01.png');
I2 = imread('../images/target_02.png');
I3 = imread('../images/target_03.png');
I4 = imread('../images/target_04.png');
I5 = imread('../images/target_05.png');

corners1 = cross_junctions(I1, bpoly1, Wpts);
corners2 = cross_junctions(I2, bpoly2, Wpts);
corners3 = cross_junctions(I3, bpoly3, Wpts);
corners4 = cross_junctions(I4, bpoly4, Wpts);
corners5 = cross_junctions(I5, bpoly5, Wpts);
