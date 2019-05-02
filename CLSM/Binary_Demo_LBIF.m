close all; clear all; clc;

addpath(genpath('InputData'));
addpath(genpath('LBIF'));
addpath(genpath('Results'));

paraBinary = Make_Default_Parameter_Binary;
I = imread('6.jpg');
[m, n, h] = size (I);
if h == 3
    Img = double(rgb2gray(I));
else
    Img = double(I);
end

G = fspecial('gaussian', paraBinary.fradius, paraBinary.fsigma);
K = fspecial('gaussian', paraBinary.kradius, paraBinary.ksigma);
figure(1);
BW = roipoly(I);
initialPhi = double(20*2*(BW - 0.5));
hold on;
[c,h] = contour(initialPhi, [0, 0], 'g', 'Linewidth', 1.5);
title('the initial contour');
hold off;


evolvePhi = initialPhi;
biasField = ones(size(Img));

for i = 1 : paraBinary.nlter
    [evolvePhi, biasField] = Evolve_LSF_Binary(Img, evolvePhi, G, K, biasField, paraBinary );
    if mod(i, 5) == 0
        pause(0.1);
        figure(2);
        imshow(I);
        hold on;
        [c,h] = contour(evolvePhi, [0, 0], 'r', 'Linewidth', 1.5);
        numIter = [num2str(i), 'iterations : '];
        title(numIter);
        hold off;
    end      
end