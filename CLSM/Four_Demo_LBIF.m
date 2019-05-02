close all; clear all; clc;

addpath(genpath('InputData'));
addpath(genpath('LBIF'));
addpath(genpath('Results'));

paraFour = Make_Default_Parameter_Four;
I = imread('9.jpg');
[m, n, h] = size (I);
if h == 3
    Img = double(rgb2gray(I));
else
    Img = double(I);
end
G = fspecial('gaussian', paraFour.fradius, paraFour.fsigma);
K = fspecial('gaussian', paraFour.kradius, paraFour.ksigma);

figure(1);
BW1 = roipoly(I);
BW2 = roipoly(I);
initialPhi1 =double(5*2*(BW1 - 0.5));
initialPhi2 =double(5*2*(BW2 - 0.5));
hold on;
[c1,h1] = contour(initialPhi1, [0, 0], 'r', 'Linewidth', 1.5);
[c2,h2] = contour(initialPhi2, [0, 0], 'g', 'Linewidth', 1.5);
title('the initial contour');
hold off;

evolvePhi1 = initialPhi1;
evolvePhi2 = initialPhi2;
biasField = double(ones(size(Img)));

for i = 1 : paraFour.nlter
    [evolvePhi1, evolvePhi2, biasField] = Evolve_LSF_Four(Img, evolvePhi1, evolvePhi2, G, K, biasField, paraFour );
    if mod(i, 10) == 0
        pause(0.1);
        figure(2);
        imshow(I);
        hold on;
        [c1,h1] = contour(evolvePhi1, [0, 0], 'r', 'Linewidth', 1.5);
        [c2,h2] = contour(evolvePhi2, [0, 0], 'g', 'Linewidth', 1.5);
        numIter = [num2str(i), 'iterations : '];
        title(numIter);
        hold off;
    end      
end