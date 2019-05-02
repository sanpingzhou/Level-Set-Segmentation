function [C1, C2] = Average_Intensity_Binary (Img, evolvePhi, K, biasField);

insideRegion = double(evolvePhi >= 0);
outsideRegion = 1 - insideRegion;

a1 = conv2(biasField, K, 'same').*Img.*insideRegion;
b1 = conv2(biasField.^2, K, 'same').*insideRegion;
C1 = sum(a1(:))./(sum(b1(:)) + 1e-5);


a2 = conv2(biasField, K, 'same').*Img.*outsideRegion;
b2 = conv2(biasField.^2, K, 'same').*outsideRegion;
C2 = sum(a2(:))./(sum(b2(:)) + 1e-5);

end