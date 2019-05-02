function [C1, C2, C3, C4] = Average_Intensity_Four (Img, evolvePhi1, evolvePhi2, K, biasField);

Region1 = double((evolvePhi1 <= 0)).*double((evolvePhi2 <= 0));
Region2 = double((evolvePhi1 <= 0)).*double((evolvePhi2 >= 0));
Region3 = double((evolvePhi1 >= 0)).*double((evolvePhi2 <= 0));
Region4 = double((evolvePhi1 >= 0)).*double((evolvePhi2 >= 0));

a1 = conv2(biasField, K, 'same').*Img.*Region1;
b1 = conv2(biasField.^2, K, 'same').*Region1;
C1 = sum(a1(:))./(sum(b1(:)) + 1e-5);

a2 = conv2(biasField, K, 'same').*Img.*Region2;
b2 = conv2(biasField.^2, K, 'same').*Region2;
C2 = sum(a2(:))./(sum(b2(:)) + 1e-5);

a3 = conv2(biasField, K, 'same').*Img.*Region3;
b3 = conv2(biasField.^2, K, 'same').*Region3;
C3 = sum(a3(:))./(sum(b3(:)) + 1e-5);

a4 = conv2(biasField, K, 'same').*Img.*Region4;
b4 = conv2(biasField.^2, K, 'same').*Region4;
C4 = sum(a4(:))./(sum(b4(:)) + 1e-5);

end

