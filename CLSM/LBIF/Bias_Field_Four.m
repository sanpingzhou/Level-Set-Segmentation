function Bias = Bias_Field_Four (Img, evolvePhi1, evolvePhi2, K, c1, c2, c3, c4);

Region1 = double((evolvePhi1 <= 0)).*double((evolvePhi2 <= 0));
Region2 = double((evolvePhi1 <= 0)).*double((evolvePhi2 >= 0));
Region3 = double((evolvePhi1 >= 0)).*double((evolvePhi2 <= 0));
Region4 = double((evolvePhi1 >= 0)).*double((evolvePhi2 >= 0));


a1 = c1.*Region1;
b1 = c1.*c1.*Region1;

a2 = c2.*Region2;
b2 = c2.*c2.*Region2;

a3 = c3.*Region3;
b3 = c3.*c3.*Region3;

a4 = c4.*Region4;
b4 = c4.*c4.*Region4;

Bias = conv2(Img.*(a1 + a2 + a3 + a4), K, 'same')./conv2((b1+b2 + b3 + b4), K, 'same');

end

