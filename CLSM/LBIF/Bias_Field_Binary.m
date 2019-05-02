function Bias = Bias_Field_Binary (Img, evolvePhi, K, c1, c2);


insideRegion = double(evolvePhi >= 0);
outsideRegion = 1 - insideRegion;

a1 = c1.*insideRegion;
b1 = c1.*c1.*insideRegion;

a2 = c2.*outsideRegion;
b2 = c2.*c2.*outsideRegion;

Bias = conv2(Img.*(a1 + a2), K, 'same')./conv2((b1+b2), K, 'same');
end

