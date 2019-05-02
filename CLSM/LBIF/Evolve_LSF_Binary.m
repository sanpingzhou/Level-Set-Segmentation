function [Phi, Bias] = Evolve_LSF_Binary(Img, evolvePhi, G, K, biasField, paraBinary );



evolvePhi = Neumann_Bound_Cond(evolvePhi);
[ux,uy] = gradient(evolvePhi);
Nx = ux./sqrt(ux.^2 + uy.^2 + 1e-10);
Ny = uy./sqrt(ux.^2 + uy.^2 + 1e-10);
[nxx, ~] = gradient(Nx);
[~, nyy] = gradient(Ny);
curvature = nxx + nyy;
distancePenaltyTerm = paraBinary.mu*(4*del2(evolvePhi) - curvature);


dertax = Dirac(evolvePhi, paraBinary.efso);
smoothImg = conv2(Img, G, 'same');
[Ix, Iy] = gradient(smoothImg);
g = 1./(1 + Ix.^2 + Iy.^2);
[vx, vy] = gradient(g);
lengthPenaltyTerm = paraBinary.nu*dertax.*(vx.*Nx + vy.*Ny + curvature);


insideRegion = double(evolvePhi >= 0);
outsideRegion = 1 - insideRegion;
[c1, c2] = Average_Intensity_Binary (Img, evolvePhi, K, biasField);
Bias = Bias_Field_Binary (Img, evolvePhi, K, c1, c2);
fittedImg = (c1.*insideRegion + c2.*outsideRegion).*Bias;
diffImg = fittedImg - Img;
correnDiffImg = exp(diffImg./paraBinary.delte);
diff  = conv2(diffImg.*correnDiffImg, K, 'same');
fittingEnergyTerm = -paraBinary.lamda*(c1 - c2).*dertax.*Bias.*diff;


Phi = evolvePhi + paraBinary.timestep*(distancePenaltyTerm + lengthPenaltyTerm +  fittingEnergyTerm);

end


