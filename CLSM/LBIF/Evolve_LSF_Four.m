function [Phi1, Phi2, Bias] = Evolve_LSF_Four(Img, evolvePhi1, evolvePhi2, G, K, biasField, paraFour);

evolvePhi1 = Neumann_Bound_Cond(evolvePhi1);
[ux1,uy1] = gradient(evolvePhi1);
Nx1 = ux1./sqrt(ux1.^2 + uy1.^2 + 1e-10);
Ny1 = uy1./sqrt(ux1.^2 + uy1.^2 + 1e-10);
[nxx1, ~] = gradient(Nx1);
[~, nyy1] = gradient(Ny1);
curvature1 = nxx1 + nyy1;
distancePenaltyTerm1 = paraFour.mu1*(4*del2(evolvePhi1) - curvature1);

evolvePhi2 = Neumann_Bound_Cond(evolvePhi2);
[ux2,uy2] = gradient(evolvePhi2);
Nx2 = ux2./sqrt(ux2.^2 + uy2.^2 + 1e-10);
Ny2 = uy2./sqrt(ux2.^2 + uy2.^2 + 1e-10);
[nxx2, ~] = gradient(Nx2);
[~, nyy2] = gradient(Ny2);
curvature2 = nxx2 + nyy2;
distancePenaltyTerm2 = paraFour.mu2*(4*del2(evolvePhi2) - curvature2);

smoothImg = conv2(Img, G, 'same');
[Ix, Iy] = gradient(smoothImg);
g = 1./(1 + Ix.^2 + Iy.^2);
[vx, vy] = gradient(g);
dertax1 = Dirac(evolvePhi1, paraFour.efso);
lengthPenaltyTerm1 = paraFour.nu1*dertax1.*(vx.*Nx1 + vy.*Ny1 + curvature1);

dertax2 = Dirac(evolvePhi2, paraFour.efso);
lengthPenaltyTerm2 = paraFour.nu2*dertax2.*(vx.*Nx2 + vy.*Ny2 + curvature2);

Region1 = double((evolvePhi1 <= 0)).*double((evolvePhi2 <= 0));
Region2 = double((evolvePhi1 <= 0)).*double((evolvePhi2 >= 0));
Region3 = double((evolvePhi1 >= 0)).*double((evolvePhi2 <= 0));
Region4 = double((evolvePhi1 >= 0)).*double((evolvePhi2 >= 0));

[c1, c2, c3, c4] = Average_Intensity_Four (Img, evolvePhi1, evolvePhi2, K, biasField);
Bias = Bias_Field_Four(Img, evolvePhi1, evolvePhi2, K, c1, c2, c3, c4);
fittedImg = (c1.*Region1 + c2.*Region2 + c3.*Region3 + c4.*Region4).*Bias;
diffImg = fittedImg - Img;
correnDiffImg = exp(diffImg./paraFour.delte);
diff  = conv2(diffImg.*correnDiffImg, K, 'same');

fittingEnergyTerm1 = paraFour.lamda1*((c1-c3).*double(evolvePhi2 <= 0) + (c2 - c4).*double(evolvePhi2 >= 0)).*dertax1.*Bias.*diff;
fittingEnergyTerm2 = paraFour.lamda1*((c1-c2).*double(evolvePhi1 <= 0) + (c3 - c4).*double(evolvePhi1 >= 0)).*dertax2.*Bias.*diff;

Phi1 = evolvePhi1 + paraFour.timestep*(distancePenaltyTerm1 + lengthPenaltyTerm1 +  fittingEnergyTerm1);
Phi2 = evolvePhi2 + paraFour.timestep*(distancePenaltyTerm2 + lengthPenaltyTerm2 +  fittingEnergyTerm2);

end

