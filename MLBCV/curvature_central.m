function k = curvature_central(phi)
% k = curvature_central(phi)compute curvature
%input:
%       phi:the level set function
%output:
%       k:the value of curvature
%Authour:Li chunming
%http://www.mathworks.com/matlabcentral/fileexchange/12711-level-set-for-im
%      age-segmentation
[ux, uy] = gradient(phi);
normDu = sqrt(ux.^2+uy.^2+1e-10);
Nx = ux./normDu; Ny = uy./normDu;
[nxx, junk] = gradient(Nx);
[junk, nyy] = gradient(Ny);
k = nxx+nyy;         
end

