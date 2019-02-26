function [phi1,phi2,By] = EVOL_MLBCV(Img, phi1_0,phi2_0, lambda1_1, lambda1_2,lambda2_1,lambda2_2, mu_1,mu_2, nu_1,nu_2, timestep, epsilon,G1,G2,By_0,numIter);
%   This function updates the level set function and the bias field according to the MLBCV model 
%   inputs: 
%       Img: input image
%       phi1_0:the level set function 1 to be updated
%       phi2_0:the level set function 2 to be updated
%       lambda1_1:  weight for c1 fitting term
%       lambda1_2:  weight for c2 fitting term
%       lambda2_1:  weight for c3 fitting term
%       lambda2_2:  weight for c4 fitting term
%       mu_1:the coefficient of the length penalizing Term 1
%       mu_2:the coefficient of the length penalizing Term 2
%       nu_1: the coefficient of the length term 1
%       nu_2: the coefficient of the length term 2
%       timestep: time step
%       epsilon: parameter for computing smooth Heaviside and dirac function
%       G1: Gaussian filter 1
%       G2: Gaussian filter 2
%       By_0: the bias field to be updated
%       numIter: number of iterations
%   outputs: 
%       phi1:updated level set function 1
%       phi2:updated level set function 2
%       By:updated bias field
%   created on 03/31/2013
%   Author: Sanping Zhou
%   email: csustzhousp@163.com


phi1=phi1_0;
phi2=phi2_0;
By=By_0;
for k=1:numIter
    %长度惩罚项计算
    phi1=NeumannBoundCond(phi1);
    [ux1,uy1]=gradient(phi1);
    Nx1=ux1./sqrt(ux1.^2 + uy1.^2 + 1e-10);
    Ny1=uy1./sqrt(ux1.^2 + uy1.^2 + 1e-10);
    [nxx1,junk1]=gradient(Nx1);
    [junk1,nyy1]=gradient(Ny1);
    K1=nxx1+nyy1;               
    penalizingTerm_1=mu_1*(4*del2(phi1)-K1);
   
    phi2= NeumannBoundCond(phi2);
    [ux2,uy2]=gradient(phi2);
    Nx2=ux2./sqrt(ux2.^2 + uy2.^2 + 1e-10);
    Ny2=uy2./sqrt(ux2.^2 + uy2.^2 + 1e-10);
    [nxx2,junk2]=gradient(Nx2);
    [junk2,nyy2]=gradient(Ny2);
    K2=nxx2+nyy2;               
    penalizingTerm_2=mu_2*(4*del2(phi2)-K2);
    
    %长度正则项计算
    dertax_1=Dirac(phi1,epsilon); 
    dertax_2=Dirac(phi2,epsilon); 
    [ux_1,uy_1]=gradient(phi1);                     %水平集函数的梯度
    [ux_2,uy_2]=gradient(phi2);                    
    Img_smooth=conv2(Img,G1,'same'); 
    [Ix,Iy]=gradient(Img_smooth);                   %图像的梯度
    f=Ix.^2+Iy.^2;
    g=1./(1+ f); 
    [vx,vy]=gradient(g);
    s_1=Ix.*ux_1+Iy.*uy_1;
    s_2=Ix.*ux_2+Iy.*uy_2;
    H_1=Heaviside_Revise(-s_1);
    H_2=Heaviside_Revise(-s_2);
    g_1=H_1./(1+ f);                                %演化曲线1的边缘指示函数
    g_2=H_2./(1+ f);                                %演化曲线2的边缘指示函数
    weightedLengthTerm_1=nu_1*dertax_1.*(vx.*Nx1 + vy.*Ny1 + g_1.*K1);
    weightedLengthTerm_2=nu_2*dertax_2.*(vx.*Nx2 + vy.*Ny2 + g_2.*K2);
    
    %拟合项计算
    [c1,c2,c3,c4] = Binaryfit_MLBCV( Img,G2,By,phi1,phi2); %计算各个区域的灰度均值
    By = By_LMBCV( Img,phi1,phi2,c1,c2,c3,c4,G2);          %计算偏移场
    [e1,e2]= Fitting_term(Img,phi1,phi2,c1,c2,c3,c4,lambda1_1,lambda1_2,lambda2_1,lambda2_2,G2,By);
    firstterm_1=dertax_1.*e1;
    firstterm_2=dertax_2.*e2;
  
    %迭代部分
    phi1=phi1+timestep*(firstterm_1+penalizingTerm_1+weightedLengthTerm_1); 
    phi2=phi2+timestep*(firstterm_2+penalizingTerm_2+weightedLengthTerm_2); 
 
end
end

