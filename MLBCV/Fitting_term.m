function  [E1,E2]= Fitting_term( Img,phi1,phi2,C1,C2,C3,C4,lambda1_1,lambda1_2,lambda2_1,lambda2_2,G,By);
% [E1,E2]= Fitting_term( Img,phi1,phi2,C1,C2,C3,C4,G) computes e1,e2
% for BBCV_M model         
%inputs:
%         Img: input image
%         phi1:the level function 1
%         phi1:the level function 2
%         C1: a constant to fit the image U in the region phi1>0&&phi2>0
%         C2: a constant to fit the image U in the region phi1<0&&phi2>0
%         C3: a constant to fit the image U in the region phi1>0&&phi2<0
%         C4: a constant to fit the image U in the region phi1<0&&phi2<0
%         lambda1_1:the coefficient of region 1
%         lambda1_2:the coefficient of region 2
%         lambda2_1:the coefficient of region 3
%         lambda2_2:the coefficient of region 4
%         G:Gaussian filtering operator
%         By:the value of bias field 
%outputs:
%         E1: the value of fitting term 1
%         E2: the value of fitting term 2
%created on 03/31/2013
%Authour:Sanping Zhou,all right reserved
%email:csustzhousp@163.com

%the partition of level set function
   Reg1_1=phi1>=0;
   Reg1_2=1-Reg1_1;
   Reg2_1=phi2>=0;
   Reg2_2=1-Reg2_1;
   
   [m,n,h]=size(Img);
   one=ones(m,n);
   aa=conv2(one,G,'same');
   by=By.^2;                                             %为了增加计算的速度
   By2=conv2(by,G,'same');
   aa1=conv2(By,G,'same');
   c1_1=Img-C1;
   c1_2=c1_1.^2.*aa;
   c1_3=2*c1_1.*aa1;
   
   c2_1=Img-C2;
   c2_2=c2_1.^2.*aa;
   c2_3=2*c2_1.*aa1;
   
   c3_1=Img-C3;
   c3_2=c3_1.^2.*aa;
   c3_3=2*c3_1.*aa1;
   
   c4_1=Img-C4;
   c4_2=c4_1.^2.*aa;
   c4_3=2*c4_1.*aa1;
   
   E1=-lambda1_1*(c1_2-c1_3+By2).*Reg2_2...
      -lambda1_1*(c2_2-c2_3+By2).*Reg2_1...
      +lambda1_2*(c3_2-c3_3+By2).*Reg2_1...
      +lambda1_2*(c4_2-c4_3+By2).*Reg2_2;
  
   E2= lambda2_2*(c1_2-c1_3+By2).*Reg1_1...
      -lambda2_1*(c2_2-c2_3+By2).*Reg1_1...
      -lambda2_1*(c3_2-c3_3+By2).*Reg1_2...
      +lambda2_2*(c4_2-c4_3+By2).*Reg1_2;
   
end

