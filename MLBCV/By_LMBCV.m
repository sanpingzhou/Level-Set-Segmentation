function By = By_LMBCV( Img,phi1,phi2,C1,C2,C3,C4,G);
% By = Bias_field( Img,phi1,phi2,C1,C2,C3,C4,G)computes bias field for MLBCV model         
%inputs:
%         Img: input image
%         phi1:the level function 1
%         phi1:the level function 2
%         C1: a constant to fit the image U in the region phi1>0&&phi2>0
%         C2: a constant to fit the image U in the region phi1<0&&phi2>0
%         C3: a constant to fit the image U in the region phi1>0&&phi2<0
%         C4: a constant to fit the image U in the region phi1<0&&phi2<0
%         G:Gaussian filtering operator
%         epsilon:parameter for function Heaviside
%outputs:
%         By: the value of bias field 
%created on 03/31/2013
%Authour:Sanping Zhou,all right reserved
%email:csustzhousp@163.com

%the partition of level set function
   Reg1_1=phi1>=0;
   Reg1_2=1-Reg1_1;
   Reg2_1=phi2>=0;
   Reg2_2=1-Reg2_1;
  
  [m,n,h]=size(Img);
  aa1=(Img-C1).*Reg1_1.*Reg2_2;  
  bb1=conv2(aa1,G,'same');
  aa2=(Img-C2).*Reg1_1.*Reg2_1; 
  bb2=conv2(aa2,G,'same');
  aa3=(Img-C3).*Reg1_2.*Reg2_1; 
  bb3=conv2(aa3,G,'same');
  aa4=(Img-C4).*Reg1_2.*Reg2_2; 
  bb4=conv2(aa4,G,'same');
  number=bb1+bb2+bb3+bb4;
  one=ones(m,n);
  denom=conv2(one,G,'same');
  By=number./denom;
 

end

