function [ C1,C2,C3,C4] = Binaryfit_MLBCV( Img,G,B_field,phi1,phi2);
% [C1,C2,C3,C4]=binaryfit(inputs) computes c1,c2,c3,c4 for MLBCV model
%inputs:
%         Img: input image
%         B_field:the value of bias field
%         G:Gaussian filtering operator
%         phi1:the level set function 1
%         phi2:the level set function2
%outputs:
%         C1: a constant to fit the image U in the region phi1>0&&phi2>0
%         C2: a constant to fit the image U in the region phi1<0&&phi2>0
%         C3: a constant to fit the image U in the region phi1>0&&phi2<0
%         C4: a constant to fit the image U in the region phi1<0&&phi2<0
%created on 03/31/2013
%Authour:Sanping Zhou
%email:csustzhousp@163.com

%the partition of level set function
   Reg1_1=phi1>=0;
   Reg1_2=1-Reg1_1;
   Reg2_1=phi2>=0;
   Reg2_2=1-Reg2_1;  
       
%computes c1,c2,c3,c4
  [m,n,h]=size(Img);
  one=ones(m,n);
  aa1=conv2(B_field,G,'same');
  aa2=conv2(one,G,'same');
  part_1=Reg1_1.*Reg2_2;
  c1=sum(part_1(:));
  if c1==0
     C1=0;
  else
  a_1=part_1.*aa2.*Img;
  b_1=part_1.*aa1;
  c_1=a_1-b_1;
  d_1=part_1.*aa2;
  number_1=sum(c_1(:));
  denom_1=sum(d_1(:));
  C1=number_1/denom_1;
  end
  

  part_2=Reg1_1.*Reg2_1;
  c2=sum(part_2(:));
  if c2==0
     C2=0;
  else
  a_2=part_2.*aa2.*Img;
  b_2=part_2.*aa1;
  c_2=a_2-b_2;
  d_2=part_2.*aa2;
  number_2=sum(c_2(:));
  denom_2=sum(d_2(:));
  C2=number_2/denom_2;
  end

  
  part_3=Reg1_2.*Reg2_1;
  c3=sum(part_3(:));
  if c3==0
     C3=0;
  else
  a_3=part_3.*aa2.*Img;
  b_3=part_3.*aa1;
  c_3=a_3-b_3;
  d_3=part_3.*aa2;
  number_3=sum(c_3(:));
  denom_3=sum(d_3(:));
  C3=number_3/denom_3;
  end
  

  part_4=Reg1_2.*Reg2_2;
  c4=sum(part_4(:));
  if c4==0
     C4=0;
  else
  a_4=part_4.*aa2.*Img;
  b_4=part_4.*aa1;
  c_4=a_4-b_4;
  d_4=part_4.*aa2;
  number_4=sum(c_4(:));
  denom_4=sum(d_4(:));
  C4=number_4/denom_4;
  end
end

