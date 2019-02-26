%the main program for BBCV model
%created on 05/01/2013
%Authour:Sanping Zhou
%email:csustzhousp@163.com
clear; clc;
close all;

%ˮƽ����ʼ��
I1=imread('14.jpg');
% I1=0.5*I1;
% I1 = imnoise(I1,'gaussian',0.15);    %��Ӹ�˹����
I=I1(:,:,1);                           %ѡȡһ��ͨ������
figure(1);
imshow(I1);
bw1=roipoly(I);
bw2=roipoly(I);
phi1=9*2*(bw1-0.5);
phi2=9*2*(bw2-0.5);
hold on;
imshow(I1);
[c,h]=contour(phi1,[0,0],'g','Linewidth',1.5);
[c,h]=contour(phi2,[0,0],'r','Linewidth',1.5);
hold off;


%��������
timestep=0.05;                          %ʱ�䲽��
lambda1_1=0.7;       lambda2_1=0.7;     %ˮƽ������1��2�У�ÿһ�����Ȩ��ϵ��
lambda1_2=0.5;       lambda2_2=0.5;
mu_1=0.2/timestep; mu_2=0.2/timestep;   %���ȳͷ���1��2ϵ��
nu_1=100;            nu_2=100;              %������ϵ��
efso=1;                                 %dertax��������
sigma=2;
numIter = 300;                          %��������

I=double(I);                            %��ͼ��ת����˫������
By=0.5*I;                               %����ƫ�Ƴ�����
By(:,:,:)=0;
G1=fspecial('gaussian',3,sigma);        %��˹�˲���Ϊƽ��ԭͼ��
G2=fspecial('gaussian',13,3);           %��˹�˲���Ϊƫ�Ƴ�

for n=1:numIter
   [phi1,phi2,By] = EVOL_MLBCV(I,phi1,phi2,lambda1_1,lambda1_2,lambda2_1,lambda2_2,mu_1,mu_2,nu_1,nu_2,timestep,efso,G1,G2,By,1);
   if mod(n,10)==0
      pause(0.1);
      figure(2);
      imshow(I1);
      hold on;
      [c,h]=contour(phi1,[0,0],'g','Linewidth',1.5);
      [c,h]=contour(phi2,[0,0],'r','Linewidth',1.5);
      iterNum=[num2str(n), ' iterations of real time'];        
      title(iterNum);       
      hold off;
   end
end

%ƫ�Ƴ�����Ч����֤
bias_b=zeros(1,256);
bias_l=zeros(1,256);
[I_m,I_n]=size(I); 
I_corrected=uint8(double(I)-By);
for i=1:I_m
    for j=1:I_n
        bias_b(I(i,j)+1)=bias_b(I(i,j)+1)+1;
    end
end
for i=1:I_m
    for j=1:I_n
        bias_l(I_corrected(i,j)+1)=bias_l(I_corrected(i,j)+1)+1;
    end
end
figure(3)
plot(bias_b,'r','Linewidth',1.5);
xlim([0 255]);
hold on
plot(bias_l,'b','Linewidth',1.5);
xlim([0 255]);
hold off
legend('Histogram Before Bias Correction','Histogram After Bias Correction');



