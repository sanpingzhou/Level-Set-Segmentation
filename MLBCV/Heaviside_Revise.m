function H = Heaviside_Revise(A); 
%  H = Heaviside_Revise(A),if the of element A is less than zero,then the
%  corresponding element of H is zero,else the  corresponding element of H
%  is one
%inputs:
%         A: input matrix
%outputs:
%         H:output matrix
%created on 05/01/2013
%Authour:Sanping Zhou
%email:csustzhousp@163.com

 [m,n]=size(A);
 for i=1:m
     for j=1:n
         if A(i,j)<0
             H(i,j)=0;
         else
            H(i,j)=1;
         end
     end
 end

end

