function dydt=singlecell_cusp_and_saddlenode_model(t,y,p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%       Differential equation corresponding to the model when the 
%       parameters are fixed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%p(1)=c
%p(2)=b
%p(3)=a
%p(4)=H
%p(5)=M

Fy = (1 + tanh(p(4)*(-y(2) + (sqrt(2) - 1).*p(5)/(2*sqrt(2)))))/2;

dydt=[Fy.*(-4.*y(1).^3 - 2.*p(3).*y(1) - p(2)) - (1 - Fy).*y(1); -y(2).*((y(2) - p(5)).^2 + p(1))];