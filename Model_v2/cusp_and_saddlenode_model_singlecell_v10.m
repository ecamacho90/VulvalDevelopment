function dydt=cusp_and_saddlenode_model_singlecell_v10(t,y,p)


%This is the equation for a single cell without signals. It is used for the
%linear noise approximation to cumpute the initial condition for the
%system.


% Tau=1;    1    
% H=20;   2          
% M=2;     3     
% m11;  4
% m12;  5
% m21;  6
% m22;  7
% m31;  8
% m32;  9
% q1;   10
% q2;   11  
% q3;   12
% gam=0.05;   13
% H1;   14
% M1;   15
% H2;   16
% M2;   17
% s1=1;         18
% s2=1;         19
% s3=1;         20
% l1=1/2;       21
% l2=1/2;       22
% l3=1/2;       23
% alpha=1/9;    24
% lambdaE = 10; 25 %Exponential decay of EGF due to post competence
% lambdaN = 10; 26 %Exponential decay of Notch due to post competence
% xiE = 10;     27 %Exponential decay of EGF due to AC ablation
% xiN = 10;     28 %Exponential decay of Notch due to AC ablation
% p=[Tau,H,M,m11,m12,m21,m22,m31,m32,q1,q2,q3,gamma,H1,M1,H2,M2,s1,s2,s3,l1,l2,l3,alpha,lambdaE,lambdaN,xiE,xiN];

dydt=(p(1))  *  [((1+tanh(p(2)*(-y(2)+(p(3)*(sqrt(2)-1)/(2*sqrt(2))))))/2)         *              (-4*y(1)^3   -   2*p(10)*y(1)   -    p(11)   )           -          (1-((1+tanh(p(2)*(-y(2)+(p(3)*(sqrt(2)-1)/(2*sqrt(2))))))/2))*y(1);
    -y(2)*((y(2)-p(3))^2 + p(12))];