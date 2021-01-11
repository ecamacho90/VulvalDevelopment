function dydt=cusp_and_saddlenode_model_v10(t,y,p)


%Seventh version of the model. Here we use a general plane like in v5 but we
%also use an exponential decay in time of EGF and Notch signal. We also
%introduce a downregulation of Notch sensitivity in 1st fated cells.


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
% n0 = -1; 14  %Notch sigmoidal function as in Siggia 2012
% n1 = 2;  15
% n1x = -1; 16
% n2x = 0; 17
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
% L = 1         29 % Downregulation strength of Notch in 1st fated cells
% H3 ;          30 % Sigmoid function EGF upregulation during competence
% M3 ;          31 % Sigmoid function EGF upregulation during competence
% C  ;          32 % Constant corresponding to max value sigmoid for
%                    postcompetence
%tAC ;          33 % AC ablation time
%t1 ;           34 % End time competence
% EGFover;      35  Over expression EGF
% D             36 Noise
% p=[Tau,H,M,m11,m12,m21,m22,m31,m32,q1,q2,q3,gamma,H1,M1,H2,M2,s1,s2,s3,l1,l2,l3,alpha,lambdaE,lambdaN,xiE,xiN,L,H3,M3,C,tAC,t1];

dydt=(p(1))*[
((1+tanh(p(2)*(-y(2)+(p(3)*(sqrt(2)-1)/(2*sqrt(2))))))/2)      *        (-4*y(1)^3      -       2*(  p(4)*p(18)*(p(13)^2)*exp(-p(25)*(t-p(33)))*exp(-p(27)*(t-p(34)))*p(32)*(1+tanh(p(30)*t+p(31)))/2+                  p(5)*p(21)*exp(-p(26)*(t-p(33)))*exp(-p(28)*(t-p(34)))*   (1-p(29)*((1+tanh(p(14)+p(15)*(p(16)*y(1)+p(17)*y(2))))/(2))) * ( p(24) *((1+tanh(p(14)+p(15)*(p(16)*y(1)+p(17)*y(2))))/(2))+    ((1+tanh(p(14)+p(15)*(p(16)*y(3)+p(17)*y(4))))/(2)))                                                        +  p(10)) *y(1)                -          ( p(6)*p(18)*(p(13)^2)*exp(-p(25)*(t-p(33)))*exp(-p(27)*(t-p(34)))*p(32)*(1+tanh(p(30)*t+p(31)))/2  +  p(7)*p(21)*exp(-p(26)*(t-p(33)))*exp(-p(28)*(t-p(34)))* (1-p(29)*((1+tanh(p(14)+p(15)*(p(16)*y(1)+p(17)*y(2))))/(2))) * ( p(24) *((1+tanh(p(14)+p(15)*(p(16)*y(1)+p(17)*y(2))))/(2)) +   ((1+tanh(p(14)+p(15)*(p(16)*y(3)+p(17)*y(4))))/(2)))                                                        +  p(11))                )               -              (1-((1+tanh(p(2)*(-y(2)+(p(3)*(sqrt(2)-1)/(2*sqrt(2))))))/2))*y(1); ...
-y(2)*        ((y(2)-p(3))^2    +     ( p(8)*p(18)*(p(13)^2)*exp(-p(25)*(t-p(33)))*exp(-p(27)*(t-p(34)))*p(32)*(1+tanh(p(30)*t+p(31)))/2  +  p(9)*p(21)*exp(-p(26)*(t-p(33)))*exp(-p(28)*(t-p(34)))* (1-p(29)*((1+tanh(p(14)+p(15)*(p(16)*y(1)+p(17)*y(2))))/(2))) * ( p(24) *((1+tanh(p(14)+p(15)*(p(16)*y(1)+p(17)*y(2))))/(2)) +   ((1+tanh(p(14)+p(15)*(p(16)*y(3)+p(17)*y(4))))/(2)))                                                        +  p(12))    ); ...
((1+tanh(p(2)*(-y(4)+(p(3)*(sqrt(2)-1)/(2*sqrt(2))))))/2)      *        (-4*y(3)^3      -       2*(  p(4)*p(19)*p(13)    *exp(-p(25)*(t-p(33)))*exp(-p(27)*(t-p(34)))*p(32)*(1+tanh(p(30)*t+p(31)))/2+                  p(5)*p(22)*exp(-p(26)*(t-p(33)))*exp(-p(28)*(t-p(34)))*   (1-p(29)*((1+tanh(p(14)+p(15)*(p(16)*y(3)+p(17)*y(4))))/(2))) * ( p(24) *((1+tanh(p(14)+p(15)*(p(16)*y(3)+p(17)*y(4))))/(2))+    ((1+tanh(p(14)+p(15)*(p(16)*y(1)+p(17)*y(2))))/(2)) + ((1+tanh(p(14)+p(15)*(p(16)*y(5)+p(17)*y(6))))/(2)))  +  p(10)) *y(3)                -          ( p(6)*p(19)*p(13)    *exp(-p(25)*(t-p(33)))*exp(-p(27)*(t-p(34)))*p(32)*(1+tanh(p(30)*t+p(31)))/2  +  p(7)*p(22)*exp(-p(26)*(t-p(33)))*exp(-p(28)*(t-p(34)))* (1-p(29)*((1+tanh(p(14)+p(15)*(p(16)*y(3)+p(17)*y(4))))/(2))) * ( p(24) *((1+tanh(p(14)+p(15)*(p(16)*y(3)+p(17)*y(4))))/(2)) +   ((1+tanh(p(14)+p(15)*(p(16)*y(1)+p(17)*y(2))))/(2)) + ((1+tanh(p(14)+p(15)*(p(16)*y(5)+p(17)*y(6))))/(2)))  +  p(11))                )               -              (1-((1+tanh(p(2)*(-y(4)+(p(3)*(sqrt(2)-1)/(2*sqrt(2))))))/2))*y(3); ...
-y(4)*        ((y(4)-p(3))^2    +     ( p(8)*p(19)*p(13)    *exp(-p(25)*(t-p(33)))*exp(-p(27)*(t-p(34)))*p(32)*(1+tanh(p(30)*t+p(31)))/2  +  p(9)*p(22)*exp(-p(26)*(t-p(33)))*exp(-p(28)*(t-p(34)))* (1-p(29)*((1+tanh(p(14)+p(15)*(p(16)*y(3)+p(17)*y(4))))/(2))) * ( p(24) *((1+tanh(p(14)+p(15)*(p(16)*y(3)+p(17)*y(4))))/(2)) +   ((1+tanh(p(14)+p(15)*(p(16)*y(1)+p(17)*y(2))))/(2)) + ((1+tanh(p(14)+p(15)*(p(16)*y(5)+p(17)*y(6))))/(2)))  +  p(12))    ); ...
((1+tanh(p(2)*(-y(6)+(p(3)*(sqrt(2)-1)/(2*sqrt(2))))))/2)      *        (-4*y(5)^3      -       2*(  p(4)*p(20)          *exp(-p(25)*(t-p(33)))*exp(-p(27)*(t-p(34)))*p(32)*(1+tanh(p(30)*t+p(31)))/2+                  p(5)*p(23)*exp(-p(26)*(t-p(33)))*exp(-p(28)*(t-p(34)))*   (1-p(29)*((1+tanh(p(14)+p(15)*(p(16)*y(5)+p(17)*y(6))))/(2))) * ( p(24) *((1+tanh(p(14)+p(15)*(p(16)*y(5)+p(17)*y(6))))/(2))+  2*((1+tanh(p(14)+p(15)*(p(16)*y(3)+p(17)*y(4))))/(2)))                                                        +  p(10)) *y(5)                -          ( p(6)*p(20)          *exp(-p(25)*(t-p(33)))*exp(-p(27)*(t-p(34)))*p(32)*(1+tanh(p(30)*t+p(31)))/2  +  p(7)*p(23)*exp(-p(26)*(t-p(33)))*exp(-p(28)*(t-p(34)))* (1-p(29)*((1+tanh(p(14)+p(15)*(p(16)*y(5)+p(17)*y(6))))/(2))) * ( p(24) *((1+tanh(p(14)+p(15)*(p(16)*y(5)+p(17)*y(6))))/(2)) + 2*((1+tanh(p(14)+p(15)*(p(16)*y(3)+p(17)*y(4))))/(2)))                                                        +  p(11))                )               -              (1-((1+tanh(p(2)*(-y(6)+(p(3)*(sqrt(2)-1)/(2*sqrt(2))))))/2))*y(5); ...
-y(6)*        ((y(6)-p(3))^2    +     ( p(8)*p(20)          *exp(-p(25)*(t-p(33)))*exp(-p(27)*(t-p(34)))*p(32)*(1+tanh(p(30)*t+p(31)))/2  +  p(9)*p(23)*exp(-p(26)*(t-p(33)))*exp(-p(28)*(t-p(34)))* (1-p(29)*((1+tanh(p(14)+p(15)*(p(16)*y(5)+p(17)*y(6))))/(2))) * ( p(24) *((1+tanh(p(14)+p(15)*(p(16)*y(5)+p(17)*y(6))))/(2)) + 2*((1+tanh(p(14)+p(15)*(p(16)*y(3)+p(17)*y(4))))/(2)))                                                        +  p(12))    )];











