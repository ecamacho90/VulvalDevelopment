function [DetSolStruct,CovMatStruct]=lnaStruct_v10(y0,t0,t1,M,CovMat0,param,dimension)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          v9 cusp-saddlenode                             %
%                                                                         %
%  This programme computes the LNA for the SDEs and returns the ode       %
%  structure of the mean and covariance matrix, so that we can then       %
%  evaluate it at the AC ablation times and use them as the initial       %
%  condition for the trajectory after AC ablation.                        %
%                                                                         %
%  CovMat0 is the initial covariance matrix at t=tspan(1).                %
%  detsolstruct is the structure of y(t).                                 %
%  tspan is the time interval discretisation.                             %
%  param are the parameters.                                              %
%  dimension is the number of rows of the covariance matrix.              %
%  The output is CovMatSol, which is the structure containing the         %
%  solution of the ODE.                                                   %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tspan=linspace(t0,t1,M+1);

    
    DetSolStruct=solution_det_eqtn_v10(y0,tspan,param); 
    
    CovMatStruct=covariancemat_direct_solution_v10(CovMat0,DetSolStruct,tspan,param,dimension);
    
%     CovMat=CovMatStruct.y(:,end);
%     
%     CovMat=reshape(CovMat,dimension,dimension);
    
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Elena 04/22/20
