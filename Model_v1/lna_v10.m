function [DetSol,CovMat]=lna_v10(y0,t0,t1,M,CovMat0,param,dimension)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          v10 cusp-saddlenode                             %
%                                                                         %
%  This programme computes the LNA for the SDEs and returns the mean      %
%  trajectory and covariance matrix at the end point of the trajectory.   %
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


    detsolstruct=solution_det_eqtn_v10(y0,tspan,param);
    
    DetSol=detsolstruct.y(:,end);
    
    CovMatStruct=covariancemat_direct_solution_v10(CovMat0,detsolstruct,tspan,param,dimension);
    
    CovMat=CovMatStruct.y(:,end);
    
    CovMat=reshape(CovMat,dimension,dimension);
    
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Elena 04/22/20
