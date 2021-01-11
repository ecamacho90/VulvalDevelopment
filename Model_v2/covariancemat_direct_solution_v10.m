function CovMatSol=covariancemat_direct_solution_v10(CovMat0,detsolstruct,tspan,param,dimension)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          v9 cusp-saddlenode                             %
%                                                                         %
%  This programme gives solves the ODE of the LNA which solution is the   %
%  covariance matrix of the Wiener process at each time.                  %
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


CovMat0=reshape(CovMat0,dimension^2,1);


CovMatSol=ode15s(@(t,CovMat) covariance_matrix_ode_v10(CovMat,detsolstruct,t,param,dimension),tspan,CovMat0);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Elena 08/11/18

