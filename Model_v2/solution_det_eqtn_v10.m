function sol=solution_det_eqtn_v10(y0,tspan,param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          v9 cusp-saddlenode                             %
%                                                                         %
%  This programme solves the deterministic differential equation of the   %
%  model for a single cell without signals.                               %
%  y0 is the initial condition.                                           %
%  tspan is the discretisation of time in which we approximate the        %
%        solution.                                                        %
%  param are the parameters                                               %
%                                                                         %
%  sol is the structure with the solution of the ODE.                     %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





sol=ode15s(@(t,y) cusp_and_saddlenode_model_singlecell_v10(t,y,param),tspan,y0);%,options);
    





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Elena 22/04/20 (New Model: v10, general plane + exponential decay + downregulation of Notch + sigmoidal EGF upregulation + Siggia's Notch Sigmoidal function)

