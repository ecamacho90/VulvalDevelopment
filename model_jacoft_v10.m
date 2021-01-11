function J=model_jacoft_v10(detsolstruct,t,param)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                          v6 cusp-saddlenode                             %
%                                                                         %
%  This programme gives the value of the jacobian of the model at time t. %
%                                                                         %
%  detsolstruct is the structure obtained from the deterministic solution %
%  The output is J, which is Jacobian(y(t),t,param).                      %                                                
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modeleqstr='cusp_and_saddlenode_model_singlecell_v10';

%Fist, we compute y(t):

    yt=deval(detsolstruct(1),t);

    
%We need to evaluate the jacobian:
    modeleqjac=strcat(modeleqstr,'_jac');
    
%Evaluate the jacobian matrix in that point:
    
    J=feval(modeleqjac,t,yt,param);
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Elena 04/22/20
        

   