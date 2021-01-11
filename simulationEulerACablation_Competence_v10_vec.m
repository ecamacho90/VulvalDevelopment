function sol=simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,parametersmodel,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%                    cusp-saddlenode                                     %
%                                                                        %
%  By using the Euler Scheme:                                            %
%  This programme compute nsimulations simulations, and for each         %
%  simulation we discretise the time interval into M subintervals.       %
%  It only saves the last time point.                                    %
%  The trajectory for P4.p will be Xprev(1:2), P5.p will be Xprev(3:4),  %
%  P6.p will be Xprev(5:6).                                              %
%                                                                        %
%  sol is a matrix in which the element i,j is the proportion of times   %
%  that cell i took fate j in the nsimulations. If j=4, that means that  %
%  the fate wasn't able to be computed.                                  %
%                                                                        %
%  In this version we take a new input variable tAC which is the AC      %
%  ablation time.                                                        %
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modeleqstr='cusp_and_saddlenode_model_v10_vec';

modeleq=str2func(modeleqstr);



%Post Competence parameters:
%---------------------------
%     saddlefate3 = param(3)-sqrt(-cPC);
%     equilibriax = equilibria_fates_1_2(bPC,DiscriminantPC);
  
    
%Function Check:
%---------------

%    dt=(t1-t0)/M;                %Time step for each simulation.
    
%     y0=[0;4;0;4;0;4];   %Initial condition for the ODE.

    
% %     %In case we want to do a computation without having to give the noise:
% %         rng(13)
% %         Z1D=randn(nsimulations,(M))*sqrt(2*D*dt);  %matrix of independent random deviations...
% %         Z2D=randn(nsimulations,(M))*sqrt(2*D*dt);  %Each row contains the random deviations of one
% %         Z3D=randn(nsimulations,(M))*sqrt(2*D*dt);  %simulation. Each column corresponds to a time step.
% %         Z4D=randn(nsimulations,(M))*sqrt(2*D*dt);
% %         Z5D=randn(nsimulations,(M))*sqrt(2*D*dt);
% %         Z6D=randn(nsimulations,(M))*sqrt(2*D*dt);
% % 
% %         Z0D=repmat(y0,1,nsimulations)+randn(6,nsimulations)*sqrt(D);  %Initial condition for the SDE.
% % 
% %         Z1DPC=randn(nsimulations,(M))*sqrt(2*D*dt);  %matrix of independent random deviations...
% %         Z2DPC=randn(nsimulations,(M))*sqrt(2*D*dt);  %Each row contains the random deviations of one
% %         Z3DPC=randn(nsimulations,(M))*sqrt(2*D*dt);  %simulation. Each column corresponds to a time step.
% %         Z4DPC=randn(nsimulations,(M))*sqrt(2*D*dt);
% %         Z5DPC=randn(nsimulations,(M))*sqrt(2*D*dt);
% %         Z6DPC=randn(nsimulations,(M))*sqrt(2*D*dt);
% % 
% %         Z0D=repmat(y0,1,nsimulations)+randn(6,nsimulations)*sqrt(D);  %Initial condition for the SDE.

    simulateddata=zeros(3,4);    %This matrix contains the proportions of times that...
                                 %each of the three cells took one of...
                                 %the 3 fates (we add a 4th in case...
                                 %we can't determine it).
                                 
% tic    
%Simulations:
%------------
steps=M;
stepsPC = MPC;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   NO ANCHOR CELL ABLATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (tAC>=t1)  %If tAC is greater or equal than t1  there is no AC ablation, thus we compute the fates with respect to the normal parameters
    
    %Initial conditions:
    %--------------------

    param=parametersmodel(1:(end-1));
    param(25:28)=zeros(1,4); %No decay
    param(32) = 1;
    param(33:34) = 0;
    
    D=parametersmodel(end);      %Diffusion constant.
    
    paramPC = parametersmodel(1:(end-1)); %Parameters Post-Competence (decay)
    paramPC(30:31) = zeros(1,2); %Sigmoidal function must be fixed
    paramPC(32) = 2*(1+tanh(parametersmodel(30)*t1+parametersmodel(31)))/2; %C is equal to what the sigmoidal function was
    paramPC(27:28) =zeros(1,2); %No ablation decay
    paramPC(33) = t1;
    param(34) = 0; %Good because paramPC(27)=0;
    


      Xprev=Z0D;    %Set the initial point for each cell at the new simulation
    
      %%%%%%%%%%%%%%%%%%%%%
      % COMPETENCE:
      %%%%%%%%%%%%%%%%%%%%%
        for timestep=1:steps
        
            dtfprev=dt*modeleq(t0+(timestep-1)*dt,Xprev,param);       %Evaluate f at the previous state
        
            Xprev=Xprev+dtfprev+[Z1D(timestep,:);Z2D(timestep,:);Z3D(timestep,:);Z4D(timestep,:);Z5D(timestep,:);Z6D(timestep,:)];  %Renew the state taking an Euler Step
    

%           if (Xprev+param(1))<0  %If x0 is so small that there is a solution such that xprev+x0 is negative, we don't accept the value of x0.
%                 sol=zeros(3,3);
%               return
%           end
          
        end
        
      %%%%%%%%%%%%%%%%%%%%%
      % POST - COMPETENCE:
      %%%%%%%%%%%%%%%%%%%%%
    
        for timestep = 1:stepsPC
        
            dtfprev=dtPC*modeleq(t1+(timestep-1)*dtPC,Xprev,paramPC);       %Evaluate f at the previous state
        
            Xprev=Xprev+dtfprev+[Z1DPC(timestep,:);Z2DPC(timestep,:);Z3DPC(timestep,:);Z4DPC(timestep,:);Z5DPC(timestep,:);Z6DPC(timestep,:)];  %Renew the state taking an Euler Step
    
        
        end
      %%%%%%%%%%%%%%%%%%%%%
      % FATES:
      %%%%%%%%%%%%%%%%%%%%%    
       for simulation=1:nsimulations
           Xprevsim = Xprev(:,simulation);
        if norm(Xprevsim)<10^10
            fatessimulated = computefates_PostCompetence(Xprevsim,param(2),param(3),aPC,bPC,cPC,DiscriminantPC,saddlefate3,equilibriax);  %Compute the fates corresponding to the end point of the simulated trajectory after competence
        else
            fatessimulated = [0 0 0 1; 0 0 0 1; 0 0 0 1];   
        end

            simulateddata=simulateddata+fatessimulated;   %Save the simulated fates.

      end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ANCHOR CELL ABLATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
else  %Otherwise there is AC ablation and we compute the fates with respect to the parameters where EGF is set to 0.
    
    %Initial conditions:
    %--------------------

    param=parametersmodel(1:(end-1));
    param(25:28)=zeros(1,4); %No decay
    param(32) = 1; % C=1
    param(33:34) = 0; % no t1 or tAC parameters
    
    D=parametersmodel(end);      %Diffusion constant.
    
    paramafterACablation = parametersmodel(1:(end-1)); %Parameters AC ablation (decay from original EGF value)
    paramafterACablation(25:26) = zeros(1,2); % No Post-competence decay
    paramafterACablation(30:31) = zeros(1,2); % No sigmoidal increase
    paramafterACablation(32) = 2*(1+tanh(parametersmodel(30)*tAC+parametersmodel(31)))/2; %EGF scaling constant corresponding to gradual increment before AC ablation
    paramafterACablation(34) = tAC;
    paramafterACablation(33) = 0;
    paramafterACablation(28) = 0; % No decay of Notch due to AC ablation
    
    paramPC = parametersmodel(1:(end-1)); %Parameters Post-Competence (decay)
    paramPC(30:31) = zeros(1,2); %Sigmoidal function must be fixed
    paramPC(32) = (2*(1+tanh(parametersmodel(30)*tAC+parametersmodel(31)))/2)*exp(-parametersmodel(27)*(t1-tAC)); %C is equal to what the sigmoidal function was
    paramPC(27:28) =zeros(1,2); %No ablation decay
    paramPC(33) = t1;
    param(34) = 0;  
    
    
 stepsbeforeAC=min([floor((tAC-t0)*M/(t1-t0)),M]);  %We consider that there is EGF until the previous time step to tAC. If tAC is bigger than t1, then the maximum number of steps is M
    
    
        Xprev=Z0D;     %Set the initial point for each cell at the new simulation
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % COMPETENCE PRE-AC Ablation:
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for timestep=1:stepsbeforeAC
            
            dtfprev=dt*modeleq(t0+(timestep-1)*dt,Xprev,param);       %Evaluate f at the previous state

            Xprev=Xprev+dtfprev+[Z1D(timestep,:);Z2D(timestep,:);Z3D(timestep,:);Z4D(timestep,:);Z5D(timestep,:);Z6D(timestep,:)];  %Renew the state taking an Euler Step
    
%             if (Xprev+param(1))<0  %If x0 is so small that there is a solution such that xprev+x0 is negative, we don't accept the value of x0.
%                 sol=zeros(3,3);
%                 return
%             end
        
        end
        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % COMPETENCE POST-AC Ablation:
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for timestep=(stepsbeforeAC+1):steps
 
            dtfprev=dt*modeleq(t0+(timestep-1)*dt,Xprev,paramafterACablation);       %Evaluate f at the previous state
        
            Xprev=Xprev+dtfprev+[Z1D(timestep,:);Z2D(timestep,:);Z3D(timestep,:);Z4D(timestep,:);Z5D(timestep,:);Z6D(timestep,:)];  %Renew the state taking an Euler Step
    
%             if (Xprev+param(1))<0  %If x0 is so small that there is a solution such that xprev+x0 is negative, we don't accept the value of x0.
%                 sol=zeros(3,3);
%                 return
%            end
        
        
        end
        
      %%%%%%%%%%%%%%%%%%%%%
      % POST - COMPETENCE:
      %%%%%%%%%%%%%%%%%%%%%
    
        for timestep = 1:stepsPC
        
            dtfprev=dtPC*modeleq(t1+(timestep-1)*dtPC,Xprev,paramPC);       %Evaluate f at the previous state
        
            Xprev=Xprev+dtfprev+[Z1DPC(timestep,:);Z2DPC(timestep,:);Z3DPC(timestep,:);Z4DPC(timestep,:);Z5DPC(timestep,:);Z6DPC(timestep,:)];  %Renew the state taking an Euler Step
    
        
        end
        
      %%%%%%%%%%%%%%%%%%%%%
      % FATES:
      %%%%%%%%%%%%%%%%%%%%%
    for simulation=1:nsimulations
        Xprevsim = Xprev(:,simulation);
        if norm(Xprevsim)<10^10
            fatessimulated = computefates_PostCompetence(Xprevsim,param(2),param(3),aPC,bPC,cPC,DiscriminantPC,saddlefate3,equilibriax);  %Compute the fates corresponding to the end point of the simulated trajectory after competence
        else
            fatessimulated = [0 0 0 1; 0 0 0 1; 0 0 0 1];   
        end
    
        simulateddata=simulateddata+fatessimulated;   %Save the simulated fates.
    
    end
    
end
    

%simulateddata(:,4)

    sol=simulateddata/nsimulations;   %Average the fates for the simulations
% toc
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Elena 29/9/16   (Now the unclassified fates are taken into account in the
%                 fourth column of sol)
%                (If the trayectory diverges, we consider that we couldn't
%                compute the fate)
%      27/9/17  Now we add the competence and post-competence period
%      07/11/17 Now we change the model to the version 4 in which we have
%               changed the Notch function to be like Siggia's one.
%      28/11/17 Now we change the model to the version 5 in which we
%               consider a more general plane or change of coordinates from Signal
%               space to control space.
%      01/12/17 Now we change the model to the version 6 in which we
%               consider a more general plane or change of coordinates from Signal
%               space to control space, and also added a signal decay for
%               Post competence period and AC ablation in EGF and Notch.
%               We have also added MPC and tPC to control how long is the
%               Post-competence period and how small is the time step.
%      18/01/18 Now we change the model to the version 7 in which we
%               consider a more general plane or change of coordinates from Signal
%               space to control space, and also added a signal decay for
%               Post competence period and AC ablation in EGF and Notch,
%               and also we have added a downregulation of Notch for cells
%               in 1st fate.
%      06/11/18 Now we change the model to the version 9 in which we
%               add a sigmoid function to get a gradual increment in EGF
%               signal when the simulation starts. We have also fixed the
%               values of the parameters after AC ablation and competence
%               to take into account that the EGF value should start from
%               the EGF value in the previous stage. We also changes the
%               AC ablation and post-competence decay functions so that
%               they start at 1 and finish close to 0.
%      22/04/20 Change it to use Model version 10, which is model 9 but the Notch
%               sigmoidal is taken as in Siggia,2012. We fit the parameters of this
%               sigmoidal function (so, in total, 16 parameters). Also fixed npar in
%               ABC_SMC function, which was set to 37 instead of 36,so when T>1, the noise was equal to
%               the distance to mutant 1!!
%