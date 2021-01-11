function ABC_SMC_Algorithm_OLCM_Template_Modelv10_v1_Fates_AbsDistance_L(mutantstofit,parfitnumbers,EpT,T,Tmax,Nmax,Numberofparticleseachjob,PreviousData,namenewdata,CovarianceMatricesParticles,weightsdistances,streamnum,jobnum)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%                         ABC SMC ALGORITHM
%                   Optimal Local Covariance Matrix 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               INPUT VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% * *mutantstofit*: Vector containing the number of the mutants that we are
% going to fit.
% * *parfitnumbers*: Vector containing the number of the parameters that we
% are going to fit.
% * *EpT*: Vector containing the thresholds that we are going to use for
% the algorithm. 
% * *T*: Number of steps that we have done in the thresholds. For computing the first
% threshold, T=0. For the second threshold, T=1.
% * *Tmax*: Number of maximum number of steps in thresholds. For computing
% one threshold, Tmax=1. For the second one, Tmax=2.
% * *Nmax*: Number of total particles that we are considering.
% * *jobnum*: Number of the parallel subset of particles that we are
% computing.
% * *Numberofparticleseachjob*: Number of particles that we will compute in
% this job.
% * *PreviousData*: Matrix containing the particles and distances from the
% previous iteration.
% * *namenewdata* : Name of the Workspace of the new data set.
% * *CovarianceMatricesParticles*: Local covariance matrix of the normal
% perturbation kernel for each particle.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cd('/cluster/elenacamacho')

%Maximum number of steps:
MaxStep = 1.e6;

%Set the stream of random number generators:
%-------------------------------------------

reset(RandStream.getGlobalStream) %reset the random stream

substreamnum = jobnum; %number of substream

cmrg = RandStream.create('mrg32k3a','NumStreams',100,'StreamIndices',streamnum);
RandStream.setGlobalStream(cmrg)
cmrg.Substream = substreamnum;



%Time for initial condition:
t0lna = 0;
t1lna = 10;
dtlna = 0.1;
tspanlna=t0lna:dtlna:t1lna;

%Number of subintervals:
Mlna = length(tspanlna)-1;

%Time step:
dt = 0.005;

%Time for solution approximation:
t0 = 0;
t1 = 1;
tspan=t0:dt:t1;

%Number of subintervals:
M = length(tspan)-1;

%Time step Post Competence:
dtPC = 0.005;

%Time for solution approximation:
t0PC = t1;
t1PC = 3;
tspanPC=t0PC:dtPC:t1PC;

%Number of subintervals:
MPC = length(tspanPC)-1;

%Number of simulations:
nsimulations = 150;%64;

%Number of parameters:
npar = 36;

%      1  2 3  4   5   6   7   8   9  10 11 12  13   14 15 16 17 18 19 20
%  p=[Tau,H,M,m11,m12,m21,m22,m31,m32,q1,q2,q3,gamma,n0,n1,n1x,n1y,s1,s2,s3,

%  21 22 23 24     25        26    27 28  29 30  31  32  33  34   35     36
%  l1,l2,l3,alpha,lambdaE,lambdaN,xiE,xiN,L, H3, M3, C, tAC, t1,  EGFover,D];

%Number of parameters to fit:
nparfit = length(parfitnumbers);

%Initial parameters:
initpar = [1.0000   25.0000    4.0000   -0.5000 -0.3    0.6000   -0.6000    0.5000 0.5000   -0.5000         0   -0.5000 0.0500    -2   2    -1 0   10   10   10 2.63 2.63 2.63    0.76 4.0000    4.0000   50.0000   50.0000 1 6  -1  0 t1 t0PC  3.5882    0.05];
minimumM = (2*sqrt(2)/(sqrt(2)-1))*(atanh(1-10^(-12))/initpar(2));

%Model name:
model = 'Vulval_Development_Modelv10_v1';

%Initial condition:
initcond = [0;initpar(3)];  %INSTEAD OF THIS POINT WE SHOULD TAKE THE POINT [0,M,0,M,0,M] SO THAT WE ARE SURE THAT IT'S IN THE THIRD FATE BASIS OF ATTRACTION


%Experimental data:
probabilitiesdataWT=[0 0 1 0;0 1 0 0;1 0 0 0];
probabilitiesdataNoNotch=[0 0 1 0;0 0 1 0;1 0 0 0];
probabilitiesdataNotchNull2AC=[0 0 1 0;1 0 0 0;1 0 0 0];
probabilitiesdataEGFover=[0.18 0.46 0.36 0;0.455 0.545 0 0;0.96 0.04 0 0];
probabilitiesdataACablation1=[0 0 1 0;0 0 1 0];
probabilitiesdataACablation2=[0.015 0.21 0.775 0;0.18 0.18 0.64 0];
probabilitiesdataACablation3=[0 0.5463 0.4537 0;0.31 0.38 0.31 0];
probabilitiesdataACablation4=[0.04 0.9 0.06 0;0.52 0.48 0 0];
probabilitiesdataACablation5=[0.01 0.99 0 0;0.65 0.35 0 0];
probabilitiesdataACablation6=[0.01 0.99 0 0;0.93 0.07 0 0];

expdata{1,1}=probabilitiesdataWT;
expdata{1,2}=probabilitiesdataWT;
expdata{1,3}=probabilitiesdataWT;
expdata{1,4}=probabilitiesdataWT;
expdata{1,5}=probabilitiesdataNotchNull2AC;
expdata{1,6}=probabilitiesdataNoNotch;
expdata{1,7}=probabilitiesdataEGFover;
expdata{1,8}=probabilitiesdataACablation1;
expdata{1,9}=probabilitiesdataACablation2;
expdata{1,10}=probabilitiesdataACablation3;
expdata{1,11}=probabilitiesdataACablation4;
expdata{1,12}=probabilitiesdataACablation5;
expdata{1,13}=probabilitiesdataACablation6;


%Number of mutants:
nmutants = length(mutantstofit);

% save(['paramFatesM5v1_',num2str(jobnum)])


%DON'T NEED TO CHANGE:
%---------------------

%Priors function handle:
priorshandle = str2func(strcat(model,'_Priors'));

%Simulated data function handle:
modelhandle = str2func(strcat(model,'_simulation'));

%Distance function handle:
distancehandle = str2func(strcat(model,'_AbsDistance_Fates'));

%Covariance matrices OLCM function handle:
covarianceOLCMhandle = str2func(strcat(model,'_CovarianceMatricesOLCM'));

%Perturbation Kernels handle:
kernelshandle = str2func(strcat(model,'_Kernels'));

%Density function prior handle:
evalpriorshandle = str2func(strcat(model,'_EvalPriors'));

%Density function kernel handle:
evalkernelhandle = str2func(strcat(model,'_EvalKernels'));

%Density function constraints handle:
constraintshandle = str2func(strcat(model,'_Constraints'));

%Relations between parameters contraints handle:
relationsconstraintshandle1 = str2func(strcat(model,'_Relations_Constraints1'));

%Relations between parameters contraints handle:
relationsconstraintshandle2 = str2func(strcat(model,'_Relations_Constraints2'));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SAMPLER T=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
if T==0
    
    tic

    
    %Update threshold number:
    %------------------------
    
    T = T+1; %This avoids that the next if T<Tmax is done.
    N = 0;
    i = 0;
    icomp = 0;
    
    %Set new threshold so that the threshold for each mutant is EpT(T):
    %------------------------------------------------------------------
    newEpT = EpT(T);
    
    %New cell in the particles matrix:
    %---------------------------------
    NewData = [];
    NewFates = zeros(3,4,nmutants,Numberofparticleseachjob); %NewFates(i,j,k,m) is the probability of cell i taking fate j for experiment k when using particle m.
    

    %Find the particles:
    %-------------------
    while (N<Numberofparticleseachjob)&&(i<MaxStep)

        if rem(i,100)==0
            disp('---------------------------------------------------')
            disp('---------------------------------------------------')
            disp(strcat('            Threshold: ',num2str(EpT(T))));
            disp(strcat('            Job number: ',num2str(jobnum)));
            disp(strcat('            Step: ',num2str(i)));
            disp(strcat('            Step Comp: ',num2str(icomp)));
            disp(strcat('            Accepted particles: ',num2str(N)));
            disp('---------------------------------------------------')
            disp('---------------------------------------------------')
        end


        %Update steps:
        %-------------

        i = i+1;

        %Sample parameter values:
        %------------------------
            
        valuepriors = zeros(1,nparfit);  % vector to save parameters values

        paramaux = initpar;              % in paramaux we will substitute the new values for the parameters we are going to fit

        parami = 0;                      % parami is a counter that will go from 0 to nparfit-1   

        priorOK = 1;                     % priorOK will become 0 when a parameter doesn't satisfy the prior distribution   

        while (parami<nparfit)&&(priorOK)  %Only if priorOK is not 0 and parami<npar the while loop will stop
            
%             minimumM = (2*sqrt(2)/(sqrt(2)-1))*(atanh(1-10^(-12))/paramaux(2));
            parami = parami+1;

            paramifit = parfitnumbers(parami); %Parameter number

            provparvalue = feval(priorshandle,paramifit,minimumM); %Sample a value

            priorOK = feval(constraintshandle,provparvalue,paramifit,minimumM); %Check that is a valid value

            paramaux(paramifit) = provparvalue; %Save the value

        end
            
        %Compute cPC and DiscriminantPC:
        %-------------------------------
        aPC = paramaux(10);
        bPC = paramaux(11);
        cPC = paramaux(12);
        DiscriminantPC = (8*aPC^3+27*bPC^2);
        
        %Check cPC<0 and DiscriminantPC<0
        %-------------------------------
        if priorOK
            
            if (cPC<0)&&(DiscriminantPC<0)
                priorOK = 1;
                
            else
                priorOK = 0;
                
            end
            
        end
        
        %Compute n1x as function of n1y
        %------------------------------
        
        paramaux(16) = -sqrt(1-paramaux(17)^2);

        
        %Check that s1=s2=s3 and l1=l2=l3:
        %---------------------------------
        
        paramaux(18:20) = paramaux(18)*ones(1,3);
        paramaux(21:23) = paramaux(21)*ones(1,3);
        
        
        %Check that EGF and Notch decay is big enough
        %--------------------------------------------
%         disp('----')
%         priorOK
        if priorOK
            
            if (paramaux(18)*exp(-paramaux(25)*(t1PC-t1))<0.001)&&(paramaux(21)*exp(-paramaux(26)*(t1PC-t1))<0.001)
                priorOK = 1;
                paramaux(27)=paramaux(25);         
            else
                priorOK = 0;
                
            end
            
        end
        
        
        %Write m31 as functions of the other ms:
        %--------------------------------
        paramaux(8) = sqrt(1-paramaux(4)^2-paramaux(6)^2);

        
        
        
%         %Check that H1=H2:
%         %-----------------
%         paramaux(13) = paramaux(11);

        
        %Check that the parameters fulfill all the constraints:
        %------------------------------------------------------
        
        if priorOK
            
            priorOK = feval(relationsconstraintshandle1,paramaux,minimumM,cPC);
            
        end

        
        if priorOK
            % Find m22 and m32 as functions of the other ms:
            %-----------------------------------------------
            [paramaux(7),paramaux(9),flag] = findm22m32(paramaux(4),paramaux(5),paramaux(6),paramaux(8));
        
           priorOK=flag;
           
        end

        %Check that the parameters fulfill all the constraints:
        %------------------------------------------------------
        
        if priorOK
            
            priorOK = feval(relationsconstraintshandle2,paramaux);
            
        end
       
        if priorOK
            icomp=icomp+1;
            
            %Post Competence parameters:
            %---------------------------
             saddlefate3 = paramaux(3)-sqrt(-cPC);
             equilibriax = equilibria_fates_1_2(bPC,DiscriminantPC);
            
            %COMPUTE THE NOISE:
            %------------------
            
            D = paramaux(end);
             
            noisecov=sqrt(2*dt*D);
             
%             Z1D = randn(nsimulations,(M))*noisecov;  %matrix of independent random deviations...
%             Z2D = randn(nsimulations,(M))*noisecov;  %Each row contains the random deviations of one
%             Z3D = randn(nsimulations,(M))*noisecov;  %simulation. Each column corresponds to a time step.
%             Z4D = randn(nsimulations,(M))*noisecov;
%             Z5D = randn(nsimulations,(M))*noisecov;
%             Z6D = randn(nsimulations,(M))*noisecov;
%             
%             Z1DPC = randn(nsimulations,(MPC))*noisecov;  %matrix of independent random deviations...
%             Z2DPC = randn(nsimulations,(MPC))*noisecov;  %Each row contains the random deviations of one
%             Z3DPC = randn(nsimulations,(MPC))*noisecov;  %simulation. Each column corresponds to a time step.
%             Z4DPC = randn(nsimulations,(MPC))*noisecov;
%             Z5DPC = randn(nsimulations,(MPC))*noisecov;
%             Z6DPC = randn(nsimulations,(MPC))*noisecov;

%Vectorial:
            Z1D = randn((M),nsimulations)*noisecov;  %matrix of independent random deviations...
            Z2D = randn((M),nsimulations)*noisecov;  %Each row corresponds to a time step.
            Z3D = randn((M),nsimulations)*noisecov;  %Each column contains the random deviations of one simulation.
            Z4D = randn((M),nsimulations)*noisecov;
            Z5D = randn((M),nsimulations)*noisecov;
            Z6D = randn((M),nsimulations)*noisecov;
            
            Z1DPC = randn((MPC),nsimulations)*noisecov;  %matrix of independent random deviations...
            Z2DPC = randn((MPC),nsimulations)*noisecov;  %Each row contains the random deviations of one
            Z3DPC = randn((MPC),nsimulations)*noisecov;  %simulation. Each column corresponds to a time step.
            Z4DPC = randn((MPC),nsimulations)*noisecov;
            Z5DPC = randn((MPC),nsimulations)*noisecov;
            Z6DPC = randn((MPC),nsimulations)*noisecov;
            
            N01 = randn(6,nsimulations);
             

            % REAL INITIAL CONDITION:
            %------------------------

            %Initial conditions Euler-Maruyama, we use lna to get the initial condition:
            
            %Initial condition:
            initcond = [0;paramaux(3)];
            
            [y0singlecell,CovMat0singlecell] = lna_v10(initcond,t0lna,t1lna,Mlna,zeros(2),[paramaux(1:17),zeros(1,6),paramaux(24),zeros(1,4),paramaux(29:end)],2);
                
            y0 = [y0singlecell;y0singlecell;y0singlecell];
            CovMat0 = blkdiag(CovMat0singlecell,CovMat0singlecell,CovMat0singlecell);

            L = chol(CovMat0,'lower'); %It's faster than simply chol, if CovMat is sparse.

            Z0D = repmat(y0,1,nsimulations) + L*N01;  %Z is a matrix with nsimulations columns and each column is a solution



            %DISTANCE BETWEEN EXPERIMENTAL AND SIMULATED DATA:
            %-------------------------------------------------

            %Vector Distance mutants:
            vdistmut = zeros(1,nmutants);
            
            %Matrix with fates of mutants:
            fatesmut = zeros(3,4,nmutants);

            distance = 0;
            
            mutantindex = 0;
            
            while (distance<=newEpT)&&(mutantindex<nmutants)
                
                mutantindex = mutantindex+1;
                
                mutant = mutantstofit(mutantindex);
 
                [distancemutant,fatesmutaux] = feval(distancehandle,paramaux,t0,t1,M,dt,MPC,dtPC,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3,mutant,weightsdistances(mutantindex),expdata{1,mutant});
                
                fatesmut(:,:,mutantindex) = fatesmutaux;
                
                vdistmut(mutantindex) = distancemutant/weightsdistances(mutantindex);
                
                distance = distance+distancemutant;

                
            end              


            %ACCEPT THE PARTICLE OR NOT DEPENDING ON THE DISTANCE:
            %-----------------------------------------------------
            if distance<=newEpT

                %Count that a new particle has been found:
                N = N+1;
                
                %Save the data in the matrix:    
                NewData = [NewData; paramaux, vdistmut, distance,1];
                NewFates(:,:,:,N) = fatesmut;
                if rem(N,100)==0
                   save(['/Users/elenacamachoaguilar/Dropbox (Personal)/Warwick/Vulval Development/5 Cusp and Saddle Node Model/Programmes v30 v2/',namenewdata,'_',num2str(jobnum),'_',num2str(N/100)],'newEpT','NewData','NewFates','streamnum','substreamnum','i','icomp')
 
                end    
                
            end
               
        end


    end
    
    eltime=toc;
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               SAMPLER T>0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
if T<Tmax
 tic
    %Update threshold number:
    %------------------------
    
    T = T+1;
    N = 0;
    i = 0;
    icomp = 0;
    
    %Set new threshold so that the threshold for each mutant is EpT(T):
    %------------------------------------------------------------------
    newEpT = EpT(T);
    
    %New cell in the particles matrix:
    %---------------------------------
    NewData = [];
      NewFates = zeros(3,4,nmutants,Numberofparticleseachjob); %NewFates(i,j,k,m) is the probability of cell i taking fate j for experiment k when using particle m.
  
    
%     %Compute the covariance matrix corresponding to each of the particles
%     %of the previous step to use a local multivariate normal perturbation
%     %kernel
%     %--------------------------------------------------------------------
% 
%     CovarianceMatricesParticles = feval(covarianceOLCMhandle,PreviousData,parfitnumbers,nparfit,Nmax,EpT(T));
     
    %Weights of each particle:
    %-------------------------
    w0 = PreviousData(:,end);
    
    %Find the particles:
    %-------------------
    while (N<Numberofparticleseachjob)&&(i<MaxStep)

        if rem(i,1000)==0
            disp('---------------------------------------------------')
            disp('---------------------------------------------------')
            disp(strcat('            Threshold: ',num2str(EpT(T))));
            disp(strcat('            Job number: ',num2str(jobnum)));
            disp(strcat('            Step: ',num2str(i)));
            disp(strcat('            Step Comp: ',num2str(icomp)));
            disp(strcat('            Accepted particles: ',num2str(N)));
            disp('---------------------------------------------------')
            disp('---------------------------------------------------')
        end


        %Update steps:
        %-------------

        i = i+1;
        

        %SAMPLE THE NEXT PARTICLE:
        %-------------------------
 
        %Computing the next particle:
        %-----------------------------
        particlenum = randsample(Nmax,1,true,w0);

        paramaux = PreviousData(particlenum,1:npar); %Sample a particle from the previous step and take the parameters that we are fitting

        newcomponentsparamaux = mvnrnd(paramaux(parfitnumbers),CovarianceMatricesParticles(:,:,particlenum)); % Use multivariate normal perturbation kernel to find new particle.


        %Check if the new components are valid:
        %-------------------------------------
        valuepriors = zeros(1,nparfit);

        parami = 0;  %Set the index of the parameter for which we will do for loop

        priorOK = 1;  %It will keep track that all the prior evaluated at every parameter is not 0

         while (parami<nparfit)&&(priorOK)  %Only if priorOK is not 0 and parami<npar the while loop will stop
             
%             minimumM = (2*sqrt(2)/(sqrt(2)-1))*(atanh(1-10^(-12))/paramaux(2));

            parami = parami+1;

            paramifit = parfitnumbers(parami);

            provparvalue = newcomponentsparamaux(parami); %Take the parami component of the new candidate vector that contains the new values of the parameters that we are fitting.
            
            priorOK = feval(evalpriorshandle,provparvalue,paramifit,minimumM); %We will have to change this so that it stops if we found one parameter which prior is 0
            
            valuepriors(parami) = priorOK;
            
            if priorOK

            priorOK = feval(constraintshandle,provparvalue,paramifit,minimumM); %Check that is a valid value

            end

            paramaux(paramifit) = provparvalue;

            

         end
         
        %Compute cPC and DiscriminantPC:
        %-------------------------------
        aPC = paramaux(10);
        bPC = paramaux(11);
        cPC = paramaux(12);
        DiscriminantPC = (8*aPC^3+27*bPC^2);
        
        %Check cPC<0 and DiscriminantPC<0
        %-------------------------------
        if priorOK
            
            if (cPC<0)&&(DiscriminantPC<0)
                priorOK = 1;
                
            else
                priorOK = 0;
                
            end
            
        end

        %Compute n1x as function of n1y
        %------------------------------
        
        paramaux(16) = -sqrt(1-paramaux(17)^2);
        
        
        
        %Check that s1=s2=s3 and l1=l2=l3:
        %---------------------------------
        
        paramaux(18:20) = paramaux(18)*ones(1,3);
        paramaux(21:23) = paramaux(21)*ones(1,3);
        
        
        %Check that EGF and Notch decay is big enough
        %--------------------------------------------
%         disp('----')
%         priorOK
        if priorOK
            
            if (paramaux(18)*exp(-paramaux(25)*(t1PC-t1))<0.001)&&(paramaux(21)*exp(-paramaux(26)*(t1PC-t1))<0.001)
                priorOK = 1;
                paramaux(27)=paramaux(25);         
            else
                priorOK = 0;
                
            end
            
        end
        
        
        %Write m31 as functions of the other ms:
        %--------------------------------
        paramaux(8) = sqrt(1-paramaux(4)^2-paramaux(6)^2);

        
        
        
%         %Check that H1=H2:
%         %-----------------
%         paramaux(13) = paramaux(11);

        
        %Check that the parameters fulfill all the constraints:
        %------------------------------------------------------
        
        if priorOK
            
            priorOK = feval(relationsconstraintshandle1,paramaux,minimumM,cPC);
            
        end
        
        
        if priorOK
            % Find m22 and m32 as functions of the other ms:
            %-----------------------------------------------
            [paramaux(7),paramaux(9),flag] = findm22m32(paramaux(4),paramaux(5),paramaux(6),paramaux(8));
        
           priorOK=flag;
           
        end
       
        %Check that the parameters fulfill all the constraints:
        %------------------------------------------------------
        
        if priorOK
            
            priorOK = feval(relationsconstraintshandle2,paramaux);
            
        end
        

        if priorOK
            icomp = icomp+1;
            
            %Post Competence parameters:
            %---------------------------
             saddlefate3 = paramaux(3)-sqrt(-cPC);
             equilibriax = equilibria_fates_1_2(bPC,DiscriminantPC);
            
            %COMPUTE THE NOISE:
            %------------------
            
            D = paramaux(end);
             
            noisecov=sqrt(2*dt*D);
             
%             Z1D = randn(nsimulations,(M))*noisecov;  %matrix of independent random deviations...
%             Z2D = randn(nsimulations,(M))*noisecov;  %Each row contains the random deviations of one
%             Z3D = randn(nsimulations,(M))*noisecov;  %simulation. Each column corresponds to a time step.
%             Z4D = randn(nsimulations,(M))*noisecov;
%             Z5D = randn(nsimulations,(M))*noisecov;
%             Z6D = randn(nsimulations,(M))*noisecov;
%             
%             Z1DPC = randn(nsimulations,(MPC))*noisecov;  %matrix of independent random deviations...
%             Z2DPC = randn(nsimulations,(MPC))*noisecov;  %Each row contains the random deviations of one
%             Z3DPC = randn(nsimulations,(MPC))*noisecov;  %simulation. Each column corresponds to a time step.
%             Z4DPC = randn(nsimulations,(MPC))*noisecov;
%             Z5DPC = randn(nsimulations,(MPC))*noisecov;
%             Z6DPC = randn(nsimulations,(MPC))*noisecov;

%Vectorial:
            Z1D = randn((M),nsimulations)*noisecov;  %matrix of independent random deviations...
            Z2D = randn((M),nsimulations)*noisecov;  %Each row corresponds to a time step.
            Z3D = randn((M),nsimulations)*noisecov;  %Each column contains the random deviations of one simulation.
            Z4D = randn((M),nsimulations)*noisecov;
            Z5D = randn((M),nsimulations)*noisecov;
            Z6D = randn((M),nsimulations)*noisecov;
            
            Z1DPC = randn((MPC),nsimulations)*noisecov;  %matrix of independent random deviations...
            Z2DPC = randn((MPC),nsimulations)*noisecov;  %Each row contains the random deviations of one
            Z3DPC = randn((MPC),nsimulations)*noisecov;  %simulation. Each column corresponds to a time step.
            Z4DPC = randn((MPC),nsimulations)*noisecov;
            Z5DPC = randn((MPC),nsimulations)*noisecov;
            Z6DPC = randn((MPC),nsimulations)*noisecov;
            
            N01 = randn(6,nsimulations);
             

            % REAL INITIAL CONDITION:
            %------------------------

            %Initial conditions Euler-Maruyama, we use lna to get the initial condition:

            %Initial condition:
            initcond = [0;paramaux(3)];
            
            [y0singlecell,CovMat0singlecell] = lna_v10(initcond,t0lna,t1lna,Mlna,zeros(2),[paramaux(1:17),zeros(1,6),paramaux(24),zeros(1,4),paramaux(29:end)],2);
                
            y0 = [y0singlecell;y0singlecell;y0singlecell];
            CovMat0 = blkdiag(CovMat0singlecell,CovMat0singlecell,CovMat0singlecell);

            L = chol(CovMat0,'lower'); %It's faster than simply chol, if CovMat is sparse.

            Z0D = repmat(y0,1,nsimulations) + L*N01;  %Z is a matrix with nsimulations columns and each column is a solution



            %DISTANCE BETWEEN EXPERIMENTAL AND SIMULATED DATA:
            %-------------------------------------------------

            %Vector Distance mutants:
            vdistmut = zeros(1,nmutants);
            
            %Matrix with fates of mutants:
            fatesmut = zeros(3,4,nmutants);

            distance = 0;
            
            mutantindex = 0;
            
            while (distance<=newEpT)&&(mutantindex<nmutants)
                
                mutantindex = mutantindex+1;
                
                mutant = mutantstofit(mutantindex);
                
                [distancemutant,fatesmutaux] = feval(distancehandle,paramaux,t0,t1,M,dt,MPC,dtPC,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3,mutant,weightsdistances(mutantindex),expdata{1,mutant});
                
                fatesmut(:,:,mutantindex) = fatesmutaux;
                
                vdistmut(mutantindex) = distancemutant/weightsdistances(mutantindex);
                
                distance = distance+distancemutant;
                
            end   


            %ACCEPT THE PARTCILE OR NOT DEPENDING ON THE DISTANCE:
            %-----------------------------------------------------

            if distance<=newEpT

                %Count that a new particle has been found:
                N = N+1;
  
                %Compute the weight of the particle:
                denominator = 0;

                for partic=1:Nmax
                   
                        denominator = denominator+w0(partic)*feval(evalkernelhandle,newcomponentsparamaux,PreviousData(partic,parfitnumbers),CovarianceMatricesParticles(:,:,partic));

                end

                w1 = prod(valuepriors)/denominator; %We multiply the priors of each component of the parameter vector to compute the value of the prior of the parameter vector

                NewData = [NewData; paramaux, vdistmut, distance, w1];
                NewFates(:,:,:,N) = fatesmut;    
                if rem(N,100)==0
                   save(['/Users/elenacamachoaguilar/Dropbox (Personal)/Warwick/Vulval Development/5 Cusp and Saddle Node Model/Programmes v30 v2/',namenewdata,'_',num2str(jobnum),'_',num2str(N/100)],'newEpT','NewData','NewFates','streamnum','substreamnum','i','icomp')
 
                end
            end
               
        end


    end
    eltime=toc;

  end
  
save(['/Users/elenacamachoaguilar/Dropbox (Personal)/Warwick/Vulval Development/5 Cusp and Saddle Node Model/Programmes v30/New_Final_Vectorial_Version/',namenewdata,'_',num2str(jobnum)],'newEpT','NewData','NewFates','streamnum','substreamnum','eltime','i','icomp')
%save([namenewdata,'_',num2str(jobnum)],'newEpT','NewData','NewFates','streamnum','substreamnum','eltime','i')
% save(['/dascratch/ec47/',namenewdata,'_',num2str(jobnum)])

reset(RandStream.getGlobalStream)


% fileID = fopen(['/scratch/ec47/FilesThatRun.txt'],'a');
% fprintf(fileID,num2str(jobnum));
% fprintf(fileID,'\n');
% fclose(fileID);
    
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%20/01/16

%Here we save the distance from each mutant. (09/12/16)

%Here we save the fates for each mutant. (20/01/17)

%Here we changed the distance function. (25/05/17)

%Here we write the parameter d as a linear function of delta1, and check that H1=H2. (04/07/17)

%Here we change H1=H2. (04/07/17)

%Here we don't set H1=H2, and use model = 'Vulval_Development_v5'. (03/08/17)

%Here we use model = 'Vulval_Development_v7', we have added competence and
%post competence time, and constraints in c post competence and
%Discriminant post competence, as well as try to fit a data set that has
%been simulated with a parameter that we know (27/09/17)

%v10 Mutants: We have added the constraints up to Mutant 4

%v12 Added more parameters to the fitting.

%v13 Added weights to the distances. 

%Modelv5_v1 Change it to use Model version 5 (competence time and weights).
%(29/11/17)

%Modelv6_v1 Change it to use Model version 6 (competence time and signal decay and weights).
%(04/12/17)

%Modelv6_v2 Change it to use Model version 6 and added the constraint that
%EGF is perpendicular to Notch. Therefore I have to compute some ms as
%functions of the known ms, and then check the relations constraints again.
%(04/12/17)

%Modelv7_v1 Change it to use Model version 7 and added the constraint that
%EGF is perpendicular to Notch. Therefore I have to compute some ms as
%functions of the known ms, and then check the relations constraints again.
%In this version, C*A>0, therefore the cusp is in the 1,2,3 region.
%(18/1/18)

%Modelv9_v1 Change it to use Model version 9 (with EGF response as sigmoidal function)
% and made dt step smaller. Also changed the decay function so that it's equal to 
% one at the beginning of the stage. Also changed the threshold for the
% decay function.
%(07/11/18)

%Changed path to save the data for the cluster

%Changed valuepriors before the if statement when T>0 to save the actual prior value

%%Modelv10_v1 Change it to use Model version 10, which is model 9 but the Notch
% sigmoidal is taken as in Siggia,2012. We fit the parameters of this
% sigmoidal function (so, in total, 16 parameters). Also fixed npar in
% ABC_SMC function, which was set to 37 instead of 36,so when T>1, the noise was equal to
%the distance to mutant 1!!
%(22/04/20)  

%ABC_SMC_Algorithm_OLCM_Template.m (c) by Elena Camacho Aguilar

%ABC_SMC_Algorithm_OLCM_Template.m is licensed under a 
%Creative Commons Attribution-ShareAlike 3.0 Unported License.

%You should have received a copy of the license along with this
%work.  If not, see <http://creativecommons.org/licenses/by-sa/3.0/>.
