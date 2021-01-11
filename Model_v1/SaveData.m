clear all

struc = Call_Parallel_function_AbsDist_Modelv10_v1(1,1);    
pathtodata = '';

%Simulations 1 to 35 were ran in the cluster
% Simulations 36 to 100 were ran in the computer
%%

% Number of parallel computations:
jobmax = struc.jobmax;

% Vector of jobs
jobsvector = struc.jobsvector;

EpTvector = struc.EpTvector;

names = struc.names;

%for loopept=4:5
loopept=struc.loopept;

% Number of particles:
Nmax = struc.Nmax;

namedata = struc.namedata;

%Model name:
model = struc.model;%'Vulval_Development_v11';

% ABC version:
ABCversion = struc.ABCversion; %'ABC_SMC_Algorithm_OLCM_Template_v11_Fates_AbsDistance'

% Weights for the distances:
weightsdistances = struc.weightsdistances;

%Path to data 
pathtodata = struc.pathtodata;

% Variables

%-----------------------
% To change:
%-----------------------

% Mutants to fit:
mutantstofit = struc.mutantstofit;

% Parameters to fit:
parfitnumbers = struc.parfitnumbers;

% New threshold: (CHANGE)
EpTsmall = EpTvector(loopept);

% Chain of thresholds, including the new: (CHANGE)
EpT = EpTvector(1:loopept);

% Name of the last data set: (CHANGE)
namedata = namedata;

% Name of the new data set:(CHANGE)
namenewdata = [namedata,names{loopept}];


%-----------------------
% Not to change:
%-----------------------

%Random number generator stream:
nthresholdsEpT = length(EpT);
streamnum = 2*nthresholdsEpT-1; 

%Number of parameters to fit:
nparfit = length(parfitnumbers);

% Number of thresholds computed so far:
T = length(EpT)-1;

% Number of the threshold that we will use: 
Tmax = length(EpT);

% Number of particles in each parallel computation:
Numberofparticleseachjob = Nmax/jobmax;

%Covariance matrices OLCM function handle:
covarianceOLCMhandle = str2func(strcat(model,'_CovarianceMatricesOLCM'));



% Check if we have previous data

if T==0
    PreviousData=[];
    CovarianceMatricesParticles = [];
    
else
    
    load([pathtodata,namedata],'ParticlesMatrix')
    PreviousData = ParticlesMatrix{1,T};
    
    %Compute the covariance matrix corresponding to each of the particles
    %of the previous step to use a local multivariate normal perturbation
    %kernel
    %--------------------------------------------------------------------

    CovarianceMatricesParticles = feval(covarianceOLCMhandle,PreviousData,parfitnumbers,nparfit,Nmax,EpT(T+1));
    
    
end

%% Save the data:
% 
load(namedata,'FatesMatrix','itotal','icomptotal','eltimetotal')
ParticlesMatrixaux = [];
FatesMatrixaux = [];
streamJobs = cell(1,jobmax);
substreamJobs = cell(1,jobmax);
itotalaux = 0;
icomptotalaux = 0;
eltimetotalaux = 0;

for jobnum = 1:jobmax
    
    load([namenewdata,'_',num2str(jobnum)],'NewData','NewFates','streamnum','substreamnum','i','icomp','eltime')
    
    ParticlesMatrixaux = [ParticlesMatrixaux;NewData];
    FatesMatrixaux = cat(4,FatesMatrixaux,NewFates);
    streamJobs{jobnum} = streamnum;
    substreamJobs{jobnum} = substreamnum;
    itotalaux = itotalaux+i;
    icomptotalaux = icomptotalaux+icomp;
    eltimetotalaux = eltimetotalaux+eltime;
    
    
end

%Normalise the weights:
ParticlesMatrixaux(:,end) = ParticlesMatrixaux(:,end)/sum(ParticlesMatrixaux(:,end));

%%
%Save it in the cell of results:
if T==0

    ParticlesMatrix{1,1} = ParticlesMatrixaux;
    FatesMatrix{1,1} = FatesMatrixaux;
    itotal{1,1} = itotalaux;
    icomptotal{1,1} = icomptotalaux;
    eltimetotal{1,1} = eltimetotalaux;

elseif T>0

    ParticlesMatrix{1,T+1} = ParticlesMatrixaux;
    FatesMatrix{1,T+1} = FatesMatrixaux;
    itotal{T+1} = itotalaux;
    icomptotal{T+1} = icomptotalaux;
    eltimetotal{T+1} = eltimetotalaux;
    
end

clear('ans')

save(namenewdata)
% 
disp('Data Saved')