%% Parallel_call_function_5000
% This script calls ABC_Algorithm in a parallel way so that the particles
% are computed in parallel subsets.
%% For loop
function Parallel_function_AbsDist_Modelv10_v1_LOCAL()
%%
%Here we fit a data set that has been simulated with a parameter vector
%that we choose.
% functionname = mfilename();
% 
% % Details for batch
% struc = Call_Parallel_function_AbsDist_Modelv10_v1(1,str2double(functionname((end-2):end)));

%%
clear all

struc = Call_Parallel_function_AbsDist_Modelv10_v1_LOCAL(1,100);
% struc = Call_Parallel_function_AbsDist_Modelv10_v1(1,1);


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



% for ii=[40,1]%100:-1:1
    ii=1
struc = Call_Parallel_function_AbsDist_Modelv10_v1_LOCAL(1,ii);
% struc = Call_Parallel_function_AbsDist_Modelv10_v1(1,1);

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



% % Check if we have previous data
% 
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
%%
% Create job in cluster:

% p = parcluster('local');
% Jobs = createJob(p);
% vecjobnum = [87,88,89,90,172,173,174,175,176,177,178,179,180,181,182,183,184,185,186,306,307,308,309,310,311,312,313,314];
% Loop to call jobs

tic
parfor jobnumaux = 1:jobsvector
    jobnum = vecjobnum(jobnumaux);
    ABC_SMC_Algorithm_OLCM_Template_Modelv10_v1_Fates_AbsDistance_L(mutantstofit,parfitnumbers,EpT,T,Tmax,Nmax,Numberofparticleseachjob,PreviousData,namenewdata,CovarianceMatricesParticles,weightsdistances,streamnum,jobnum)
 
%    createTask(Jobs,str2func(ABCversion),0,{mutantstofit,parfitnumbers,EpT,T,Tmax,Nmax,Numberofparticleseachjob,PreviousData,namenewdata,CovarianceMatricesParticles,weightsdistances,streamnum,jobnum});
%   
%    disp(['Job ',num2str(jobnum),' created']);
end
toc

% disp('Submiting job.......................')
% submit(Jobs);
% disp('Submitted!')
% disp('Waiting for it to finish............')
% 
% wait(Jobs);

disp('Finished!')

% delete(Jobs)
% clear Jobs
% 
% exit
% 

% end

disp('Finished All!')
    %%
% Save the data:
% 
load(namedata,'FatesMatrix')
ParticlesMatrixaux = [];
FatesMatrixaux = [];
streamJobs = cell(1,jobmax);
substreamJobs = cell(1,jobmax);
itotal = 0;
icomptotal = 0;

for jobnum = 1:jobmax
    
    load([namenewdata,'_',num2str(jobnum)],'NewData','NewFates','streamnum','substreamnum','i','icon')
    
    ParticlesMatrixaux = [ParticlesMatrixaux;NewData];
    FatesMatrixaux = cat(4,FatesMatrixaux,NewFates);
    streamJobs{jobnum} = streamnum;
    substreamJobs{jobnum} = substreamnum;
    itotal = itotal+i;
    icomptotal = icontotal+icon;
    
    
end

%Normalise the weights:
ParticlesMatrixaux(:,end) = ParticlesMatrixaux(:,end)/sum(ParticlesMatrixaux(:,end));
%%
%Save it in the cell of results:
if T==0

    ParticlesMatrix{1,1} = ParticlesMatrixaux;
    FatesMatrix{1,1} = FatesMatrixaux;

elseif T>0

    ParticlesMatrix{1,T+1} = ParticlesMatrixaux;
    FatesMatrix{1,T+1} = FatesMatrixaux;
    
end

clear('ans')

save(namenewdata)
% 
disp('Data Saved')
% exit