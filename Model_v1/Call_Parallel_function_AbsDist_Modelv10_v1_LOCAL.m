function struc = Call_Parallel_function_AbsDist_Modelv10_v1_LOCAL(Mut,nbatch)

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW PRIORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Mut==1

    EpTvector = [5,3.3,2.68,2.147,1.836,1.625,1.431,1.25,1.05,0.91,0.72,0.583,0.497,0.4286]% Maximum threshold is (6*7+4*2)/9 = 50/9 = 5.5 =approx 5

    names = cell(1,1);
    names{1} = {'','_Eps3p3','_Eps2p68','_Eps2p147','_Eps1p836','_Eps1p625','_Eps1p431','_Eps1p25','_Eps1p05','_Eps0p91','_Eps0p72','_Eps0p583','_Eps0p497','_Eps0p4286'};

    loopept=14

    namedata = {'20_07_08_9Mut_SENSITIVITYNEWPRIORS_ABCMv10v2_20000part_150sim_22par_Eps5_Eps3p3_Eps2p68_Eps2p147_Eps1p836_Eps1p625_Eps1p431_Eps1p25_Eps1p05_Eps0p91_Eps0p72_Eps0p583_Eps0p497'};
    
    pathtodata = 'pathtodata';

%     pathtodata ='/Users/elenacamachoaguilar/Dropbox (Personal)/Warwick/Vulval Development/5 Cusp and Saddle Node Model/Programmes v30 v2/';
    weightsdistances = ones(1,9)/9;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Number of particles:
    Nmax = 20000;

    % Number of jobs max:
%     jobmax = 5000;
%     jobmax = 2000;
    jobmax = 1000;

    % N jobs per batch:
    njobsperbatch = 10;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Mutants to fit:
    mutantstofit = [1:8,10]%;13;

    % Parameters to fit:
    parfitnumbers = [1,3,4,5,6,10,11,12,13,14,15,17,18,21,24,25,26,29,30,31,35,36];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Model:
    model = 'Vulval_Development_Modelv10_v1'; % !!! Change in ABC_SMC_Algorithm... too !!!

    % ABC version:
    ABCversion = 'ABC_SMC_Algorithm_OLCM_Template_Modelv10_v1_Fates_AbsDistance_L';

    % Vector of jobs for the batch:
    jobsvector = ((nbatch-1)*njobsperbatch + 1): njobsperbatch*nbatch;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    struc = struct('pathtodata',pathtodata,'EpTvector',EpTvector,'names',names,'loopept',loopept,'Nmax',Nmax,'namedata',namedata,'mutantstofit',mutantstofit,'parfitnumbers',parfitnumbers,'model',model,'ABCversion',ABCversion,'jobsvector',jobsvector,'jobmax',jobmax,'weightsdistances',weightsdistances);
    
end

