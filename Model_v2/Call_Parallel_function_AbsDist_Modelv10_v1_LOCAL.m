function struc = Call_Parallel_function_AbsDist_Modelv10_v1_LOCAL(Mut,nbatch)

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEW PRIORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Mut==1

    EpTvector = [5,3.3,2.8,2.32,1.9967,1.765,1.5751,1.3398,1.156,0.979,0.8336,0.6035,0.5242,0.4619]% Maximum threshold is (6*7+4*2)/9 = 50/9 = 5.5 =approx 5

    names = cell(1,1);
    names{1} = {'','_Eps3p3','_Eps2p8','_Eps2p32','_Eps1p9967','_Eps1p765','_Eps1p5751','_Eps1p3398','_Eps1p156','_Eps0p979','_Eps0p8336','_Eps0p6035','_Eps0p5242','_Eps0p4619'};

    loopept=14

    namedata = {'20_07_27_9Mut_SENSITIVITYNEWPRIORS_ABCMv10v1_20000part_150sim_22par_Eps5_Eps3p3_Eps2p8_Eps2p32_Eps1p9967_Eps1p765_Eps1p5751_Eps1p3398_Eps1p156_Eps0p979_Eps0p8336_Eps0p6035_Eps0p5242'};
    
    pathtodata = '/Users/elenacamachoaguilar/Dropbox (Personal)/Warwick/Vulval Development/5 Cusp and Saddle Node Model/Programmes v30/New_Final_Vectorial_Version/';

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

