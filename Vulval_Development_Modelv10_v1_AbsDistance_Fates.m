function [distance,fates]=Vulval_Development_Modelv10_v1_AbsDistance_Fates(p,t0,t1,M,dt,MPC,dtPC,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3,mutantnumber,weightmutant,expdata)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                        %
%                    v5 cusp-saddlenode                                  %
%                                                                        %
%                                                                        %
%   Given a parameter set and an initial condition it computes the       %
%   fates of the simulations.
%                                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



if mutantnumber==1

    %WT:
    %---
    tAC = 10;
    
    fates = simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,p,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3);

    distance = weightmutant*sum(reshape(abs(fates-expdata),1,12));
    
elseif mutantnumber==2
    
    %No EGF receptors in P5p
    %-----------------------
    p(19) = 0;    %no EGF signal in P5p (s2=0)
    
    tAC = 10;
    
    fates = simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,p,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3);

    distance = weightmutant*sum(reshape(abs(fates-expdata),1,12));
    
elseif mutantnumber==3
    
    %Half dose lin 3 (Half EGF)
    %--------------------------
    p(18:20)=p(18:20)/2;
    
    tAC = 10;
    
    fates = simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,p,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3);

    distance = weightmutant*sum(reshape(abs(fates-expdata),1,12));
    
elseif mutantnumber==4
    
    %Half dose lin 12
    %----------------
    p(21:23)=p(21:23)/2;
    
    tAC = 10;
    
    fates = simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,p,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3);

    distance = weightmutant*sum(reshape(abs(fates-expdata),1,12));
    
elseif mutantnumber==5
    
    %Notch null, 2ACs
    %----------------

    %s1,s2,s3 are double 
    p(18:20) = 2*p(18:20);
    
    %no notch
    p(21:23) = zeros(1,3);
    
    tAC = 10;
    
    fates = simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,p,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3);

    distance = weightmutant*sum(reshape(abs(fates-expdata),1,12));
    
elseif mutantnumber==6
    
    %No Notch Signal: l1=l2=l3=0
    %-----------------------
    p(21:23)=zeros(1,3);
    
    tAC = 10;
    
    fates = simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,p,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3);

    distance = weightmutant*sum(reshape(abs(fates-expdata),1,12));

elseif mutantnumber==7
    
    %EGF overexpression: s1=EGFover*s1...
    %-----------------------------
    p(18:20)=p(35)*p(18:20);
    
    tAC = 10;
    
    fates = simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,p,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3);

    distance = weightmutant*sum(reshape(abs(fates-expdata),1,12));

elseif mutantnumber==8

    tAC = 0.2;
    
    fates = simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,p,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3);

    distance = weightmutant*sum(reshape(abs(fates(2:3,:)-expdata),1,8));

elseif mutantnumber==9
    
    tAC = 0.32;
    
    fates = simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,p,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3);
    
    distance = weightmutant*sum(reshape(abs(fates(2:3,:)-expdata),1,8));
    
elseif mutantnumber==10
    
    tAC = 0.44;
    
    fates = simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,p,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3);
    
    distance = weightmutant*sum(reshape(abs(fates(2:3,:)-expdata),1,8));
    
elseif mutantnumber==11
    
    tAC = 0.56;
    
    fates = simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,p,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3);
    
    distance = weightmutant*sum(reshape(abs(fates(2:3,:)-expdata),1,8));

elseif mutantnumber==12
    
    tAC = 0.68;
    
    fates = simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,p,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3);
    
    distance = weightmutant*sum(reshape(abs(fates(2:3,:)-expdata),1,8));
    
elseif mutantnumber==13
    
    tAC = 0.8;
    
    fates = simulationEulerACablation_Competence_v10_vec(t0,t1,tAC,dt,dtPC,M,MPC,p,nsimulations,Z0D,Z1D,Z2D,Z3D,Z4D,Z5D,Z6D,Z1DPC,Z2DPC,Z3DPC,Z4DPC,Z5DPC,Z6DPC,aPC,bPC,cPC,DiscriminantPC,equilibriax,saddlefate3);
    
    distance = weightmutant*sum(reshape(abs(fates(2:3,:)-expdata),1,8));

end




