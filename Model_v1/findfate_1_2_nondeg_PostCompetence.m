function [fate,ysol]=findfate_1_2_nondeg_PostCompetence(y,c,b,a,H,M,saddley,equilibriax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%       Given the coordinates of the end point of the trajectory of the
%       cell, it looks for the attractor towards which the trajectory will
%       converge in case we take that point as initial condition.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% % % % To try as an example:
% % % y=[2,0.3];
% % % c=1
% % % b=-1
% % % a=-1.5
% % % D=8*a^3+27*b^2
% % % H=25
% % % M=5


x0=y(1);

%We choose how close to the attractors we need the trajectory to be in
%order to score the fate.

if b<0 %We know that the saddle is closer to fate 1 and the stable manifold is bended towards fate 1.
    
    mindistancef2 = (equilibriax(3)-equilibriax(2))*0.9;
    
    if c<=0 
        mindistancef1 = min(saddley,equilibriax(2)-equilibriax(1))*1.e-3;  %Minimum between the distance between the saddle of fate 3 and the distance between the equilibria 1 and the saddle
        
    else 
        
        mindistancef1 = (equilibriax(2)-equilibriax(1))*1.e-3;  %There is no saddle f3 so we take the distance between the equilibria 1 and the saddle as minimum distance
    end
    
elseif b>0 %We know the same but for fate 2
    
    mindistancef1 = (equilibriax(2)-equilibriax(1))*0.9;
    
    if c<=0 
        mindistancef2 = min(saddley,equilibriax(3)-equilibriax(2))*1.e-3;  %Minimum between the distance between the saddle of fate 3 and the distance between the equilibria 1 and the saddle
        
    else 
        
        mindistancef2 = (equilibriax(3)-equilibriax(2))*1.e-3;  %There is no saddle f3 so we take the distance between the equilibria 1 and the saddle as minimum distance
    end
    
elseif b==0  %The stable manifold of the saddle is perpendicular to x=0, so we now the basins of attraction
    
    mindistancef1 = (equilibriax(2)-equilibriax(1))*0.9;  %There is no saddle f3 so we take the distance between the equilibria 1 and the saddle as minimum distance
    mindistancef2 = (equilibriax(3)-equilibriax(2))*0.9;  %There is no saddle f3 so we take the distance between the equilibria 2 and the saddle as minimum distance
    
    
end
    




if norm(y-[equilibriax(1);0])-mindistancef1<=0

    fate = 1;
    
elseif norm(y-[equilibriax(3);0])-mindistancef2<=0

    fate = 2;

else

    
    options = odeset('Events',@(t,y) myEventsFcnNonDeg(t,y,equilibriax,mindistancef1,mindistancef2));

    p = [c,b,a,H,M];

    tspan = 0:10:10000;

    [tsol,ysol,te,ye,ie] = ode15s(@(t,y) singlecell_cusp_and_saddlenode_model(t,y,p),tspan,y,options);

    if isempty(ie) || not(length(ie)==1)

        fate = 4;

    elseif ie==1

        fate = 1;

    else

        fate = 2;

    end


end


