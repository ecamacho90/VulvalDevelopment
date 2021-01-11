function sol=computefates_PostCompetence(y,H,M,a,b,c,D,saddlef3y,equilibriax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                         %
%                           cusp-saddlenode                               %
%                                                                         %
%  This programme assumes that the cells have moved in a potential        %
%  landscape without signals after competence. This landscape is fixed    %
%  for all the cells, and must have the three attractors present (c>0 and %
%  D<0).                                                                  %
%                                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Coordinates in the bifurcation set:
%-----------------------------------

%c=p(6)*cos(p(5))+p(4); % All the cells lay on the same landscape after competence.
                       % The landscape corresponds to the one with no
%b=p(7);                % signals.

%a=-p(6)*sin(p(5));

%D=(8*a^3+27*b^2);

% The landscape needs to have the three fates, therefore:
        % c<0 (Third fate present)
        % D<0 (First and second fates present).

%Fates:
%------

sol=zeros(3,4);

    %sol(i,j) is 1 if cell i has taken fate j. If j=4 it means that the
    %fate could not be determined.
    
    
%Computation of the fates:
%-------------------------

for cell=1:3
   
    i=2*cell-1;
    j=2*cell;
    
    notfoundfate=1;
    
    
      % We assume that c<0 and D<0:
      % ---------------------------
      
        %HERE WE COMPUTE THE SADDLE IN X=0 and check if the point is in the
        %basin of attraction of fate 3.
        

        if y(j)>saddlef3y   %Check if it's in basin of attraction of fate 3 when saddle+eq
            
            sol(cell,3)=1;
            
        elseif (y(j)==saddlef3y)  %Check if saddle+eq but the point is in the stable manifold of the saddle OR degenerate and in basin of attraction of 3
            
            sol(cell,4)=1;
            
        elseif abs(D)<1.e-10   %Degenerate in 1-2
            
        fate = findfate_1_2_deg_PostCompetence(y(i:j),c,b,a,H,M,saddlef3y,equilibriax);
        
        sol(cell,fate) = 1;
        
            
        elseif D<0  %There are 2 equilibria with y=0 and a saddle

        fate = findfate_1_2_nondeg_PostCompetence(y(i:j),c,b,a,H,M,saddlef3y,equilibriax);
        
        sol(cell,fate) = 1;
            
        else
            
            sol(cell,4)=1;
            
            disp('Computefates_PostCompetence: Could not find any fate')
            
        end
        
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Elena 27/08/17
% We focus only in the cases where there are three equilibria, since the
% cells have a post-competence period in which they converge to one of the
% three fates.