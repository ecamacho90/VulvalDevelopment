function priorOK=Vulval_Development_Modelv10_v1_Relations_Constraints1(paramaux,minimumM,cPC)

priorOK = 0;




 if paramaux(8)~=0 % C not equal to 0

                            if paramaux(3)>=minimumM % All mutants: Bifurcation step function:

                            criticaly = atanh(1-10^(-12))/paramaux(2)+(paramaux(3)*(sqrt(2)-1))/(2*sqrt(2));

                            criticalC = -(criticaly-paramaux(3))^2;

                            if cPC > criticalC % All mutants: Values of H,M,c valids for the step function
%                                             disp('1')
% Not changed for version 5 yet       cbif = - paramaux(4)/cos(paramaux(5));
% 
%                                     if (2*paramaux(16)*paramaux(10)*cos(paramaux(8))+paramaux(6)>cbif) % Mutant 5: Third fate not present for Notch null, 2WT EGF for P5p
% 
                                            priorOK=1;
% 
%                                     end
                            end


                            end
end
