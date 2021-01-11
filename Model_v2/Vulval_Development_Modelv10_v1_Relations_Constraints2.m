function priorOK=Vulval_Development_Modelv10_v1_Relations_Constraints2(paramaux)

priorOK = 0;

A = (paramaux(8)*paramaux(7)-paramaux(6)*paramaux(9));
B = (paramaux(4)*paramaux(9)-paramaux(8)*paramaux(5));
C = (paramaux(5)*paramaux(6)-paramaux(4)*paramaux(7));
D = A*paramaux(11)+B*paramaux(12)+C*paramaux(13);


 if A ~= 0 % A not equal to 0
%       disp('1')
     if C ~= 0 % C not equal to 0
         
         % if paramaux(8)>0 (already in ABC script)
             
%              if paramaux(9)>0 (already in findm22m32 script)
%         disp('2')
%          if
%          (paramaux(1)*paramaux(2)+paramaux(3)*paramaux(4)+paramaux(5)*paramaux(6))          == 0 % EGF axis perpendicular to Notch (Already taken into
%          account for computing m32 and m22)
%             disp('3')

%              if paramaux(12) < 0 %Origin with fate 3 present (Already in
%              the ABC script)
                 
%                  if 8*paramaux(11)^3+27*paramaux(12)^2 < 0 %Origin with
%                  fates 1 and 2 present (already in the ABC script)

%   Version 1:
                     if C*A>0 % Cusp in 1,2,3 region
%                         disp('4') 
                         if D*A<0 %B1 intersection B2 intersection the plane is not empty
%                              disp('5')

                            if ((paramaux(8)*paramaux(18)+paramaux(12))>0) %Constraint for Mutant 6 (P6.p)
                                
                                if ((paramaux(8)*2*paramaux(18)*paramaux(13)+paramaux(12))>0) %Constraint for Mutant 5 (P5.p)
                                    
                                     if B==0
                                         priorOK = 1;

                                     elseif B~=0
%                                          AA = -8*B^3;
%                                          BB = 24*D*B^2+27*A^3;
%                                          CC = -24*B*D^2;
%                                          DD = 8*D^3;
%                                          Q = (3*AA*CC-BB^2)/(9*AA^2);
%                                          R = (9*AA*BB*CC-27*AA^2*DD-2*BB^3)/(54*AA^3);
%                                          DiscPi = Q^3+R^2;
% 
%                                          if DiscPi<=0
        %                                  disp('6')

        % Not changed for version 5 yet       cbif = - paramaux(4)/cos(paramaux(5));
        % 
        %                                     if (2*paramaux(16)*paramaux(10)*cos(paramaux(8))+paramaux(6)>cbif) % Mutant 5: Third fate not present for Notch null, 2WT EGF for P5p
        % 
                                                        priorOK=1;
        % 
                                         end


                                     end
                                end
                            end
                         end
                     end
%                  end
%              end
%         end
%              end
%          end
     end
end

%v10: Competence time added, so we don't need constraints for mutants 6 and
%7