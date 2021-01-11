function   sol=equilibria_fates_1_2(b,D)

if D==0
    
    S = nthroot(b,3);
    
   if b>0
       
       sol = [-S,S/2];
       
   else
       
       sol = [S/2,-S];   
       
   end

    
else %D<0
    
    Disc = D/27;
    
    u = ((-b+sqrt(Disc))^(1/3))/2;
    
    A = real(u);
    
    B = sqrt(3)*imag(u);
    
    sol = sort([2*A,-A-B,-A+B]);
    
end