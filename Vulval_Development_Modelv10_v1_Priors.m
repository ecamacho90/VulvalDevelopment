function pvalue=Vulval_Development_Modelv10_v1_Priors(paramnumber,minimumM)

if paramnumber==1
    %1/tau
    %pvalue is taken from a gamma distribution obtained from the particles
    %that fit biologically the experimental data on the previous version
    %with 5000 particles
    pvalue = random('Gamma',2,2);
    
elseif paramnumber==2
    %H
    pvalue = random('Gamma',2,2);
    
elseif paramnumber==3
    %M
    %pvalue is taken from a gamma distribution obtained from the particles
    %that fit biologically the experimental data on the previous version
    %with 5000 particles
    pvalue = minimumM+random('Gamma',2,2);
        
elseif paramnumber==4   %Modelv7_v2
    %m11
    pvalue = random('Uniform',0,1);
            
elseif paramnumber==5   %Modelv7_v2
    %m12
    pvalue = random('Uniform',0,1);
    
elseif paramnumber==6
    %m21
    pvalue = random('Uniform',-1,1);
        
elseif paramnumber==7
    %m22
    pvalue = random('Uniform',-1,1);
            
elseif paramnumber==8
    %m31
    pvalue = random('Uniform',0,1);

elseif paramnumber==9
    %m32
    pvalue = random('Uniform',0,1);
    
elseif paramnumber==10
    %q1
    pvalue = random('Normal',0,2);
    
elseif paramnumber==11
    %q2
    pvalue = random('Normal',0,2);
    
elseif paramnumber==12
    %q3
    pvalue = -random('Gamma',2,2);

elseif paramnumber==13
    %gamma
    pvalue = random('Uniform',0,1);

elseif paramnumber==14
    %n0
    pvalue = random('Normal',0,1);
    
elseif paramnumber==15
    %n1
    pvalue = random('Gamma',2,2);

elseif paramnumber==16
    %n1x
    pvalue = random('Uniform',-1,0);

elseif paramnumber==17
    %n1y
    pvalue = random('Uniform',-1,1);
    
elseif paramnumber==18
    %s1=s2=s3
    pvalue = random('Gamma',2,2);

elseif paramnumber==19
    %s2

elseif paramnumber==20
    %s3

elseif paramnumber==21
    %l1=l2=l3
    %pvalue is taken from an exponential distribution:
    pvalue = random('Gamma',2,2);

elseif paramnumber==22
    %l2

elseif paramnumber==23
    %l3

elseif paramnumber==24
    %alpha
    %pvalue is taken from an exponential distribution:
    pvalue = random('Exponential',1);
    
elseif paramnumber==25
    %lambdaE
    %pvalue is taken from an exponential distribution:
pvalue = random('Gamma',3/2,4);

elseif paramnumber==26
    %lambdaN
    %pvalue is taken from an exponential distribution:
pvalue = random('Gamma',3/2,4);

elseif paramnumber==27
    %xiE
    %pvalue is taken from an exponential distribution:
pvalue = random('Gamma',3/2,4);

elseif paramnumber==28
    %xiN
    %pvalue is taken from an exponential distribution:
pvalue = random('Gamma',3/2,4);

elseif paramnumber==29
    %L
    %pvalue is taken from an exponential distribution:
    pvalue = random('Uniform',0,1);
    
elseif paramnumber==30
    %H3
    pvalue = random('Gamma',3/2,4);
    
elseif paramnumber==31
    %M3
    %pvalue is taken from an exponential distribution:
    pvalue = random('Normal',0,2);
    
elseif paramnumber==32
    %C
    %pvalue is taken from an exponential distribution:
    
elseif paramnumber==33
    %tAC
    %pvalue is taken from an exponential distribution:
    
elseif paramnumber==34
    %t1
    %pvalue is taken from an exponential distribution:
    
elseif paramnumber==35
    %EGFover
    %pvalue is taken from an exponential distribution:
    pvalue = random('Gamma',2,2);
    
elseif paramnumber==36
    %D
    pvalue = random('Uniform',0.05,1);
    
end

%New EGF prior
