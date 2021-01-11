function parameterOK=Vulval_Development_Modelv10_v1_Constraints(paramvalue,paramnumber,minimumM)


%This function is to check that the parameters take the allowed values.

parameterOK=1;

if paramnumber==1
%tau
    if paramvalue==0
       parameterOK=0; 
    end
    
elseif paramnumber==2
    %H
elseif paramnumber==3
    %M
    if paramvalue<=minimumM
        parameterOK = 0;
    end
        
elseif paramnumber==4
    %m11
    if (paramvalue>0)||(paramvalue<-1)  %Modelv7_v1
        parameterOK = 0;
    end
    
elseif paramnumber==5
    %m12
    if paramvalue>0||(paramvalue<-1)    %Modelv7_v1
        parameterOK = 0;
    end
    
elseif paramnumber==6
    %m21
    if abs(paramvalue)>1
        parameterOK = 0;
    end
        
elseif paramnumber==7
     %m22
    if abs(paramvalue)>1
        parameterOK = 0;
    end
     
elseif paramnumber==8
    % m31
    if (paramvalue<0)||(paramvalue>1)
        parameterOK = 0;
    end

elseif paramnumber==9
    % m32
    if (paramvalue<0)||(paramvalue>1)
        parameterOK = 0;
    end

elseif paramnumber==10
    % q1

    
elseif paramnumber==11
    % q2

    
elseif paramnumber==12
    % q3
    if paramvalue>=0
        parameterOK = 0;
    end

elseif paramnumber==13
    %gamma
    if (paramvalue<=0)||(paramvalue>=1)
        parameterOK=0;
    end

elseif paramnumber==14
    %n0
    
elseif paramnumber==15
    %n1
    if paramvalue<=0
        parameterOK = 0;
        
    end
elseif paramnumber==16
    %n1x
    if paramvalue>0
        parameterOK = 0;
        
    end
elseif paramnumber==17
    %n1y
    if abs(paramvalue)>1
        parameterOK = 0;
        
    end
elseif paramnumber==18
    %s1
    if paramvalue<=0
        parameterOK = 0;
        
    end

elseif paramnumber==19
    %s2
elseif paramnumber==20
    %s3
elseif paramnumber==21
    %l1
    if paramvalue<=0
        parameterOK = 0;
        
    end

elseif paramnumber==22
    %l2
elseif paramnumber==23
    %l3
elseif paramnumber==24
    %alpha
    if paramvalue<0
        parameterOK = 0;     
    end
elseif paramnumber==25
    %lambdaE
    if paramvalue<0
        parameterOK = 0;     
    end
    
elseif paramnumber==26
    %lambdaN
    if paramvalue<0
        parameterOK = 0;     
    end
    
elseif paramnumber==27
    %xiE
    if paramvalue<0
        parameterOK = 0;     
    end
    
elseif paramnumber==28
    %xiN
    if paramvalue<0
        parameterOK = 0;     
    end
elseif paramnumber==29
    %L
    if (paramvalue<0)||(paramvalue>1)
        parameterOK = 0;     
    end
    
elseif paramnumber==30
    %H3
    if paramvalue<=0
        parameterOK = 0;
        
    end
elseif paramnumber==31
    %M3
    
elseif paramnumber==32
    %C
    
elseif paramnumber==33
    %tAC  
     if paramvalue<=0
        parameterOK = 0;
        
    end
elseif paramnumber==34
    %t1  
     if paramvalue<=0
        parameterOK = 0;
        
    end
elseif paramnumber==35
    %EGFover
    if paramvalue<0
        parameterOK = 0;     
    end

elseif paramnumber==36
    %D
    if paramvalue>2
        parameterOK = 0;     
    end
end

%New constraint EGF over