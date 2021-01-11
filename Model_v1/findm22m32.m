function [m22,m32,flag] = findm22m32(m11,m12,m21,m31)

[valmax,indmax] = max([m21,m31]);

if indmax == 1
    
    m32aux(1) = (-m11*m12*m31 - sqrt((-m11^2)*(m12^2)*(m21^2) + (m21^2)*(m31^2) - (m12^2)*(m21^2)*(m31^2) + (m31^4) - (m12^2)*(m21^4)))/(m21^2 + m31^2);
    m22aux(1) = (-m11*m12-m31*m32aux(1))/m21;
    
    m32aux(2) = (-m11*m12*m31 + sqrt((-m11^2)*(m12^2)*(m21^2) + (m21^2)*(m31^2) - (m12^2)*(m21^2)*(m31^2) + (m31^4) - (m12^2)*(m21^4)))/(m21^2 + m31^2);
    m22aux(2) = (-m11*m12-m31*m32aux(2))/m21;
    
    if isreal([m22aux,m32aux])
        
        if (abs(m22aux(1))<=1)&&(m32aux(1)>0)&&(m32aux(1)<=1) % |m22| must be <=1 and m32 must be positive and also less than 1
            
            m22=m22aux(1);
            m32=m32aux(1);
            flag=1;
            
        elseif (abs(m22aux(2))<=1)&&(m32aux(2)>0)&&(m32aux(2)<=1)
            
            m22=m22aux(2);
            m32=m32aux(2);
            flag=1;
            
        else
            m22=0;
            m32 = 0;
            flag=0;
        end
    else
        
        m22=0;
        m32 = 0;
        flag=0;
        
    end
    
    
else
    
    m22aux(1) = (-m11*m12*m21 - sqrt((-m11^2)*(m12^2)*(m31^2) + (m21^2)*(m31^2) - (m12^2)*(m21^2)*(m31^2) + (m31^4) - (m12^2)*(m31^4)))/(m21^2 + m31^2);
    m32aux(1) = (-m11*m12-m21*m22aux(1))/m31;
    
    m22aux(2) = (-m11*m12*m21 + sqrt((-m11^2)*(m12^2)*(m31^2) + (m21^2)*(m31^2) - (m12^2)*(m21^2)*(m31^2) + (m31^4) - (m12^2)*(m31^4)))/(m21^2 + m31^2);
    m32aux(2) = (-m11*m12-m21*m22aux(2))/m31;
    
    if isreal([m22aux,m32aux])
        
        if (abs(m22aux(1))<=1)&&(m32aux(1)>0)&&(m32aux(1)<=1)
            
            m22=m22aux(1);
            m32=m32aux(1);
            flag=1;
            
        elseif (abs(m22aux(2))<=1)&&(m32aux(2)>0)&&(m32aux(1)<=1)
            
            m22=m22aux(2);
            m32=m32aux(2);
            flag=1;
            
        else
            m22=0;
            m32 = 0;
            flag=0;
        end
    else
        
        m22=0;
        m32 = 0;
        flag=0;
        
    end
end

