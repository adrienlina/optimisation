function [F,G,ind] = OraclePG(qc,ind)
    //AdTInv = inv(AdT);
    //B = [-AdTInv*AdC;eye(n-md,n-md)]
    //q0 = [AdTInv*fd;zeros(n-md,1)]
    v = q0+B*qc;
    u = r.*v.*abs(v);
    if ind==2 then 
        F = 1/3*sum(u.*v)+sum(pr.*(Ar*v));
        G = %nan;
    elseif ind==3 then
        G = B'*(r.*v.*abs(v))+(Ar*B)'*pr;
        F = %nan;
    elseif ind==4 then
        F = 1/3*sum(u.*v)+sum(pr.*(Ar*v));
        G = B'*(r.*v.*abs(v))+(Ar*B)'*pr;
    else
        F = %nan;
        G = %nan;
    end
endfunction
