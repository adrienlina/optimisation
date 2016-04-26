function [F,G,ind] = OracleDG(lambda,ind)
    K = -(Ar'*pr+Ad'*lambda)./r;
    q = sign(K).*sqrt(abs(K));
    rqq = -(Ar'*pr+Ad'*lambda);
    if ind==2 then 
        F = -(1/3*q'*rqq+pr'*Ar*q+lambda'*(Ad*q-fd));
        G = %nan;
    elseif ind==3 then
        G =-(Ad*q-fd);
        F = %nan;
    elseif ind==4 then
        F = -(1/3*q'*rqq+pr'*Ar*q+lambda'*(Ad*q-fd));
        G = -(Ad*q-fd);
    else
        F = %nan;
        G = %nan;
    end
endfunction
