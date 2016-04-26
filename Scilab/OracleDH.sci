function [F,G,H,ind] = OracleDH(lambda,ind)
    K = -((Ar'*pr)+Ad'*lambda)./r;
    q = sign(K).*sqrt(abs(K));
    H = .5 * Ad * diag(ones(n,1) ./ (r.*abs(q)))* Ad';
    rqq = -(Ar'*pr+Ad'*lambda);
    F = %nan;
    G = %nan;
    H = %nan;
    if ind==2 then 
        F = -(1/3*q'*rqq+pr'*Ar*q+lambda'*(Ad*q-fd));
    elseif ind==3 then
        G = -(Ad*q-fd);
    elseif ind==4 then
        F = -(1/3*q'*rqq+pr'*Ar*q+lambda'*(Ad*q-fd));
        G = -(Ad*q-fd);
    elseif ind == 5 then
        H = 0.5*Ad*diag(ones(n,1)./(r.*abs(q)))*Ad';
    elseif ind == 6 then
        G = -(Ad*q-fd);
        H = 0.5*Ad*diag(ones(n,1)./(r.*abs(q)))*Ad';
    elseif ind == 7 then
        F = -(1/3*q'*rqq+pr'*Ar*q+lambda'*(Ad*q-fd));
        G = -(Ad*q-fd);
        H = 0.5*Ad*diag(ones(n,1)./(r.*abs(q)))*Ad';
        //R = (ones(sizeB(2),1)*r');
        //V = (ones(sizeB(2),1)*v');
        //H = (B'.*R.*V)*B;
    end
endfunction
