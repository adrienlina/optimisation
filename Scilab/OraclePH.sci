function [F,G,H,ind] = OraclePH(qc,ind)
    //AdTInv = inv(AdT);
    //B = [-AdTInv*AdC;eye(n-md,n-md)]
    //q0 = [AdTInv*fd;zeros(n-md,1)]
    v = q0+B*qc;
    u = r.*v.*abs(v);
    F = %nan;
    G = %nan;
    H = %nan;
    if ind==2 then 
        F = 1/3*sum(u.*v)+sum(pr.*(Ar*v));
    elseif ind==3 then
        G = B'*(r.*v.*abs(v))+(Ar*B)'*pr;
    elseif ind==4 then
        F = 1/3*sum(u.*v)+sum(pr.*(Ar*v));
        G = B'*(r.*v.*abs(v))+(Ar*B)'*pr;
    elseif ind == 5 then
        sizeB = size(B);
        RV = diag(abs(v).*r);
        H = 2*B'*(RV*B);
        //H = 2*(B'*repmat(r',(size(r))))*(B.*abs(repmat(v',(size(v)))));
    elseif ind == 6 then
        G = B'*(r.*v.*abs(v))+(Ar*B)'*pr;
        RV = diag(abs(v).*r);
        H = 2*B'*(RV*B);
    elseif ind == 7 then
        F = 1/3*sum(u.*v)+sum(pr.*(Ar*v));
        G = B'*(r.*v.*abs(v))+(Ar*B)'*pr;
        RV = diag(abs(v).*r);
        H = 2*B'*(RV*B);
        //disp(size(F),size(G),size(H))
        //R = (ones(sizeB(2),1)*r');
        //V = (ones(sizeB(2),1)*v');
        //H = (B'.*R.*V)*B;
    end
endfunction
