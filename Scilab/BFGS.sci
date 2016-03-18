function [fopt,xopt,gopt]=BFGS(Oracle,xini)


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//         Methode de gradient a pas variable                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// ------------------------
// Parametres de la methode
// ------------------------

   titre = "Parametres de la méthode Polake-Ribière";
   labels = ["Nombre maximal d''iterations";...
             "Valeur du pas de gradient";...
             "Seuil de convergence sur ||G||"];
   typ = list("vec",1,"vec",1,"vec",1);
   default = ["5000";"1";"0.000001"];
   [ok,iter,alphai,tol] = getvalue(titre,labels,typ,default);

// ----------------------------
// Initialisation des variables
// ----------------------------
   
   logG = [];
   logP = [];
   Cout = [];

   timer();

// -------------------------
// Boucle sur les iterations
// -------------------------

   x = xini;
   alpha = alphai;
   Gk= %nan;
   Dk = %nan;
   xk = %nan;
   kstar = iter;
   for k = 1:iter

//    - valeur du critere et du gradient

      ind = 4;
      [F,G] = Oracle(x,ind);

//    - test de convergence
      //disp(G);
      //disp(F);
      //disp(x);
      if norm(G) <= tol then
         kstar = k;
         break
      end

//    - calcul de la direction de descente
      if k == 1 then
         D = -G;
         Wk = eye(length(G),length(G));
      else
         deltag = G - Gk;
         deltau = x-xk;
          if k >= 169 then
          disp(norm(eye(length(G),length(G))-(deltau*deltau')/(deltag'*deltau)));
          disp(norm(deltau*deltag'/(deltag'*deltau) ));
          disp(norm(Wk));
          disp("fin");
      end
         Wk = (eye(length(G),length(G))-(deltau*deltag'/(deltag'*deltau)))*Wk*(eye(length(G),length(G))-(deltag*deltau'/(deltag'*deltau)))+(deltau*deltau')/(deltag'*deltau);
         D = -Wk*G;
      end
     xk = x;
     
     Dk = D;
     Gk = G;
//    - mise a jour des variables
      if k >= 169 then
          //disp(deltag);
          //disp(deltau);
          //disp(Wk);
      end
      
      [alpha,ok] = Wolfe(alphai,x,D,Oracle);
      x = x + (alpha*D);

//    - evolution du gradient, du pas et du critere

      logG = [ logG ; log10(norm(G)) ];
      logP = [ logP ; log10(alpha) ];
      Cout = [ Cout ; F ];

   end

// ---------------------------
// Resultats de l'optimisation
// ---------------------------

   fopt = F;
   xopt = x;
   gopt = G;

   tcpu = timer();

   cvge = ['Iteration         : ' string(kstar);...
           'Temps CPU         : ' string(tcpu);...
           'Critere optimal   : ' string(fopt);...
           'Norme du gradient : ' string(norm(gopt))];
   disp('Fin de la methode de gradient a pas fixe')
   disp(cvge)

// - visualisation de la convergence

   Visualg(logG,logP,Cout);

endfunction
