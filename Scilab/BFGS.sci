function [fopt,xopt,gopt]=BFGS(Oracle,xini)


///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//         Methode de BFGS a pas variable                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////


// ------------------------
// Parametres de la methode
// ------------------------

   titre = "Parametres de la méthode BFGS";
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
   ind = 4;
   [F,G] = Oracle(x,ind);
   Id = eye(length(G),length(G)); 
   for k = 1:iter

//    - valeur du critere et du gradient
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
         Wk = Id;
      else
         deltag = G - Gk;
         deltau = x-xk;
         Wk = (Id-(deltau*deltag'/(deltag'*deltau)))*Wk*(Id-(deltag*deltau'/(deltag'*deltau)))+(deltau*deltau')/(deltag'*deltau);
         D = -Wk*G;
      end
     xk = x;
     
     Dk = D;
     Gk = G;
//    - mise a jour des variables

      
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
   disp('Fin de la methode de gradient de BFGS')
   disp(cvge)

// - visualisation de la convergence

   Visualg(logG,logP,Cout);

endfunction
