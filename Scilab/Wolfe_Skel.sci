function [alphan,ok]=Wolfe(alpha,x,D,Oracle)


//////////////////////////////////////////////////////////////
//                                                          //
//   RECHERCHE LINEAIRE SUIVANT LES CONDITIONS DE WOLFE     //
//                                                          //
//                                                          //
//  Arguments en entree                                     //
//  -------------------                                     //
//    alpha  : valeur initiale du pas                       //
//    x      : valeur initiale des variables                //
//    D      : direction de descente                        //
//    Oracle : nom de la fonction Oracle                    //
//                                                          //
//  Arguments en sortie                                     //
//  -------------------                                     //
//    alphan : valeur du pas apres recherche lineaire       //
//    ok     : indicateur de reussite de la recherche       //
//             = 1 : conditions de Wolfe verifiees          //
//             = 2 : indistinguabilite des iteres           //
//                                                          //
//                                                          //
//    omega1 : coefficient pour la 1-ere condition de Wolfe //
//    omega2 : coefficient pour la 2-eme condition de Wolfe //
//                                                          //
//////////////////////////////////////////////////////////////


// -------------------------------------
// Coefficients de la recherche lineaire
// -------------------------------------

   omega1 = 0.001;
   omega2 = 0.99;

   alphamin = 0.0;
   alphamax = %inf;

   ok = 0;
   dltx = 0.00000001;

// ---------------------------------
// Algorithme de Fletcher-Lemarechal
// ---------------------------------

   // Appel de l'oracle au point initial
   
   ind = 4;
   [F,G] = Oracle(x,ind);

   // Initialisation de l'algorithme

   alphan = alpha;
   xn     = x;

   // Boucle de calcul du pas
   //
   // xn represente le point pour la valeur courante du pas,
   // xp represente le point pour la valeur precedente du pas.
   [J0,GJ0] = Oracle(x,4);
   while ok == 0
      
      xp = xn;
      xn = x + (alphan*D);

      // Calcul des conditions de Wolfe
      [J,GJ] = Oracle(xn,4);
      condition1 = (J<=J0+omega1*alphan*(GJ0'*D));
      condition2 = ((GJ'*D)>=omega2*(GJ0'*D));
      disp(J);
      disp(J0);
      disp(GJ'*D);
      disp(GJ0'*D);
      disp("fin")
      // Test de la valeur de alphan :
      // - si les deux conditions de Wolfe sont verifiees,
      //   faire ok = 1 : on sort alors de la boucle while
      // - sinon, modifier la valeur de alphan : on reboucle.
      if ~condition2 then disp(condition2); end
      if ~(condition1) then
          alphamax = alphan;
          alphan = 1/2*(alphamin+alphamax);
      elseif ~(condition2) then
          alphamin = alphan;
          if alphamax == %inf then alphan = 2 * alphamin;
          else alphan = 1/2*(alphamax+alphamin);
          end
          disp("condition 2")
      else
          ok = 1;
      end
      

      // Test d'indistinguabilite

      if norm(xn-xp) < dltx then
        ok = 2;
      end

   end

endfunction
