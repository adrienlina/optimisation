///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  MONITEUR D'ENCHAINEMENT POUR LE CALCUL DE L'EQUILIBRE D'UN RESEAU D'EAU  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// --------------------------------------
// Dimensionnement de l'espace de travail
// --------------------------------------

   stacksize(10000000);

// ------------------------------------------
// Fonctions fournies dans le cadre du projet
// ------------------------------------------

   // Donnees du problemes

   exec('Probleme_S.sce');
   exec('Structures_S.sce');
   
   // Affichage des resultats

   exec('Visualg.sci');
   
   // Verification  des resultats

   exec('HydrauliqueP.sci');
   exec('HydrauliqueD.sci');
   exec('Verification.sci');

// ------------------------------------------
// Fonctions a ecrire dans le cadre du projet
// ------------------------------------------

   // ---> Charger les fonctions  associees a l'oracle du probleme,
   //      aux algorithmes d'optimisation et de recherche lineaire.
   //
   // Exemple : la fonction "optim" de Scilab
   //
   exec('OraclePH.sci');
   exec('OraclePG.sce');
   exec('OracleDG.sci');
   exec('OracleDH.sci');
   exec('Wolfe_Skel.sci');
   exec('Optim_Scilab.sci');
   exec('Gradient_F.sci');
   exec('Gradient_V.sci');
   exec('polak-ribiere.sce');
   exec('BFGS.sci');
   exec('Newton.sci');
   
   titrgr = "Fonction optim de Scilab sur le probleme primal";

   // -----> A completer...
   // -----> A completer...
   // -----> A completer...

// ------------------------------
// Initialisation de l'algorithme
// ------------------------------

   // La dimension (n-md) est celle du probleme primal

   xini = 0.1 * rand(n-md,1);

// ----------------------------
// Minimisation proprement dite
// ----------------------------

   // Exemple : la fonction "optim" de Scilab
   //[fopt,xopt,gopt] = Optim_Scilab(OraclePG,xini);

   // Exemple : le gradient à pas fixe
   //[fopt,xopt,gopt] = Gradient_F(OraclePG,xini);
   
   // Exemple : le gradient à pas variable
   //[fopt,xopt,gopt] = Gradient_V(OraclePG,xini);
   
   // Exemple : la méthode de Polak-Ribière
   //[fopt,xopt,gopt] = polak_ribiere(OraclePG,xini);
   
   // Exemple : la méthode de BFGS
   //[fopt,xopt,gopt] = BFGS(OraclePG,xini);
   
   // Exemple : la méthode de Newton
   //[fopt,xopt,gopt] = Newton(OraclePH,xini);
 
    //[q,z,f,p] = HydrauliqueP(xopt);
 
// ----------------------------
// Minimisation problème dual
// ----------------------------

lambdaini = 0.1 * rand(md,1);

// Exemple : la fonction "optim" de Scilab
//   [fopt,xopt,gopt] = Optim_Scilab(OracleDG,lambdaini);

// Exemple : le gradient à pas fixe
//   [fopt,xopt,gopt] = Gradient_F(OracleDG,lambdaini);

// Exemple : le gradient à pas variable
   //[fopt,xopt,gopt] = Gradient_V(OracleDG,lambdaini);
   
// Exemple : la méthode de Polak-Ribière
//   [fopt,xopt,gopt] = polak_ribiere(OracleDG,lambdaini);
   
// Exemple : la méthode de BFGS
   //[fopt,xopt,gopt] = BFGS(OracleDG,lambdaini);
   
     // Exemple : la méthode de Newton
   [fopt,xopt,gopt] = Newton(OracleDH,lambdaini);


// --------------------------
// Verification des resultats
// --------------------------

   [q,z,f,p] = HydrauliqueD(xopt);

   Verification(q,z,f,p);

//
