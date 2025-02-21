// a simple test macro for THERMUS Predictions

// root[0] .L test/prediction.C++
// root[1] prediction();

// or simply:
// root[0] .x test/prediction.C

// Includes for Compilation
#if !defined(__CLING__) || defined(__ROOTCLING__)
#include "TTMParticleSet.h"
#include "TTMYield.h"
#include "TTMParameterSetBSQ.h"
#include "TThermalFitBSQ.h"
#endif

// **************************************************
// Always better to start with a verbose mode
Bool_t debugMode=1; // verbose with debugMode=1

void all_predictions(){

  Bool_t constrainMuQ = 1, quantRes = 1;
  // **************************************************
  // First, definition the particle list
  // and their properties (and decays):
//  TTMParticleSet set(THERMUS+"/share/Thermus/particles/PartList_PPB2014_CBHN.txt",true); // -> still OK !!!
//  set.InputDecays(THERMUS+"/share/Thermus/particles");  // here true means the decays are scaled to sum(BR) = 100%
  
   TTMParticleSet set("local/share/Thermus/particles/PartList_PPB2014_CBHN.txt",true); // -> still OK !!!

  set.InputDecays("/local/share/Thermus/particles");  // here true means the decays are scaled to sum(BR) = 100%
  // **************************************************
  // Second, choice of formalism:
  // - We want here a Grand Canonical Treatment so we use
  //  a TTMThermalFitBSQ as a FIT instance with starting parameters
  //  other FIT classes are:
  // -> TTMThermalFitBQ for Strangeness Canonical;
  // -> TTMThermalFitCanBSQ for Full Canonical \n"); 
  //
  // Third, choice of starting parameters (nucl-ex/0403014)
  // Double_t Tch = 0.157,  muB = 0.0022, muS = 0.0038, muQ = 0.0;// GeV
  Double_t Tch = 0.156,  muB = 0.0001, muS = 0.00, muQ = 0.0;// GeV
  Double_t gammas = 1.00; // 0.86
  Double_t radius = 10.5; //units = fm but no dependence for Grand Canonical.
  // R=10.9 fm give 5425 fm^3 and R=10.5 fm give 4850 fm^3
  Double_t volume = (4*3.1416/3.)*radius*radius*radius;
  //  double B=197, Q=79, Bover2Q = B/(2*Q); if(debugMode) printf("INFO: expected colliding system: Au-Au\n");
  double B=208, Q=82, Bover2Q = B/(2*Q); if(debugMode) printf("INFO: expected colliding system: Pb-Pb\n");
  Double_t muC = 0, gammac = 1;
  TTMParameterSetBSQ par(Tch,muB,muS,muQ,gammas,radius,muC,gammac);
  if(debugMode) par.List();
  // **************************************************
  // Fourth, choice of parameter to fit or to fix
  // -> Default is fixed at zero or at unity so important to check !! 
  //
  //--- Tch:
  par.FixT(Tch);
  //
  //--- muB:
  par.FixMuB(muB);
  //
  //--- gamma s:
  par.FixGammas(gammas);
  //
  //--- gamma c&b:
  par.FixGammac(1.);
  //par.FixGammab(0.);
  //
  //--- radius:
  par.FixRadius(radius);
  //
  // **************************************************
  // Fifth, option of adding some constraints e.g. B/2Q
  if (constrainMuQ){
    par.ConstrainMuQ(Bover2Q);
    printf("INFO: constraining B/2Q = %.4f\n",Bover2Q);
  }
  if (radius&&debugMode) printf("INFO: radius = %.1f fm, volume = %.2f fm^3\n",radius,volume);
  // **************************************************
  // Fifth, Create an instance as if a fit is performed
  // -> Branch the fake experimental values

  TTMThermalFitBSQ fit(&set,&par,"local/share/doc/Thermus/tests/prediction.txt");

  // -> Turn off default quantum statistics and resonance width treatment
  if(!quantRes){
    fit.SetQStats(kFALSE);
    fit.SetWidth(kFALSE);
  }
  if(debugMode) printf("INFO: now generate Yields\n");
  fit.GenerateYields();
  
  if(debugMode) fit.ListYields();
  if(quantRes) printf("INFO: quantum Stat and Resonance Width are included \n");
  else printf("INFO: quantum Stat and Resonance Width are not included \n");
  radius = fit.GetParameterSet()->GetRadiusPar()->GetValue();
  volume = (4*3.1416/3.)*radius*radius*radius;
  if (radius&&debugMode) printf("INFO: radius = %.1f fm, volume = %.2f fm^3\n",radius,volume);

  printf("\n\n predicted values\n");
  printf(" id \t name \t\t mod. value \n");
  TTMYield *lYield = 0;
  Int_t idYields[30]={211,-211,321,-321,310,2212,-2212,3122,3312,-3312,3334,-3334,313,333,3124,
    3222,3224,3214,3114,3324,
    1000010020,-1000010020,1000020030,1010010030,-1010010030,411,421,4122,443,553};
  for(Int_t iYields = 0;iYields<25;iYields++)
    {
      lYield = fit.GetYield(idYields[iYields],0,"PREDICTION");
      printf(" %d \t %s \t %6f \n",lYield->GetID1(),lYield->GetTMName().Data(),lYield->GetModelValue());
    }
  
  TTMYield *lRatio = 0;
  lRatio = fit.GetYield(-211,211,"PREDICTION");
  printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  lRatio = fit.GetYield(-321,321,"PREDICTION");
  printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  lRatio = fit.GetYield(-2212,2212,"PREDICTION");
  printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  
  Int_t idRatio[7]={-321,-2212,3122,-3122,3312,-3312,3334};
  for(Int_t iRatio = 0;iRatio<7;iRatio++)
  {
    lRatio = fit.GetYield(idRatio[iRatio],-211,"PREDICTION");
    printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  }
  lRatio = fit.GetYield(-3334,3334,"PREDICTION");
  printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  lRatio = fit.GetYield(3122,310,"PREDICTION");
  printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  lRatio = fit.GetYield(2224,2212,"PREDICTION");
  printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  lRatio = fit.GetYield(333,-321,"PREDICTION");
  printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  lRatio = fit.GetYield(313,-321,"PREDICTION");
  printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  lRatio = fit.GetYield(3124,3122,"PREDICTION");
  printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  lRatio = fit.GetYield(3124,3122,"PREDICTION");
  printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  lRatio = fit.GetYield(3224,3122,"PREDICTION");
  printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  lRatio = fit.GetYield(3224,3222,"PREDICTION");
  printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  lRatio = fit.GetYield(3324,3312,"PREDICTION");
  printf(" %d \t %d \t %s \t %6f \n",lRatio->GetID1(),lRatio->GetID2(),lRatio->GetTMName().Data(),lRatio->GetModelValue());
  
}
