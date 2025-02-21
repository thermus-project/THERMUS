// Macro for THERMUS Fits
// each yields can be included or excluded.

// Includes for Compilation

#if !defined(__CLING__) || defined(__ROOTCLING__)

#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TFile.h"

#include "TTMParticleSet.h"
#include "TTMYield.h"
#include "TTMParameterSetBSQ.h"
#include "TThermalFitBSQ.h"

#include "TTMDensObj.h"

#include "TCanvas.h"
#include "TPad.h"
#include "TH2.h"
#include "TStopwatch.h"


#endif

/*
TGraph *myOneSigmaContour;
TGraph *myTwoSigmaContour;
TGraph *myThreeSigmaContour;
*/
//Int_t myDrawContours(TTMThermalFitBSQ &fitresult);
Int_t mySaveFitParameters(char * fileOutName, TTMThermalFitBSQ &fitresult);


Int_t lhc5020_fit_charm(Bool_t rWrite = 1, Bool_t gsfixed = 1, Bool_t fitMuQandS = 0, Bool_t constrainMuQ = 1, Bool_t volCor = 0,  Bool_t quantRes = 1){
    
    TStopwatch t;
    t.Start();
    
    Int_t  gDebugMode=1;
  // **************************************************
  // First, definition the particle list and their properties (and decays):
  // TTMParticleSet set("THERMUS/particles/PartList_PPB2014_CBHN.txt",true); // -> still OK !!!
  // set.InputDecays("THERMUS/particles");  // here true means the decays are scaled to sum(BR) = 100%
    TTMParticleSet set("particles/PartList_PPB2014_CBHN.txt",true); // -> still OK !!!
    set.InputDecays("particles");  // here true means the decays are scaled to sum(BR) = 100%
  
  if (volCor) set.SetRadii(0.3); // to be discussed, see personal notes
  
  // **************************************************
  // Second, choice of formalism: here, GC
  // **************************************************
  // Third, choice of starting parameters
  Double_t Tch = 0.153,  muB = 0.001, muS = 0.00, muQ = 0.0;// GeV
  Double_t gammas = 1.00; // 0.86
  Double_t radius = 10; //units = fm but no dependence for Grand Canonical.
  radius = 10.5; //units = fm but no dependence for Grand Canonical for Pb-Pb
  //radius = 3; //for p-Pb
  // R=10.9 fm give 5425 fm^3 and R=10.5 fm give 4850 fm^3
  // HBT volume ~4850 fm^3 // PLB696 (2011) 328
  Double_t volume = (4*3.1416/3.)*(radius*radius*radius);
  if (radius&&gDebugMode) printf("INFO: radius = %.1f fm, volume = %.2f fm^3\n",radius,volume);
  // Double_t B=197, Q=79, Bover2Q = B/(2*Q); if(debugMode) printf("INFO: expected colliding system: Au-Au\n");
   Double_t B=208, Q=82, Bover2Q = B/(2*Q); if(gDebugMode) printf("INFO: expected colliding system: Pb-Pb\n");
  // Double_t B=209, Q=83, Bover2Q = B/(2*Q); if(gDebugMode) printf("INFO: expected colliding system: p-Pb\n");
  Double_t muC = 0, gammac = 1;
  TTMParameterSetBSQ par(Tch,muB,muS,muQ,gammas,radius,muC,gammac);
  // **************************************************
  // Fourth, choice of parameter to fit or to fix
  // -> Default is fixed at zero or at unity so important to check !!
  //
  //--- Tch:
  //par.FixT(0.140);
  par.FitT(Tch,0.130,0.190,0.0001);
  //
  //--- muB:
  par.FixMuB(0.0001);
  //par.FitMuB(muB,0.00001,0.10);
  //
  //--- gamma s: 
  if (gsfixed) par.FixGammas(1.);
  else par.FitGammas(gammas,0.75,1.25,0.0001);
  //
  //--- gamma c&b:
  //par.FitGammac(30,10,50,0.0001);
  par.FixGammac(40.);
  par.FixGammab(40.);
  //
  //--- radius:
  //radius = 11.22;
  par.FitRadius(radius,1,15,0.01);
  //radius = pow(7000 *3. /4. /3.1416,1/3.);
  //par.FixRadius(radius); // check nucl-th/0110035 i.e. Phys.Rev.C65:027901,2002
  //
  // **************************************************
  // Fifth, option of adding some constraints e.g. B/2Q  
  if (fitMuQandS) {
    par.FitMuQ(muQ);
    //par.FitMuS(muS);
  }
  if (constrainMuQ){
    par.ConstrainMuQ(Bover2Q);
    printf("INFO: constraining B/2Q = %.4f\n",Bover2Q);
  }
  if(gDebugMode) par.List();

  // **************************************************
  // Fifth, Create an instance of the Fit
  // -> Branch the experimental values
  //  TTMThermalFitBSQ fit(&set,&par,"rhic200_exp.txt");
  char fileName[60] = "./tests/lhc5020_final_0_single_charm.txt";
  //char fileName[60] = "./lhcAP_cleaned_no_anti.txt";
  
  printf("INFO: ***** this is the input filename:  %s ***** \n",fileName);

  TTMThermalFitBSQ fit(&set,&par,fileName);
  // Specify excluded volume, quantum statistics and resonance width treatments
  fit.SetExclVol(volCor);
  fit.SetQStats(quantRes);
  fit.SetWidth(quantRes);
  // -> Switch condition for ratios and yields exclusions

    //save
    char fileOutName[100] = "./current_fit_lhc5020_fit_test.txt";
    
  // specific conditions for Pb-Pb
   fit.GetYield(313,0,"ALICE")->Predict();                   //K*
   fit.GetYield(333,0,"ALICE")->Predict();                   //phi

  // fit.GetYield(2212,0,"ALICE")->Predict();                  //p
  // fit.GetYield(-2212,0,"ALICE")->Predict();
  // fit.GetYield(211,0,"ALICE")->Predict();                   //pi+
  // fit.GetYield(-211,0,"ALICE")->Predict();

  // fit.GetYield(321,0,"ALICE")->Predict();                   //K+
  // fit.GetYield(-321,0,"ALICE")->Predict();
  // fit.GetYield(310,0,"ALICE")->Predict();                   //K0S
  // fit.GetYield(3122,0,"ALICE")->Predict();                  //Lambda
  // fit.GetYield(-3122,0,"ALICE")->Predict();

  // fit.GetYield(3312,0,"ALICE")->Predict();                //Xi
  // fit.GetYield(-3312,0,"ALICE")->Predict();
  // fit.GetYield(3334,0,"ALICE")->Predict();                //Omega
  // fit.GetYield(-3334,0,"ALICE")->Predict();

  // fit.GetYield(1000010020,0,"ALICE")->Predict();            //d
  // fit.GetYield(-1000010020,0,"ALICE")->Predict();
  // fit.GetYield(1000020030,0,"ALICE")->Predict();            //HE3
    
   fit.GetYield(411,0,"ALICE")->Predict();            //D+
   fit.GetYield(421,0,"ALICE")->Predict();            //D0
   fit.GetYield(413,0,"ALICE")->Predict();            //D*+

  
  if(gDebugMode) printf("INFO: now generate Yields\n");
  fit.GenerateYields();
  
  if(gDebugMode) fit.ListYields();
  
  if (gDebugMode==2) return gDebugMode;
  
  if(gDebugMode) printf("INFO: this is the chi2 %.3f and quadratic deviation %.3f \n",fit.GetChiSquare(),fit.GetQuadDev());
  printf("INFO: finally we performe a chi2 fit\n"); // means FitData(0) (for Quadratic dev fit -> 1)
  
  fit.FitData(0);
  
  fit.GetParameterSet()->List();
  if(gDebugMode){
    printf("INFO: print the output of MINUIT !! \n");
    printf("  ****************************************************************************** \n");
    fit.ListMinuitInfo();
  }
    
    //Printing Yields******************************************************
    printf("  ****************************************************************************** \n");
    printf("INFO: Listing yields from the model \n");
    printf(" %-15s  %-20s  mod. value \t Std.Dev.\n","id","name");
    TTMYield *lYield = 0;
    Int_t idYields[]={211,321,310,2212,3122,3312,3334,313,333,1000010020,1010010030,-1010010030,1000020030,411,421,413};
    Int_t nCharges = sizeof(idYields)/sizeof(Int_t);
    
    for(Int_t iYields = 0;iYields<nCharges;iYields++)
    {
        lYield = fit.GetYield(idYields[iYields],0,"ALICE");
        printf(" %-15d  %-20s  %6f \t %.2f \n",lYield->GetID1(),lYield->GetTMName().Data(),lYield->GetModelValue(),lYield->GetStdDev());
    }
    //Printing Yields******************************************************
    
  printf("  ****************************************************************************** \n");

  if(gsfixed) printf("INFO: gamma_s was fixed to unity \n");
  else printf("INFO: gamma_s was a free parameter \n");
  if(volCor) printf("INFO: volume corrections were included \n");
  else printf("INFO: volume corrections were not included \n");
  if(quantRes) printf("INFO: quantum Stat and Resonance Width were included \n");
  else printf("INFO: quantum Stat and Resonance Width were not included \n");
    
    radius = fit.GetParameterSet()->GetRadiusPar()->GetValue();
    Float_t radiusError = fit.GetParameterSet()->GetRadiusPar()->GetError();
    volume = (4*3.1416/3.)*radius*radius*radius;
    Float_t volume_hi = (4*3.1416/3.)*(radius+radiusError)*(radius+radiusError)*(radius+radiusError);
    Float_t volume_lo = (4*3.1416/3.)*(radius-radiusError)*(radius-radiusError)*(radius-radiusError);
    if (radius) printf("INFO: radius = %.2f +- %.2f fm, volume = %.0f +%.0f -%.0f fm^3 \n",radius,radiusError,volume,volume_hi-volume,volume-volume_lo);

    
    if(rWrite){
    Bool_t saveOk = mySaveFitParameters(fileOutName, fit);
    
  if(saveOk) printf("***** done with saving fit information \n\n");
    }
  //Bool_t plotOk = myDrawContours(fit);
  //if(plotOk) printf("***** done with plotting contours \n\n");
    
    t.Stop();
    t.Print();
    
  return gDebugMode;
}

/*
Int_t myDrawContours(TTMThermalFitBSQ &fitresult){
  TCanvas *myCanvasContour = new TCanvas("myCanvasContour","ROOT Contour Fit/Calculation",200,5,700,700);
  myCanvasContour->Draw();

  // Set the Pads
  TPad *myPadContour = new TPad("myPadContour","myPadContour",0.0,0.0,1.0,1.0,0);
  myPadContour->SetLeftMargin(0.17);
  myPadContour->SetTopMargin(0.02);
  myPadContour->SetRightMargin(0.03);
  myPadContour->SetBottomMargin(0.13);
  myPadContour->Draw();
  myPadContour->SetTickx(1);
  myPadContour->SetTicky(1);
  myPadContour->cd();

  TMinuit *access = (TMinuit*)fitresult.GetMinuit();
  Int_t accessOk = 0;
  
  if(access){
    accessOk = 1;
    TH2F *hMyBlankContour = new TH2F("hMyBlankContour","; radius (fm); T (GeV)",10,0,20,10,0.140,0.180);
    hMyBlankContour->SetNdivisions(505,"x");
    hMyBlankContour->SetNdivisions(505,"y");
    hMyBlankContour->SetTitleOffset(1.3,"y");
    hMyBlankContour->SetStats(0);
    hMyBlankContour->Draw();
    
    access->SetErrorDef(1);
    myOneSigmaContour = (TGraph*)access->Contour(50,1,0);
    if (myOneSigmaContour) {
      printf("INFO: ***** let's check the contour 1 sigma... \n");
      myOneSigmaContour->Draw("SAME L");
      myOneSigmaContour->SetLineColor(2);
      myOneSigmaContour->SetLineStyle(1);
      myOneSigmaContour->SetLineWidth(1);
    }
    access->SetErrorDef(4);
    myTwoSigmaContour = (TGraph*)access->Contour(50,1,0);
    if (myTwoSigmaContour) {
      printf("INFO: ***** let's check the contour 2 sigma... \n");
      myTwoSigmaContour->Draw("SAME L");
      myTwoSigmaContour->SetLineColor(6);
      myTwoSigmaContour->SetLineStyle(1);
      myTwoSigmaContour->SetLineWidth(1);
    }
    access->SetErrorDef(9);
    myThreeSigmaContour = (TGraph*)access->Contour(50,1,0);
    if (myThreeSigmaContour) {
      printf("INFO: ***** let's check the contour 3 sigma... \n");
      myThreeSigmaContour->Draw("SAME L");
      myThreeSigmaContour->SetLineColor(9);
      myThreeSigmaContour->SetLineStyle(1);
      myThreeSigmaContour->SetLineWidth(1);
    }
  }
  return accessOk;
}
*/
Int_t mySaveFitParameters(char * fileOutName,
                          TTMThermalFitBSQ &fitresult){
    Int_t Ndf = 0;
    Float_t volCorvalue = 0.3;
    
  printf("\n \nINFO: ***** this is the output filename:  %s ***** \n",fileOutName);

  TTMYield *lYield = 0;
    Int_t idYields[16] = {211,321,310,2212,3122,3312,3334,313,333,1000010020,1010010030,-1010010030,1000020030,411,421,413};
  Int_t nCharges = 16;

  printf("  ****************************************************************************** \n");
    printf("INFO: Generating output file: info and tables for plotting with single charges \n");
    printf("  ****************************************************************************** \n");
    ofstream fileOut(fileOutName);
 //   fileOut<<"Grand Canonical:\n";
    fileOut<<"Ndf\t= "<<Form("%d", Ndf)<<endl;
    fileOut<<"chi2\t= "<<Form("%.1f",fitresult.GetChiSquare())<<endl;
    fileOut<<"VolCor\t= "<<Form("%.1f",volCorvalue)<<endl;
    fileOut<<"T\t= "<<fitresult.GetParameterSet()->GetTPar()->GetValue()<<" +- "<<fitresult.GetParameterSet()->GetTPar()->GetError()<<endl;
    fileOut<<"gammas\t= "<<fitresult.GetParameterSet()->GetGammasPar()->GetValue()<<" +- "<<fitresult.GetParameterSet()->GetGammasPar()->GetError()<<endl;
    Float_t radius = fitresult.GetParameterSet()->GetRadiusPar()->GetValue();
    Float_t radiusError = fitresult.GetParameterSet()->GetRadiusPar()->GetError();
    Float_t volume = (4*3.1416/3.)*radius*radius*radius;
    Float_t volume_hi = (4*3.1416/3.)*(radius+radiusError)*(radius+radiusError)*(radius+radiusError);
    Float_t volume_lo = (4*3.1416/3.)*(radius-radiusError)*(radius-radiusError)*(radius-radiusError);
 //   if (radius){
      fileOut<<"radius\t= "<<Form("%.2f",radius)<<" +- "<<Form("%.2f",radiusError)<<endl;
      fileOut<<"volume\t = "<<Form("%.0f",volume)<<" + "<<Form("%.0f",volume_hi-volume)<<" - "<<Form("%.0f",volume-volume_lo)<<endl;
//      if (volCor) fileOut<<" with EVC\n"; else fileOut<<" no EVC\n";
//    }
    
    fileOut<<"Model Value"<<endl;
    
    for(Int_t iYields = 0;iYields<nCharges;iYields++)
    {
      lYield = fitresult.GetYield(idYields[iYields],0,"ALICE");
        fileOut << Form(" %-15d",lYield->GetID1())
        <<Form("%-20s",lYield->GetTMName().Data())
        <<"   "<<Form("%7f",lYield->GetModelValue())
        <<" \t"<<Form("%7f",lYield->GetStdDev())<<endl;
        /*
         if (TMath::Abs(lYield->GetID1())<5000)
         fileOut << Form(" %d",lYield->GetID1())
         <<" \t \t"<<Form("%7f",lYield->GetModelValue())
         <<" \t"<<Form("%7f",lYield->GetStdDev())<<endl;
         else
         fileOut << Form(" %d",lYield->GetID1())
         <<" \t"<<Form("%7f",lYield->GetModelValue())
         <<" \t"<<Form("%7f",lYield->GetStdDev())<<endl;
         */
    }
    
    fileOut<<"Copy and Paste Following: "<<endl;
    
    fileOut<<"Exp.Value.\t= {";
    for(Int_t iYields = 0;iYields<nCharges;iYields++)
    {
        lYield = fitresult.GetYield(idYields[iYields],0,"ALICE");
        fileOut<<Form("%g,",lYield->GetExpValue());
    }
    fileOut<<"\b};\n";
    fileOut<<"Exp.Error.\t= {";
    for(Int_t iYields = 0;iYields<nCharges;iYields++)
    {
        lYield = fitresult.GetYield(idYields[iYields],0,"ALICE");
        fileOut<<Form("%g,",lYield->GetExpError());
    }
    fileOut<<"\b};\n";
    fileOut<<"Mod.Value.\t= {";
    for(Int_t iYields = 0;iYields<nCharges;iYields++)
    {
        lYield = fitresult.GetYield(idYields[iYields],0,"ALICE");
        fileOut<<Form("%.6f,",lYield->GetModelValue());
    }
    fileOut<<"\b};\n";
    fileOut<<"Std.Dev.\t= {";
    for(Int_t iYields = 0;iYields<nCharges;iYields++)
    {
        lYield = fitresult.GetYield(idYields[iYields],0,"ALICE");
        fileOut<<Form("%.2f,",lYield->GetStdDev());
    }
    fileOut<<"\b};\n";
    
    
    fileOut.close();
  
  return true;
}
