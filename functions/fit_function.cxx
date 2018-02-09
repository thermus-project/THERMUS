/************************************************************************
 * Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 *                                                                      *
 * Author: The THERMUS Project (A Thermal Model Package for ROOT).      *
 * Contributors (UCT-IPHC) are mentioned in the code where appropriate. *
 *                                                                      *
 * Permission to use, copy, modify and distribute this software and its *
 * documentation strictly for non-commercial purposes is hereby granted *
 * without fee, provided that the above copyright notice appears in all *
 * copies and that both the copyright notice and this permission notice *
 * appear in the supporting documentation. The authors make no claims   *
 * about the suitability of this software for any purpose. It is        *
 * provided "as is" without express or implied warranty.                *
 ************************************************************************/

// Author: Spencer Wheaton 14 July 2004         //

#include <TROOT.h>
#include <TObject.h>
#include <TMath.h>
#include <TMinuit.h>
#include <TList.h>

#ifndef Thermal_TTMThermalFit
#include <TThermalFit.h>
#endif

// **** Globals **** //

TTMThermalFit *gFit;
Int_t j;                        // # of fit parameters 
Double_t min_f = 10000000.;     // minimum chi-square or quad-deviation
Bool_t Chi;                     // true if chi-square fit
Bool_t QDev;                    // true if quad deviation fit

// ***************** //

//_______________________________________________________________________//
void Minuit_fcn(Int_t & npar, Double_t * gin, Double_t & f,
                Double_t * par, Int_t iflag);

void fit_function(TTMThermalFit * fit, Int_t flag = 0)
{
  gFit = fit;	 

  Int_t imax = 0;
  TString FitDescr = gFit->GetDescriptor();

  if(FitDescr == "GCanonical"){
    imax = 9;
  }else if(FitDescr == "SCanonical"){
    imax = 5;
  }else if(FitDescr == "BSQCanonical"){
    imax = 5;
  }


  if (flag == 0) {
    Chi = true;
    QDev = false;
  } else if (flag == 1) {
    Chi = false;
    QDev = true;
  }
   
  // Check # of fit parameters //

  j = 0;

  for (Int_t i = 0; i <= imax; i++) {
    TTMParameter *current;
    current = (gFit->GetParameterSet())->GetParameter(i);
    if (current->GetFlag() == 0) {
      j++;
    }
  }
  TMinuit *lMinuit = new TMinuit(j);
  lMinuit->SetFCN(Minuit_fcn);
  lMinuit->SetErrorDef(1.);  		// 1-sigma errors


  // Set fit parameters //

  Int_t k = 0;
  Int_t error_flag = 0;

  for (Int_t i = 0; i <= imax; i++) {
    TTMParameter *current;
    current = (gFit->GetParameterSet())->GetParameter(i);
    if (current->GetFlag() == 0) {
      lMinuit->mnparm(k, current->GetName(),
                      current->GetStart(), current->GetStep(),
                      current->GetMin(), current->GetMax(), 
                      error_flag); // error_flag = 0 if no problems
      k++;
    }
  }

  // List Experimental Yields //

  gFit->ListYields();

  // ************************************************ //
  // ************************************************ //

  // arglist[0] : max # of iterations
  // arglist[1] : convergence when EDM < arglist[1] * 10^-3
  // check error_flag afterwards to see if command executed normally

  // Other options:
  //
  // Double_t temp = 1000;
  // lMinuit->mnexcm("MIGRAD", &temp, 1, error_flag);
  // (sets max iterations to 1000)
  //            OR
  // lMinuit->mnexcm("MIGRAD", &temp, 0, error_flag);
  // (uses defaults)

  Double_t arglist[10];
  arglist[0] = 15000;
  arglist[1] = 1;
   
  lMinuit->mnexcm("MIGRAD", arglist, 2, error_flag);

   
  // ************************************************ //

  // Return the current status of the minimization
  // Check icstat for state of covariance matrix

  Double_t amin, edm, errdef;
  Int_t nvpar, nparx, icstat;
  lMinuit->mnstat(amin, edm, errdef, nvpar, nparx, icstat);

  // ************************************************ //

  // Outputs the results depending on the first flag

  lMinuit->mnprin(1, amin);

  // ************************************************ //
  // ************************************************ //

  // Now after fit update parameter set with best-fit //
  // values, re-calculate densities and display       //

  Int_t kk = 0;

  for (Int_t i = 0; i <= imax; i++) {
    TTMParameter *current;
    current = (gFit->GetParameterSet())->GetParameter(i);
    if (current->GetFlag() == 0) {
      Double_t val, error;
      lMinuit->GetParameter(kk, val, error);
      current->SetValue(val);
      current->SetError(error);
      current->SetStatus("(FITTED!)");
      kk++;
    }
  }
  
  gFit->GenerateYields();
  gFit->SetMinuit(lMinuit);
}
