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

#ifndef Thermal_TTMThermalFit
#include <TThermalFit.h>
#endif

extern TTMThermalFit *gFit;
extern Int_t j;                 // # of fit parameters 
extern Double_t min_f;          // minimum chi-square or quad-deviation
extern Bool_t Chi;              // true if chi-square fit
extern Bool_t QDev;             // true if quad-dev fit

void Minuit_fcn(Int_t & npar, Double_t * gin, Double_t & f,
                Double_t * par, Int_t iflag)
{
  Int_t k = 0;

  // Change fit parameters //

  Int_t imax = 0;
  TString FitDescr = gFit->GetDescriptor();

  if(FitDescr == "GCanonical"){
    imax = 9;
  }else if(FitDescr == "SCanonical"){
    imax = 5;
  }else if(FitDescr == "BSQCanonical"){
    imax = 5;
  }

  for (Int_t i = 0; i <= imax; i++) {
    TTMParameter *current;
    current = (gFit->GetParameterSet())->GetParameter(i);
    if (current->GetFlag() == 0) {
      current->SetStatus("(** FITTING **)");
      current->SetValue(par[k]);
      k++;
    }
  }

  // Calculate new densities and chi-square or quad-dev //

  gFit->GenerateYields();

  cout <<"  *********************************** FITTING ************"
       <<"**********************"
       << endl << endl;

  // *** //

  for (Int_t i = 0; i <= imax; i++)
    {
      TTMParameter *current;
      current = (gFit->GetParameterSet())->GetParameter(i);
      cout << "\t\t";
      cout.width(12);
      cout << current->GetName() << " = \t";
      cout.width(10);
      cout << current->GetValue() << "\t\t";
      cout.width(15);
      cout << current->GetStatus() << endl;
    }
  cout << endl << "\t\t" << (gFit->GetParameterSet())
    ->GetConstraintInfo() << endl << endl;
  cout <<"  ****************************************************"
       <<"**************************"
       << endl << endl;

  // *** //

  cout << "\t ******************** ";

  if (Chi) {
    f = gFit->GetChiSquare();
    cout << "ChiSquare = " << f;
  } else if (QDev) {
    f = gFit->GetQuadDev();
    cout << "Quad Dev  = " << f;
  }
  cout << " ******************* " << endl << endl;
  TTMThermalModel *model = gFit->GenerateThermalModel(gFit->GetParticleSet());
  model->GenerateParticleDens();
  cout << "\t\t S/V \t = \t" << model->GetStrange() << endl;
  cout << "\t\t B/2Q \t = \t" << model->GetBaryon() / 2.
    / model->GetCharge() << endl << endl;
  cout << "\t\t C/V \t = \t" << model->GetCharm() << endl;
  cout << "\t\t b/V \t = \t" << model->GetBeauty() << endl;

  // Check for minimum and display //

  if (f < min_f) {
    min_f = f;
    cout << "\t New Minimum!" << endl << endl;
    for (Int_t i = 0; i <= imax; i++)
      {
        TTMParameter *current;
        current = (gFit->GetParameterSet())->GetParameter(i);
        cout << "\t\t";
        cout.width(12);
        cout << current->GetName() << " = \t";
        cout.width(10);
        cout << current->GetValue() << "\t\t";
        cout.width(15);
        cout << current->GetStatus() << endl;
      }
    cout << endl << "\t\t" << (gFit->GetParameterSet())
      ->GetConstraintInfo() << endl << endl;
    cout <<"  ****************************************************"
         <<"**************************"
         << endl << endl;
    gFit->ListYields();

  } else {
    cout <<"  ****************************************************"
         <<"**************************"
         << endl << endl;
  }

}
