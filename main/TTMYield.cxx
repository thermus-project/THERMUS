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

// Author: Spencer Wheaton 14 July 2004 //
//__________________________________________________________________________
// Object for storage of experimental yields and model predictions.
//

#include <TTMYield.h>

ClassImp(TTMYield)

//__________________________________________________________________________
TTMYield::TTMYield()
{
   fName = "";
   fID1 = 0;
   fID2 = 0;
   fExpValue = 0.;
   fExpError = 0.;
   fModelValue = 0.;
   fModelError = 0.;
   fSet1 = (TTMParticleSet *) 0;
   fSet2 = (TTMParticleSet *) 0;
   fFit = true;
}

//__________________________________________________________________________
TTMYield::TTMYield(TString name, Double_t val, Double_t err, Int_t id1,
                   Int_t id2, Bool_t fit)
{
   fName = name;
   fExpValue = val;
   fExpError = err;
   fID1 = id1;
   fID2 = id2;
   fModelValue = 0.;
   fModelError = 0.;
   fSet1 = (TTMParticleSet *) 0;
   fSet2 = (TTMParticleSet *) 0;
   fFit = fit;
}

//__________________________________________________________________________
void TTMYield::List()
{
   cout.width(20);
   cout << fName << ":" << endl;
   if(fFit){cout << "\t\t\t" << "FIT YIELD   " << endl;}
   else{cout << "\t\t\t" << "PREDICTED YIELD   " << endl;}
   cout << "\t\t\t" << "Experiment: ";
   cout.width(8);
   cout << fExpValue << "\t";
   cout.width(3);
   cout << " +- ";
   cout.width(6);
   cout << fExpError << endl;

   if (fModelValue) {
      cout << "\t\t\t" << "Model:      ";
      cout.width(8);
      cout << fModelValue << "\t";
      cout.width(3);
      cout << " +- ";
      cout.width(6);
      cout << fModelError << endl;
      cout << "\t\t\t" << "Std.Dev.: " << GetStdDev() << "  Quad.Dev.: "
          << GetQuadDev() << endl << endl;
   }
}

//__________________________________________________________________________
TTMYield TTMYield::operator=(TTMYield obj)
{
   fName = obj.GetName();
   fID1 = obj.GetID1();
   fID2 = obj.GetID2();
   fFit = obj.GetFit();
   fSet1 = obj.GetPartSet1();
   fSet2 = obj.GetPartSet2();
   fExpValue = obj.GetExpValue();
   fExpError = obj.GetExpError();
   fModelValue = obj.GetModelValue();
   fModelError = obj.GetModelError();
   return *this;
}
