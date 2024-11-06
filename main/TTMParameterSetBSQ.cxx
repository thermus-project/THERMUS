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

// Author: Spencer Wheaton 26 April 2010 //
//__________________________________________________________________________
// The parameter set to be used when treating B, S, Q, C and b grand-canonically.
// S/V, B/2Q, C/V and b/V may be used to constrain muS, muQ, muC and mub respectively 
// when this parameter set is linked to a particle set in a TThermalModelBSQ 
// object.
//

#include <TTMParameterSetBSQ.h>

ClassImp(TTMParameterSetBSQ)

//__________________________________________________________________________
TTMParameterSetBSQ::TTMParameterSetBSQ(Double_t temp, Double_t mub,
                                         Double_t mus, Double_t muq,
                                         Double_t gs, Double_t r, Double_t muc,
                                         Double_t gc, Double_t mubeauty, Double_t gb, Double_t b2q, Double_t s,
                                         Double_t c, Double_t beauty, Double_t temp_error,
                                         Double_t mub_error,
                                         Double_t mus_error,
                                         Double_t muq_error,
                                         Double_t gs_error, Double_t r_error,
					 Double_t muc_error, Double_t gc_error, Double_t mubeauty_error, 
					 Double_t gb_error)
{
  // Sets all parameters and their errors as well as B/2Q, S/V, C/V and b/V.
  // All parameters are set as "fixed type".
  //
   
  fB2Q = b2q;
  fSDens = s;
  fCDens = c;
  fbDens = beauty;

  fPar = fParArray;

  fParArray[0].SetParameter("T", temp, temp_error);
  fParArray[1].SetParameter("muB", mub, mub_error);
  fParArray[2].SetParameter("muS", mus, mus_error);
  fMuSConstrain = false;
  fParArray[3].SetParameter("muQ", muq, muq_error);
  fMuQConstrain = false;
  fParArray[4].SetParameter("gammas", gs, gs_error);
  fParArray[5].SetParameter("radius", r, r_error);
  fParArray[6].SetParameter("muC", muc, muc_error);
  fMuCConstrain = false;
  fParArray[7].SetParameter("gammac", gc, gc_error);
  fParArray[8].SetParameter("mub", mubeauty, mubeauty_error);
  fMubConstrain = false;
  fParArray[9].SetParameter("gammab", gb, gb_error);

  fConstraintInfo = "Parameters unconstrained";
}

//__________________________________________________________________________
TTMParameterSetBSQ::TTMParameterSetBSQ()
{
  // Sets parameter names and values and errors to 0.
  //
	 
  fB2Q = 0.;
  fSDens = 0.;
  fCDens = 0.;
  fbDens = 0.;

  fPar = fParArray;

  fParArray[0].SetParameter("T", 0., 0.);
  fParArray[1].SetParameter("muB", 0., 0.);
  fParArray[2].SetParameter("muS", 0., 0.);
  fMuSConstrain = false;
  fParArray[3].SetParameter("muQ", 0., 0.);
  fMuQConstrain = false;
  fParArray[4].SetParameter("gammas", 0., 0.);
  fParArray[5].SetParameter("radius", 0., 0.);
  fParArray[6].SetParameter("muC", 0., 0.);
  fMuCConstrain = false;
  fParArray[7].SetParameter("gammac", 0., 0.);
  fParArray[8].SetParameter("mub", 0., 0.);
  fMubConstrain = false;
  fParArray[9].SetParameter("gammab", 0., 0.);
  fConstraintInfo = "Parameters unconstrained";
}

//__________________________________________________________________________
void TTMParameterSetBSQ::ConstrainMuS(Double_t x)
{
  // Changes muS to a constrained type parameter. x is the initial 
  // strangeness density in the system.
  // 	
   
  fSDens = x;
  fMuSConstrain = true;
  fParArray[2].Constrain();
}

//__________________________________________________________________________
void TTMParameterSetBSQ::ConstrainMuQ(Double_t x)
{
  // Changes muQ to a constrained type parameter. x is the initial B/2Q 
  // ratio in the system.
  // 	
	
  fB2Q = x;
  fMuQConstrain = true;
  fParArray[3].Constrain();
}

//__________________________________________________________________________
void TTMParameterSetBSQ::ConstrainMuC(Double_t x)
{
  // Changes muC to a constrained type parameter. x is the initial charm 
  // density in the system.
  // 	
	
  fCDens = x;
  fMuCConstrain = true;
  fParArray[6].Constrain();
}

//__________________________________________________________________________
void TTMParameterSetBSQ::ConstrainMub(Double_t x)
{
  // Changes mub to a constrained type parameter. x is the initial beauty 
  // density in the system.
  // 	
	
  fbDens = x;
  fMubConstrain = true;
  fParArray[8].Constrain();
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitT(Double_t start, Double_t min,
                              Double_t max, Double_t step)
{
  // Changes T to fit type parameter
  // 	
   
  fParArray[0].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitMuB(Double_t start, Double_t min,
                                Double_t max, Double_t step)
{
  // Changes muB to fit type parameter
  // 	

  fParArray[1].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitMuS(Double_t start, Double_t min,
                                Double_t max, Double_t step)
{
  // Changes muS to fit type parameter
  // 	
	
  fMuSConstrain = false;
  fParArray[2].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitMuQ(Double_t start, Double_t min,
                                Double_t max, Double_t step)
{
  // Changes muQ to fit type parameter
  // 	
   
  fMuQConstrain = false;
  fParArray[3].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitGammas(Double_t start, Double_t min,
                                   Double_t max, Double_t step)
{
  // Changes gammas to fit type parameter
  // 	

  fParArray[4].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitRadius(Double_t start, Double_t min,
                                   Double_t max, Double_t step)
{
  // Changes radius to fit type parameter
  // 	
	
  fParArray[5].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitMuC(Double_t start, Double_t min,
                                Double_t max, Double_t step)
{
  // Changes muC to fit type parameter
  // 	

  fMuCConstrain = false;
  fParArray[6].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitGammac(Double_t start, Double_t min,
                                   Double_t max, Double_t step)
{
  // Changes gammac to fit type parameter
  // 	

  fParArray[7].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitMub(Double_t start, Double_t min,
                                Double_t max, Double_t step)
{
  // Changes mub to fit type parameter
  // 	

  fMubConstrain = false;
  fParArray[8].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FitGammab(Double_t start, Double_t min,
                                   Double_t max, Double_t step)
{
  // Changes gammab to fit type parameter
  // 	

  fParArray[9].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixT(Double_t value, Double_t error)
{
  // Fixes T
  //
	
  fParArray[0].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixMuB(Double_t value, Double_t error)
{

  // Fixes muB
  //
   
  fParArray[1].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixMuS(Double_t value, Double_t error)
{
  // Fixes muS
  //
	
  fMuSConstrain = false;
  fParArray[2].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixMuQ(Double_t value, Double_t error)
{
  // Fixes muQ
  //
	 
  fMuQConstrain = false;
  fParArray[3].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixGammas(Double_t value, Double_t error)
{
  // Fixes gammas
  //	
   
  fParArray[4].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixRadius(Double_t value, Double_t error)
{
  // Fixes radius
  //
	 	
  fParArray[5].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixMuC(Double_t value, Double_t error)
{

  // Fixes muC
  //

  fMuCConstrain = false;
  fParArray[6].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixGammac(Double_t value, Double_t error)
{
  // Fixes gammac
  //	
   
  fParArray[7].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixMub(Double_t value, Double_t error)
{

  // Fixes mub
  //

  fMubConstrain = false;
  fParArray[8].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::FixGammab(Double_t value, Double_t error)
{
  // Fixes gammab
  //	
   
  fParArray[9].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBSQ::List()
{	
  cout << "  ***************************** Thermal Parameters ";
  cout << "**************************** " << endl << endl;
  for (Int_t i = 0; i <= 9; i++) {
    fParArray[i].List();
    if (i == 2 && fMuSConstrain) {
      cout << "\t\t\t\t\t\t\t\t" << " S/V:    " << fSDens << endl << endl;
    } else if (i == 3 && fMuQConstrain) {
      cout << "\t\t\t\t\t\t\t\t" << " B/2Q: " << fB2Q << endl << endl;
    } else if (i == 6 && fMuCConstrain) {
      cout << "\t\t\t\t\t\t\t\t" << " C/V: " << fCDens << endl << endl;
    } else if (i == 8 && fMubConstrain) {
      cout << "\t\t\t\t\t\t\t\t" << " b/V: " << fbDens << endl << endl;
    }
  }

  cout << "\t\t\t " << fConstraintInfo << endl << endl;
  cout << "  **************************************************";
  cout << "****************************" << endl << endl;
}

//__________________________________________________________________________
TTMParameterSetBSQ& TTMParameterSetBSQ::operator=(const TTMParameterSetBSQ& obj)
{
  if (this == &obj) return *this;

  fB2Q = obj.GetB2Q();
  fSDens = obj.GetSDens();
  fCDens = obj.GetCDens();
  fbDens = obj.GetbDens();
  fConstraintInfo = obj.GetConstraintInfo();
  fMuSConstrain = obj.GetMuSConstrain();
  fMuQConstrain = obj.GetMuQConstrain();
  fMuCConstrain = obj.GetMuCConstrain();
  fMubConstrain = obj.GetMubConstrain();
  fPar = fParArray;

  fParArray[0] = obj.GetParameter(0);
  fParArray[1] = obj.GetParameter(1);
  fParArray[2] = obj.GetParameter(2);
  fParArray[3] = obj.GetParameter(3);
  fParArray[4] = obj.GetParameter(4);
  fParArray[5] = obj.GetParameter(5);
  fParArray[6] = obj.GetParameter(6);
  fParArray[7] = obj.GetParameter(7);
  fParArray[8] = obj.GetParameter(8);
  fParArray[9] = obj.GetParameter(9);

  return *this;
}
