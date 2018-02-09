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
// The parameter set to be used when treating B and Q grand-canonically and
// S canonically. B/2Q may be used to constrain muQ when this parameter
// set is linked to a particle set in a ThermalModelBQ object
//

#include <TTMParameterSetBQ.h>

ClassImp(TTMParameterSetBQ)

//__________________________________________________________________________
TTMParameterSetBQ::TTMParameterSetBQ(Double_t temp, Double_t mub,
                                       Double_t muq, Double_t gs,
                                       Double_t can_r, Double_t r,
                                       Double_t b2q, Double_t S, Double_t temp_error,
                                       Double_t mub_error,
                                       Double_t muq_error, Double_t gs_error,
                                       Double_t can_r_error,
                                       Double_t r_error)
{
  // Sets all parameters and their errors as well as B/2Q.
  // All parameters are set as "fixed type".
  //

  fB2Q = b2q;
  fS = S;

  fPar = fParArray;

  fParArray[0].SetParameter("T", temp, temp_error);
  fParArray[1].SetParameter("muB", mub, mub_error);
  fParArray[2].SetParameter("muQ", muq, muq_error);
  fMuQConstrain = false;
  fParArray[3].SetParameter("gammas", gs, gs_error);
  fParArray[4].SetParameter("Can. radius", can_r, can_r_error);
  fCorrRConstrain = false;
  fParArray[5].SetParameter("radius", r, r_error);
  fConstraintInfo = "Parameters unconstrained";
}

//__________________________________________________________________________
TTMParameterSetBQ::TTMParameterSetBQ()
{
  // Sets parameter names and values and errors to 0.
  //

  fB2Q = 0.;
  fS = 0.;

  fPar = fParArray;

  fParArray[0].SetParameter("T", 0., 0.);
  fParArray[1].SetParameter("muB", 0., 0.);
  fParArray[2].SetParameter("muQ", 0., 0.);
  fMuQConstrain = false;
  fParArray[3].SetParameter("gammas", 0., 0.);
  fParArray[4].SetParameter("Can. radius", 0., 0.);
  fCorrRConstrain = false;
  fParArray[5].SetParameter("radius", 0., 0.);
  fConstraintInfo = "Parameters unconstrained";
}

//__________________________________________________________________________
void TTMParameterSetBQ::ConserveSGlobally()
{
  // Changes the correlation radius to a constrained type parameter. 
  // 	
   
  fCorrRConstrain = true;
  fParArray[4].Constrain();
}

//__________________________________________________________________________
void TTMParameterSetBQ::ConstrainMuQ(Double_t x)
{
  // Changes muQ to a constrained type parameter. x is the initial B/2Q 
  // ratio in the system.
  // 	
   
  fB2Q = x;
  fMuQConstrain = true;
  fParArray[2].Constrain();
}

//__________________________________________________________________________
void TTMParameterSetBQ::FitT(Double_t start, Double_t min,
                             Double_t max, Double_t step)
{
  fParArray[0].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBQ::FitMuB(Double_t start, Double_t min,
                               Double_t max, Double_t step)
{
  fParArray[1].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBQ::FitMuQ(Double_t start, Double_t min,
                               Double_t max, Double_t step)
{
  fMuQConstrain = false;
  fParArray[2].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBQ::FitGammas(Double_t start, Double_t min,
                                  Double_t max, Double_t step)
{
  fParArray[3].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBQ::FitCanRadius(Double_t start, Double_t min,
                                     Double_t max, Double_t step)
{
  fCorrRConstrain = false;
  fParArray[4].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBQ::FitRadius(Double_t start, Double_t min,
                                  Double_t max, Double_t step)
{
  fParArray[5].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetBQ::FixT(Double_t value, Double_t error)
{
  fParArray[0].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBQ::FixMuB(Double_t value, Double_t error)
{
  fParArray[1].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBQ::FixMuQ(Double_t value, Double_t error)
{
  fMuQConstrain = false;
  fParArray[2].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBQ::FixGammas(Double_t value, Double_t error)
{
  fParArray[3].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBQ::FixCanRadius(Double_t value, Double_t error)
{
  fCorrRConstrain = false;
  fParArray[4].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBQ::FixRadius(Double_t value, Double_t error)
{
  fParArray[5].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetBQ::List()
{
  cout << "  ***************************** Thermal Parameters ";
  cout << "**************************** " << endl << endl;
  cout << "\t\t" << " Strangeness inside Canonical Volume = " << fS << endl << endl;

  for (Int_t i = 0; i <= 5; i++) {
    fParArray[i].List();
    if (i == 2 && fMuQConstrain) {
      cout << "\t\t\t\t\t\t\t\t" << " B/2Q: " << fB2Q << endl << endl;
    }

  }

  cout << "\t\t\t " << fConstraintInfo << endl << endl;
  cout << "  ******************************************************";
  cout << "************************" << endl << endl;
}

//__________________________________________________________________________
TTMParameterSetBQ& TTMParameterSetBQ::operator=(const TTMParameterSetBQ& obj)
{
  if (this == &obj) return *this; 

  fB2Q = obj.GetB2Q();
  fS = obj.GetS();
  fConstraintInfo = obj.GetConstraintInfo();
  fMuQConstrain = obj.GetMuQConstrain();
  fCorrRConstrain = obj.GetCorrRConstrain();  
  fParArray[0] = obj.GetParameter(0);
  fParArray[1] = obj.GetParameter(1);
  fParArray[2] = obj.GetParameter(2);
  fParArray[3] = obj.GetParameter(3);
  fParArray[4] = obj.GetParameter(4);
  fParArray[5] = obj.GetParameter(5);
  return *this;
}
