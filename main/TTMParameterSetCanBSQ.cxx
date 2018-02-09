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
// The parameter set to be used when treating B, S and Q canonically. 
//

#include <TTMParameterSetCanBSQ.h>

ClassImp(TTMParameterSetCanBSQ)

//__________________________________________________________________________
TTMParameterSetCanBSQ::TTMParameterSetCanBSQ(Double_t temp, Int_t b,
                                               Double_t s, Double_t q, Double_t gs,
                                               Double_t r,
                                               Double_t temp_error,
                                               Double_t b_error, Double_t s_error,
                                               Double_t q_error, Double_t gs_error,
                                               Double_t r_error)
{
  // Sets all parameters and their errors.
  // All parameters are set as "fixed type".
  //

  fPar = fParArray;

  fParArray[0].SetParameter("T", temp, temp_error);
  fParArray[1].SetParameter("B", b, b_error);
  fParArray[2].SetParameter("S", s, s_error);
  fParArray[3].SetParameter("Q", q, q_error);
  fParArray[4].SetParameter("gammas", gs, gs_error);
  fParArray[5].SetParameter("radius", r, r_error);
  fConstraintInfo = "";
}

//__________________________________________________________________________
TTMParameterSetCanBSQ::TTMParameterSetCanBSQ()
{
  // Sets parameter names and values and errors to 0.
  //

  fPar = fParArray;

  fParArray[0].SetParameter("T", 0., 0.);
  fParArray[1].SetParameter("B", 0, 0);
  fParArray[2].SetParameter("S", 0, 0);
  fParArray[3].SetParameter("Q", 0, 0);
  fParArray[4].SetParameter("gammas", 0., 0.);
  fParArray[5].SetParameter("radius", 0., 0.);
}

//__________________________________________________________________________
void TTMParameterSetCanBSQ::FitT(Double_t start, Double_t min,
                                 Double_t max, Double_t step)
{
  fParArray[0].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetCanBSQ::FitB(Int_t start, Int_t min,
                                 Int_t max, Int_t step)
{
  fParArray[1].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetCanBSQ::FitS(Double_t start, Double_t min,
                                 Double_t max, Double_t step)
{
  fParArray[2].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetCanBSQ::FitQ(Double_t start, Double_t min,
                                 Double_t max, Double_t step)
{
  fParArray[3].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetCanBSQ::FitGammas(Double_t start, Double_t min,
                                      Double_t max, Double_t step)
{
  fParArray[4].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetCanBSQ::FitRadius(Double_t start, Double_t min,
                                      Double_t max, Double_t step)
{
  fParArray[5].Fit(start, min, max, step);
}

//__________________________________________________________________________
void TTMParameterSetCanBSQ::FixT(Double_t value, Double_t error)
{
  fParArray[0].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetCanBSQ::FixB(Int_t value, Double_t error)
{
  fParArray[1].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetCanBSQ::FixS(Double_t value, Double_t error)
{
  fParArray[2].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetCanBSQ::FixQ(Double_t value, Double_t error)
{
  fParArray[3].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetCanBSQ::FixGammas(Double_t value, Double_t error)
{
  fParArray[4].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetCanBSQ::FixRadius(Double_t value, Double_t error)
{
  fParArray[5].Fix(value, error);
}

//__________________________________________________________________________
void TTMParameterSetCanBSQ::List()
{
  cout << "  ***************************** Thermal Parameters ";
  cout << "**************************** " << endl << endl;
  for (Int_t i = 0; i <= 5; i++) {
    fParArray[i].List();
  }

  cout << "\t\t\t " << fConstraintInfo << endl << endl;
  cout << "  ******************************************************";
  cout << "************************" << endl << endl;
}

//__________________________________________________________________________
TTMParameterSetCanBSQ& TTMParameterSetCanBSQ::operator=(const TTMParameterSetCanBSQ& obj)
{
  if (this == &obj) return *this;

  fConstraintInfo = obj.GetConstraintInfo();
  fParArray[0] = obj.GetParameter(0);
  fParArray[1] = obj.GetParameter(1);
  fParArray[2] = obj.GetParameter(2);
  fParArray[3] = obj.GetParameter(3);
  fParArray[4] = obj.GetParameter(4);
  fParArray[5] = obj.GetParameter(5);
  return *this;
}
