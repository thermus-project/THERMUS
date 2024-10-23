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
//
// TTMParameter stores the name, value and error of a variable as well as the
// "variable type" (i.e. whether it is to be fixed, constrained or fit). 
// When of "fit-type", the fit parameters can be set  
//

#include <TTMParameter.h>

ClassImp(TTMParameter)

//__________________________________________________________________________
TTMParameter::TTMParameter()
{
  fName = "unnamed";
  fValue = 0.;
  fError = 0.;
  fFlag = 2;
  fStatus = "(UNINITIALISED)";
  fStart = 0.;
  fMin = 0.;
  fMax = 0.;
  fStep = 0.;
}

//__________________________________________________________________________
TTMParameter::TTMParameter(TString name, Double_t value, Double_t error)
{
  // Simple constructor: sets the name, value and error of a variable,
  // making it of "FIXED" type. Use Constrain() or Fit(...) to change 
  // its type

  SetParameter(name, value, error);
}

//__________________________________________________________________________
void TTMParameter::SetParameter(TString name, Double_t value,
                                Double_t error)
{
  // Sets the name, value and error of a variable, making it of 
  // "FIXED" type. Use Constrain() or Fit(...) to change its   
  // type

  fName = name;
  fValue = value;
  fError = error;
  fFlag = 1;

  fStatus = "(FIXED)";

  fStart = 0.;
  fMin = 0.;
  fMax = 0.;
  fStep = 0.;
}

//__________________________________________________________________________
void TTMParameter::Constrain()
{
  // Changes the variable type to "CONSTRAIN" and sets error to 0.
	
  fFlag = -1;
  fStatus = "(to be CONSTRAINED)";
  fError = 0.;
}

//__________________________________________________________________________
void TTMParameter::Fit(Double_t start, Double_t min, Double_t max,
                       Double_t step)
{
  // Changes the variable type to "FIT", sets error to 0 and sets
  // the fit parameters:
  //		start: initial value for fit  
  //  		min & max : bounds for fit 
  //		step	: step size

  fFlag = 0;
  fStatus = "(to be FITTED)";
  fStart = start;
  fMin = min;
  fMax = max;
  fStep = step;
  fValue = start;
  fError = 0.;
}

//__________________________________________________________________________
void TTMParameter::Fix(Double_t value, Double_t error)
{
  // Changes the variable type to "FIXED" and sets its value and error.
     
  fFlag = 1;
  fStatus = "(FIXED)";
  fValue = value;
  fError = error;
}

//__________________________________________________________________________
void TTMParameter::List()
{
  // Outputs the variable's properties	
	
  cout.width(8);
  cout << fName << "\t = \t ";
  cout.width(7);
  cout << fValue;
  if (fError != 0.) {
    cout.width(7);
    cout << "\t +- \t";
    cout.width(10);
    cout << fError;
  } else {
    cout.width(7);
    cout << "\t    \t";
    cout.width(10);
    cout << " ";
  }
  if (fFlag == 0) {
    cout << "\t" << fStatus << endl;
    cout << "\t\t\t\t\t\t\t\t" << " start: " << fStart << endl;
    cout << "\t\t\t\t\t\t\t\t" << " range: " << fMin << " -- " << fMax <<
      endl;
    cout << "\t\t\t\t\t\t\t\t" << " step:  " << fStep << endl << endl;
  } else {
    cout << "\t" << fStatus << endl << endl;
  }
}

//__________________________________________________________________________
TTMParameter& TTMParameter::operator=(const TTMParameter& obj)
{
  if (this == &obj) return *this;

  fName = obj.GetName();
  fValue = obj.GetValue();
  fError = obj.GetError();
  fFlag = obj.GetFlag();
  fStatus = obj.GetStatus();
  fStart = obj.GetStart();
  fMin = obj.GetMin();
  fMax = obj.GetMax();
  fStep = obj.GetStep();
  return *this;
}
