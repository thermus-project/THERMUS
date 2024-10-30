/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>
#include <iostream>

using namespace std;

#ifndef Thermal_TTMParameter
#define Thermal_TTMParameter

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMParameter                  		                        //
//                                                                  	//
// A Parameter Object 				     	            	//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMParameter:public TNamed {

 private:

  Double_t fValue;   // value of variable
  Double_t fError;   // error in variable
  Int_t fFlag;       // -1:constrained ; 0:fitted ; 1:fixed ; 2:uninitialised
  TString fStatus;   // reflects intended treatment or action taken  
  Double_t fStart;   // initial value used in fit
  Double_t fMin;     // lower bound of fit range
  Double_t fMax;     // upper bound of fit range	
  Double_t fStep;    // step size used in fit

 public:

  TTMParameter();
  TTMParameter(TString name, Double_t value, Double_t error = 0.);
  ~TTMParameter() {}
   
  void SetParameter(TString name, Double_t value, Double_t error = 0.);

  void SetTMName(TString x) {fName = x;}
  void SetValue(Double_t x) {fValue = x;}
  void SetError(Double_t x) {fError = x;}
  void SetStatus(TString x) {fStatus = x;}

  void Constrain();
  void Fit(Double_t start, Double_t min, Double_t max, Double_t step); 
  void Fix(Double_t value, Double_t error = 0.);
  
  Double_t GetValue() const {return fValue;}
  Double_t GetError() const {return fError;}
  Int_t GetFlag() const {return fFlag;}
  Double_t GetStart() const {return fStart;}
  Double_t GetMin() const {return fMin;}
  Double_t GetMax() const {return fMax;}
  Double_t GetStep() const {return fStep;}
  TString GetStatus() const {return fStatus;}

  void List();

  TTMParameter& operator=(const TTMParameter& obj);

  ClassDef(TTMParameter,1) // A Parameter Object 

};

#endif
