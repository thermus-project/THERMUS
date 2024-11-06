/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>
#include <TMath.h>

#ifndef Thermal_TTMParameter
#include <TTMParameter.h>
#endif

#ifndef Thermal_TTMParameterSet
#define Thermal_TTMParameterSet

//////////////////////////////////////////////////////////////////////////
//                                                                   	//
// TTMParameterSet                                           	     	//
//                                                                   	//
// Base class for parameter set objects. All derived classes must 	//
// define Double_t GetRadius().  					//
//                                                                   	//
//////////////////////////////////////////////////////////////////////////

class TTMParameterSet:public TObject {

 protected:

  TTMParameter *fPar;	       //!Pointer to first Parameter Object
  TString fConstraintInfo;     // Constraint information

 public:

  TTMParameterSet();
  ~TTMParameterSet() { }
 
  TString GetConstraintInfo() const {return fConstraintInfo;}
  TTMParameter* GetParameter(Int_t i) {return (fPar+i);}
  TTMParameter GetParameter(Int_t i) const {return *(fPar+i);}

  void SetConstraintInfo(TString x) {fConstraintInfo = x;}

  virtual Double_t GetRadius() const = 0 ;
   
  Double_t GetVolume() const {return 4. * TMath::Pi() / 3. * pow(GetRadius(), 3);}

  ClassDef(TTMParameterSet, 1) // Base class for parameter set objects.

};

#endif
