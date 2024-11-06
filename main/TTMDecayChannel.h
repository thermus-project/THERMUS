/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>
#include <fstream>
#include <iostream>

using namespace std;

#ifndef Thermal_TTMDecayChannel
#define Thermal_TTMDecayChannel

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMDecayChannel                                                  	//
//                                                                  	//
// Stores daughter id's and Branching ratio in an object for		// 
// storage in a ROOT container class				 	//
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMDecayChannel:public TObject {

 private:

  Double_t fBRatio;            // Branching ratio (fraction i.e. NOT %!)
  TList *fDaughters;	       // list of daughter id's

 public:

  TTMDecayChannel(); 
  TTMDecayChannel(Double_t fraction, TList *list);
  ~TTMDecayChannel() { }
   
  Double_t GetBRatio() const {return fBRatio;}
  TList* GetDaughterList() const {return fDaughters;}

  void SetBRatio(Double_t x) {fBRatio = x;}
  void SetDaughterList(TList *x) {fDaughters = x;}

  void List();

  TTMDecayChannel& operator=(const TTMDecayChannel& obj);
 
  ClassDef(TTMDecayChannel,1) // Decay channel object for storage in a ROOT container class 

};

#endif
