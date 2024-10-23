/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>
#include <iostream>

using namespace std;

#ifndef Thermal_TTMIDObj
#define Thermal_TTMIDObj

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMIDObj                                          			//
//                                                                  	//
// Object containing ID for storage in container class. 		// 
//                                                                  	//
//////////////////////////////////////////////////////////////////////////

class TTMIDObj:public TObject {

 private:

  Int_t fID;                   // particle ID

 public:

  TTMIDObj();
  TTMIDObj(Int_t x);
  ~TTMIDObj() { } 
   
  void SetID(Int_t x);

  Int_t GetID() const {return fID;}
   
  TTMIDObj& operator=(const TTMIDObj& obj);

  ClassDef(TTMIDObj,1) // ID Object for storage in ROOT container class

};

#endif
