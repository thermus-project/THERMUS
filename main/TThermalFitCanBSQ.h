/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <fstream>
#include <iostream>

#ifndef Thermal_TTMThermalFit
#include <TThermalFit.h>
#endif

#ifndef Thermal_TTMParameterSetCanBSQ
#include <TTMParameterSetCanBSQ.h>
#endif

#ifndef Thermal_TTMThermalModelCanBSQ
#include <TThermalModelCanBSQ.h>
#endif

#ifndef Thermal_TTMThermalFitCanBSQ
#define Thermal_TTMThermalFitCanBSQ

using namespace std;

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalFitCanBSQ							//
//            	                                 			//
// Fit class applying the canonical formalism	 			// 
//									//
//////////////////////////////////////////////////////////////////////////

class TTMThermalFitCanBSQ:public TTMThermalFit {

 protected:
  
  TTMParameterSetCanBSQ *fParm;	
  Bool_t fWidth;
  TTMThermalModelCanBSQ *fModel;	 	

  TTMThermalModelCanBSQ* GenerateThermalModel(TTMParticleSet *set);

 public:

  TTMThermalFitCanBSQ();
  TTMThermalFitCanBSQ(TTMParticleSet *set, TTMParameterSetCanBSQ *par, char *file);
   
  ~TTMThermalFitCanBSQ();

  TTMParameterSetCanBSQ* GetParameterSet(){return fParm;}

  void SetParameterSet(TTMParameterSetCanBSQ *par){fParm = par;}
  void SetWidth(Bool_t x){fWidth = x;}

  Bool_t GetWidth(){return fWidth;}

  ClassDef(TTMThermalFitCanBSQ,1) // Canonical fit object

};

#endif
