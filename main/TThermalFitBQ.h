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

#ifndef Thermal_TTMParameterSetBQ
#include <TTMParameterSetBQ.h>
#endif

#ifndef Thermal_TTMThermalModelBQ
#include <TThermalModelBQ.h>
#endif

#ifndef Thermal_TTMThermalFitBQ
#define Thermal_TTMThermalFitBQ

using namespace std;

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalFitBQ							//
//            	                                 			//
// Fit class applying the canonical treatment of strangeness with 	//
// B and Q treated grand canonically					// 
//									//
//////////////////////////////////////////////////////////////////////////

class TTMThermalFitBQ:public TTMThermalFit {

 protected:
  
  TTMParameterSetBQ *fParm;
  Bool_t fWidth;
  Bool_t fNonStrangeQStats;
  TTMThermalModelBQ *fModel;

  TTMThermalModelBQ* GenerateThermalModel(TTMParticleSet *set);

 public:

  TTMThermalFitBQ();
  TTMThermalFitBQ(TTMParticleSet *set, TTMParameterSetBQ *par, char *file); 
   
  ~TTMThermalFitBQ();

  TTMParameterSetBQ* GetParameterSet(){return fParm;}
   
  void SetWidth(Bool_t x){fWidth = x;}
  void SetNonStrangeQStats(Bool_t x){fNonStrangeQStats = x;}

  Bool_t GetWidth(){return fWidth;}
  Bool_t GetNonStrangeQStats(){return fNonStrangeQStats;}

  ClassDef(TTMThermalFitBQ,1) // Class for fitting in the strangeness-canonical formalism 

};

#endif
