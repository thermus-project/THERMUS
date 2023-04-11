/* Copyright(c) 2004-2018, THERMUS Project,        All rights reserved. *
 * See cxx source for full Copyright notice                             */

//Authors: Natasha Sharma and Jean Cleymans 2 February 2021 (modification of the corresponding
//BQ file, TThermalFitBQ.h,  by Spencer Wheaton 14 July 2004)

#include <TROOT.h>
#include <TObject.h>
#include <fstream>
#include <iostream>

#ifndef Thermal_TTMThermalFit
#include <TThermalFit.h>
#endif

#ifndef Thermal_TTMParameterSetSQ
#include <TTMParameterSetSQ.h>
#endif

#ifndef Thermal_TTMThermalModelSQ
#include <TThermalModelSQ.h>
#endif

#ifndef Thermal_TTMThermalFitSQ
#define Thermal_TTMThermalFitSQ

using namespace std;

//////////////////////////////////////////////////////////////////////////
//                                                                 	//
// TTMThermalFitSQ							//
//            	                                 			//
// Fit class applying the canonical treatment of strangeness with 	//
// B and Q treated grand canonically					// 
//									//
//////////////////////////////////////////////////////////////////////////

class TTMThermalFitSQ:public TTMThermalFit {

 protected:
  
  TTMParameterSetSQ *fParm;
  Bool_t fWidth;
  Bool_t fNonBaryonQStats;
  TTMThermalModelSQ *fModel;

  TTMThermalModelSQ* GenerateThermalModel(TTMParticleSet *set);

 public:

  TTMThermalFitSQ();
  TTMThermalFitSQ(TTMParticleSet *set, TTMParameterSetSQ *par, char *file); 
   
  ~TTMThermalFitSQ();

  TTMParameterSetSQ* GetParameterSet(){return fParm;}
   
  void SetWidth(Bool_t x){fWidth = x;}
  void SetNonStrangeQStats(Bool_t x){fNonBaryonQStats = x;}

  Bool_t GetWidth(){return fWidth;}
  Bool_t GetNonStrangeQStats(){return fNonBaryonQStats;}

  ClassDef(TTMThermalFitSQ,1) // Class for fitting in the strangeness-canonical formalism 

};

#endif
