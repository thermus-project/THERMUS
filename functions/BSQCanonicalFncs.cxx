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

// Author: Spencer Wheaton 25 March 2005        //

#include <TROOT.h>
#include <TObject.h>
#include <TMath.h>

using namespace std;

Double_t Norm(Double_t *x, Double_t *mpar){

  return 1. / pow(2. * TMath::Pi(), 2.) *
    exp(2. * (mpar[0] * (cos(x[1])-1.) + 
              mpar[1] * (cos(x[0] + x[1])-1.) + mpar[2] * (cos(x[0])-1.)));
}

Double_t FcnAbsOmega(Double_t *x, Double_t *bpar){
	    
  Double_t OmegaRe = bpar[0] * cos(x[1]) +
    bpar[1] +
    bpar[2] * cos(-x[0]) +
    bpar[3] * cos(-x[0] + x[1]) +
    bpar[4] * cos(-x[0] - x[1]) +
    bpar[5] * cos(-x[1]) +
    bpar[6] * cos(2. * x[1]) +
    bpar[7] * cos(-2. * x[0] - x[1]) +
    bpar[8] * cos(-2 * x[0]) +
    bpar[9] * cos(-3 * x[0] - x[1]);
	    
  Double_t OmegaIm = bpar[0] * sin(x[1]) +
    bpar[2] * sin(-x[0]) +
    bpar[3] * sin(-x[0] + x[1]) +
    bpar[4] * sin(-x[0] - x[1]) +
    bpar[5] * sin(-x[1]) +
    bpar[6] * sin(2. * x[1]) +
    bpar[7] * sin(-2. * x[0] - x[1]) +
    bpar[8] * sin(-2 * x[0]) +
    bpar[9] * sin(-3 * x[0] - x[1]);

  return TMath::Sqrt(pow(OmegaRe, 2) + pow(OmegaIm, 2));

}

Double_t FcnArgOmega(Double_t *x, Double_t *bpar){
	    
  Double_t OmegaRe = bpar[0] * cos(x[1]) +
    bpar[1] +
    bpar[2] * cos(-x[0]) +
    bpar[3] * cos(-x[0] + x[1]) +
    bpar[4] * cos(-x[0] - x[1]) +
    bpar[5] * cos(-x[1]) +
    bpar[6] * cos(2. * x[1]) +
    bpar[7] * cos(-2. * x[0] - x[1]) +
    bpar[8] * cos(-2 * x[0]) +
    bpar[9] * cos(-3 * x[0] - x[1]);
	    
  Double_t OmegaIm = bpar[0] * sin(x[1]) +
    bpar[2] * sin(-x[0]) +
    bpar[3] * sin(-x[0] + x[1]) +
    bpar[4] * sin(-x[0] - x[1]) +
    bpar[5] * sin(-x[1]) +
    bpar[6] * sin(2. * x[1]) +
    bpar[7] * sin(-2. * x[0] - x[1]) +
    bpar[8] * sin(-2 * x[0]) +
    bpar[9] * sin(-3 * x[0] - x[1]);

  Double_t arg;

  if (OmegaIm == 0. && OmegaRe == 0.){
    arg = 0.;
  }else{
    arg = TMath::ATan2(OmegaIm, OmegaRe);
  }

  return arg;

}

Double_t Partition(Double_t *x, Double_t *par){

  Double_t mpar[3];
  for(Int_t i = 0 ; i < 3 ; i++){
    mpar[i] = par[3+i];
  }

  Double_t bpar[10];
  for(Int_t j = 0 ; j < 10 ; j++){
    bpar[j] = par[6+j];
  }

  Double_t ArgOmega = FcnArgOmega(x,bpar);
  Double_t AbsOmega = FcnAbsOmega(x,bpar);
  Double_t cfactor = Norm(x,mpar);

  return cfactor * TMath::BesselI((Int_t)TMath::Abs(par[0]),2.*AbsOmega)*
    cos(par[1]*x[0] + par[2]*x[1] - par[0]*ArgOmega);

}

Double_t CorrPiplus(Double_t *x, Double_t *par){

  Double_t OldPar2 = par[2];
  par[2] = (par[2]-1.);

  Double_t result = Partition(x,par);

  par[2] = OldPar2;

  return result;

}


Double_t CorrPiminus(Double_t *x, Double_t *par){

  Double_t OldPar2 = par[2];
  par[2] = (par[2]+1.);

  Double_t result = Partition(x,par);

  par[2] = OldPar2;

  return result;

}
	    
Double_t CorrKminus(Double_t *x, Double_t *par){

  Double_t OldPar1 = par[1];
  Double_t OldPar2 = par[2];
  par[1] = (par[1]+1.);
  par[2] = (par[2]+1.);

  Double_t result = Partition(x,par);

  par[1] = OldPar1;
  par[2] = OldPar2;

  return result;

}
	    
Double_t CorrKplus(Double_t *x, Double_t *par){

  Double_t OldPar1 = par[1];
  Double_t OldPar2 = par[2];
  par[1] = (par[1]-1.);
  par[2] = (par[2]-1.);

  Double_t result = Partition(x,par);

  par[1] = OldPar1;
  par[2] = OldPar2;

  return result;

}
	    
Double_t CorrKzero(Double_t *x, Double_t *par){

  Double_t OldPar1 = par[1];
  par[1] = (par[1]-1.);

  Double_t result = Partition(x,par);

  par[1] = OldPar1;

  return result;

}
	    
Double_t CorrAKzero(Double_t *x, Double_t *par){

  Double_t OldPar1 = par[1];
  par[1] = (par[1]+1.);

  Double_t result = Partition(x,par);

  par[1] = OldPar1;

  return result;

}
	    
Double_t CorrProton(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]-1.);
  par[2] = (par[2]-1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[2] = OldPar2;

  return result;

}
	    
Double_t CorrAProton(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]+1.);
  par[2] = (par[2]+1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[2] = OldPar2;

  return result;

}
	    
Double_t CorrNeutron(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  par[0] = (par[0]-1.);
    
  Double_t result = Partition(x,par);

  par[0] = OldPar0;

  return result;

}
	   
Double_t CorrANeutron(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  par[0] = (par[0]+1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;

  return result;

}
	   
Double_t CorrLa(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar1 = par[1];
  par[0] = (par[0]-1.);
  par[1] = (par[1]+1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[1] = OldPar1;

  return result;

}
	    
Double_t CorrALa(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar1 = par[1];
  par[0] = (par[0]+1.);
  par[1] = (par[1]-1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[1] = OldPar1;

  return result;

}
	    
Double_t CorrSigmaplus(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar1 = par[1];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]-1.);
  par[1] = (par[1]+1.);
  par[2] = (par[2]-1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[1] = OldPar1;
  par[2] = OldPar2;

  return result;

}
	    
Double_t CorrASigmaplus(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar1 = par[1];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]+1.);
  par[1] = (par[1]-1.); 
  par[2] = (par[2]+1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[1] = OldPar1;
  par[2] = OldPar2;

  return result;

}
	   
Double_t CorrSigmaminus(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar1 = par[1];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]-1.);
  par[1] = (par[1]+1.);
  par[2] = (par[2]+1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[1] = OldPar1;
  par[2] = OldPar2;

  return result;

}
	   
Double_t CorrASigmaminus(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar1 = par[1];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]+1.);
  par[1] = (par[1]-1.);
  par[2] = (par[2]-1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[1] = OldPar1;
  par[2] = OldPar2;

  return result;

}
	   
Double_t CorrDeltaminus(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]-1.);
  par[2] = (par[2]+1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[2] = OldPar2;

  return result;

}
	    
Double_t CorrADeltaminus(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]+1.);
  par[2] = (par[2]-1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[2] = OldPar2;

  return result;

}
	    
Double_t CorrDeltaplusplus(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]-1.);
  par[2] = (par[2]-2.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[2] = OldPar2;

  return result;

}
	    
Double_t CorrADeltaplusplus(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]+1.);
  par[2] = (par[2]+2.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[2] = OldPar2;

  return result;

}
	    
Double_t CorrKsiminus(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar1 = par[1];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]-1.);
  par[1] = (par[1]+2.); 
  par[2] = (par[2]+1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[1] = OldPar1;
  par[2] = OldPar2;

  return result;

}
	    
Double_t CorrAKsiminus(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar1 = par[1];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]+1.);
  par[1] = (par[1]-2.);
  par[2] = (par[2]-1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[1] = OldPar1;
  par[2] = OldPar2;

  return result;

}
	    
Double_t CorrKsi0(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar1 = par[1];
  par[0] = (par[0]-1.);
  par[1] = (par[1]+2.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[1] = OldPar1;

  return result;

}
	    
Double_t CorrAKsi0(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar1 = par[1];
  par[0] = (par[0]+1.);
  par[1] = (par[1]-2.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[1] = OldPar1;

  return result;

}
	    
Double_t CorrOmega(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar1 = par[1];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]-1.);
  par[1] = (par[1]+3.); 
  par[2] = (par[2]+1.);

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[1] = OldPar1;
  par[2] = OldPar2;

  return result;

}

 
Double_t CorrAOmega(Double_t *x, Double_t *par){

  Double_t OldPar0 = par[0];
  Double_t OldPar1 = par[1];
  Double_t OldPar2 = par[2];
  par[0] = (par[0]+1.);
  par[1] = (par[1]-3.); 
  par[2] = (par[2]-1.); 

  Double_t result = Partition(x,par);

  par[0] = OldPar0;
  par[1] = OldPar1;
  par[2] = OldPar2;

  return result;

}	   
