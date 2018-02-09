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

// Author: Spencer Wheaton 14 July 2004         //

class TTMParticle;
class TTMThermalModel;
class TTMThermalModelBSQ;
class TTMThermalFit;
class TF1;
class TF2;

Double_t Norm(Double_t *x, Double_t *par);

Double_t Partition(Double_t *x, Double_t *par);
Double_t CorrPiplus(Double_t *x, Double_t *par);
Double_t CorrPiminus(Double_t *x, Double_t *par);
Double_t CorrKminus(Double_t *x, Double_t *par);
Double_t CorrKplus(Double_t *x, Double_t *par);
Double_t CorrKzero(Double_t *x, Double_t *par);
Double_t CorrAKzero(Double_t *x, Double_t *par);
Double_t CorrProton(Double_t *x, Double_t *par);
Double_t CorrAProton(Double_t *x, Double_t *par);
Double_t CorrNeutron(Double_t *x, Double_t *par);
Double_t CorrANeutron(Double_t *x, Double_t *par);
Double_t CorrLa(Double_t *x, Double_t *par);
Double_t CorrALa(Double_t *x, Double_t *par);
Double_t CorrSigmaplus(Double_t *x, Double_t *par);
Double_t CorrASigmaplus(Double_t *x, Double_t *par);
Double_t CorrSigmaminus(Double_t *x, Double_t *par);
Double_t CorrASigmaminus(Double_t *x, Double_t *par);
Double_t CorrDeltaminus(Double_t *x, Double_t *par);
Double_t CorrADeltaminus(Double_t *x, Double_t *par);
Double_t CorrDeltaplusplus(Double_t *x, Double_t *par);
Double_t CorrADeltaplusplus(Double_t *x, Double_t *par);
Double_t CorrKsiminus(Double_t *x, Double_t *par);
Double_t CorrAKsiminus(Double_t *x, Double_t *par);
Double_t CorrKsi0(Double_t *x, Double_t *par);
Double_t CorrAKsi0(Double_t *x, Double_t *par);
Double_t CorrOmega(Double_t *x, Double_t *par);
Double_t CorrAOmega(Double_t *x, Double_t *par);

Double_t FcnAbsOmega(Double_t *x, Double_t *par);
Double_t FcnArgOmega(Double_t *x, Double_t *par);

Double_t FcnDistr(Double_t *x, Double_t *par);
Double_t FcnDens(Double_t *x, Double_t *par);
Double_t FcnEnergyDens(Double_t *x, Double_t *par);
Double_t FcnEntropyDens(Double_t *x, Double_t *par);
Double_t FcnPressure(Double_t *x, Double_t *par);
Double_t FcnDensWidth(Double_t *x, Double_t *par);
Double_t FcnEnergyDensWidth(Double_t *x, Double_t *par);
Double_t FcnEntropyDensWidth(Double_t *x, Double_t *par);
Double_t FcnPressureWidth(Double_t *x, Double_t *par);
Double_t FcnDensNormWidth(Double_t *x, Double_t *par);
Double_t FcnDensBoltzmannWidth(Double_t *x, Double_t *par);
Double_t FcnEnergyBoltzmannWidth(Double_t *x, Double_t *par);
Double_t FcnEntropyBoltzmannWidth(Double_t *x, Double_t *par);

Double_t IntegrateLegendre32(TF1* func, Double_t a, Double_t b);
Double_t IntegrateLegendre20(TF1* func, Double_t a, Double_t b);
Double_t IntegrateLegendre40(TF1* func, Double_t a, Double_t b);

Double_t IntegrateLaguerre32(TF1* func);
Double_t Integrate2DLaguerre32Legendre32(TF2* func, Double_t ay, Double_t by);

Double_t FindExclVolPressure(TTMThermalModelBSQ *model, Double_t limit);

void Minuit_fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, 
		Int_t iflag);
void fit_function(TTMThermalModel *mod,char *file,Int_t flag = 0);

Double_t FcnEnergyFlucBoltzmannWidth(Double_t *x, Double_t *par);
Double_t GenericQStatDeriv(Double_t *x, Double_t *par);
Double_t GenericQStatDerivWidth(Double_t *x, Double_t *par);
