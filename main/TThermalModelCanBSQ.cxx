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
// Canonical thermal model class (exact B, S & Q conservation).
//  

#include <TThermalModelCanBSQ.h>
#include <FncsThermalModel.h>
#include <FncsConstrain.h>
#include <TH2F.h>
#include <TF2.h>

ClassImp(TTMThermalModelCanBSQ)

//__________________________________________________________________________
TTMThermalModelCanBSQ::TTMThermalModelCanBSQ(TTMParticleSet *particles,
                                               TTMParameterSetCanBSQ 
                                               *parameters, Bool_t width)
{
  fDescriptor = "BSQCanonical";
  fPartSet = particles;
  fParm = parameters;
  fWidth = width;
  fDensTable->Delete();
  fSplus = fSminus = fBplus = fBminus = fQplus = fQminus = fCplus = fCminus = fbplus = fbminus = 0.;
  fStrange = fBaryon = fCharge = fDensity = fWroblewski = fCharm = fBeauty = 0.;
  fEnergy = fEntropy = fPressure = 0.;

  fListCorrFactorOutput = false;
  fTF2IntCorrFactorCalc = false;

  for(Int_t i = 0 ; i < 17 ; i++){
    par[i] = 0.;
  }

  fTol = 1.e-20;
  fRadius = TMath::Pi();
  fMin = 10;
  fStepScale = 0.5;
  fOrder = 8;
}

//__________________________________________________________________________
TTMThermalModelCanBSQ::TTMThermalModelCanBSQ()
{
  fDescriptor = "BSQCanonical";
  fPartSet = (TTMParticleSet *) 0;
  fParm = (TTMParameterSetCanBSQ *) 0;
  fWidth = true;
  fDensTable->Delete();
  fSplus = fSminus = fBplus = fBminus = fQplus = fQminus = fCplus = fCminus = fbplus = fbminus = 0.;
  fStrange = fBaryon = fCharge = fDensity = fWroblewski = fCharm = fBeauty = 0.;
  fEnergy = fEntropy = fPressure = 0.;

  fListCorrFactorOutput = false;
  fTF2IntCorrFactorCalc = false;

  for(Int_t i = 0 ; i < 17 ; i++){
    par[i] = 0.;
  }

  fTol = 1.e-20;
  fRadius = TMath::Pi();
  fMin = 10;
  fStepScale = 0.5;
  fOrder = 8;
}

//__________________________________________________________________________
void TTMThermalModelCanBSQ::UpdatePartitionFcnParameters()
{

  for(Int_t i = 0 ; i < 17 ; i++){
    par[i] = 0.;
  }

  TIter next(GetParticleTable());
  TTMParticle *part;

  Double_t Z_pi0 = 0.;
  Double_t Z_pip = 0.;
  Double_t Z_kp = 0.;
  Double_t Z_k0 = 0.;
  Double_t Z_proton = 0.;
  Double_t Z_neutron = 0.;
  Double_t Z_lambda = 0.;
  Double_t Z_sigmap = 0.;
  Double_t Z_sigmam = 0.;
  Double_t Z_deltam = 0.;
  Double_t Z_deltapp = 0.;
  Double_t Z_ksim = 0.;
  Double_t Z_ksi0 = 0.;
  Double_t Z_omega = 0.;
   
  while ((part = (TTMParticle *) next())) {
      
    TTMThermalParticleCanBSQ *ptr = new TTMThermalParticleCanBSQ(part, fParm, 1.);

    Double_t GCPartDens;      // uncorrected densities (GC)

    if(part->GetWidth() == 0. || !fWidth){
      GCPartDens = ptr->DensityBoltzmannNoWidth();
    }else{
      GCPartDens = ptr->DensityBoltzmannWidth();     
    }

    delete ptr;      
    
    Double_t volume = fParm->GetVolume();
    
    Double_t Z = GCPartDens * volume;
    Double_t B = part->GetB();
    Double_t S = part->GetS();
    Double_t Q = part->GetQ();
    
    if (S == 0 && B == 0 && Q == 0) {  
      Z_pi0 += Z;
    } else if (S == 0 && B == 0 && Q == +1) {
      Z_pip += Z;
    } else if (S == +1 && B == 0 && Q == +1) {
      Z_kp += Z;
    } else if (S == 1 && B == 0 && Q == 0) {
      Z_k0 += Z;
    } else if (S == 0 && B == 1 && Q == +1) { 
      Z_proton += Z;
    } else if (S == 0 && B == 1 && Q == 0) { 
      Z_neutron += Z;
    } else if (S == -1 && B == 1 && Q == 0) {
      Z_lambda += Z;
    } else if (S == -1 && B == 1 && Q == +1) {
      Z_sigmap += Z;
    } else if (S == -1 && B == 1 && Q == -1) {
      Z_sigmam += Z;
    } else if (S == 0 && B == 1 && Q == -1) {
      Z_deltam += Z;
    } else if (S == 0 && B == 1 && Q == +2) {
      Z_deltapp += Z;
    } else if (S == -2 && B == 1 && Q == -1) {
      Z_ksim += Z;
    } else if (S == -2 && B == 1 && Q == 0) {
      Z_ksi0 += Z;
    } else if (S == -3 && B == 1 && Q == -1) {
      Z_omega += Z; 
    }

  }

  Int_t Btot = fParm->GetB();
  Double_t Stot = fParm->GetS();
  Double_t Qtot = fParm->GetQ();

  par[0] = Btot;
  par[1] = Stot;
  par[2] = Qtot;
  par[3] = Z_pip;
  par[4] = Z_kp;
  par[5] = Z_k0;
  par[6] = Z_proton;
  par[7] = Z_neutron;
  par[8] = Z_lambda;
  par[9] = Z_sigmap;
  par[10] = Z_sigmam;
  par[11] = Z_deltam;
  par[12] = Z_deltapp;
  par[13] = Z_ksim;
  par[14] = Z_ksi0;
  par[15] = Z_omega;
  par[16] = Z_pi0;

}

//__________________________________________________________________________
Int_t TTMThermalModelCanBSQ::PrimPartDens()
{
  // Calculates the primordial particle densities and populates the density
  // hash table. The parameters are not constrained first!. 
  // This is the function used by GenerateParticleDens().
  //

  fDensTable->Delete();
  TIter next(GetParticleTable());
  TTMParticle *part;

  fSplus = fSminus = fBplus = fBminus = fQplus = fQminus = fCplus = fCminus = fbplus = fbminus = 0.;
  fStrange = fBaryon = fCharge = fDensity = fWroblewski = fCharm = fBeauty = 0.;
  fEnergy = fEntropy = fPressure = 0.;

  UpdatePartitionFcnParameters();

  Double_t xg4[4] = { -0.861136311594053, -0.339981043584856,  0.339981043584856,  0.861136311594053};
  
  Double_t wg4[4] = {  0.347854845137454,  0.652145154862546,  0.652145154862546,  0.347854845137454};
  
  Double_t xg8[8] = { -0.960289856497536, -0.796666477413627, -0.525532409916329, -0.183434642495650,
                      0.183434642495650,  0.525532409916329,  0.796666477413627,  0.960289856497536};

  Double_t wg8[8] = {  0.101228536290376,  0.222381034453374,  0.313706645877887,  0.362683783378362,
                       0.362683783378362,  0.313706645877887,  0.222381034453374,  0.101228536290376};

  Double_t xg16[16] = { -0.989400934991649, -0.944575023073232, -0.865631202387831, -0.755404408355003,
                        -0.617876244402643, -0.458016777657227, -0.281603550779258, -0.095012509837637,
                        0.095012509837637,  0.281603550779258,  0.458016777657227,  0.617876244402643,
                        0.755404408355003,  0.865631202387831,  0.944575023073232,  0.989400934991649};

  Double_t wg16[16] = {  0.027152459411754,  0.062253523938647,  0.095158511682492,  0.124628971255533,
                         0.149595988816576,  0.169156519395002,  0.182603415044923,  0.189450610455068,
                         0.189450610455068,  0.182603415044923,  0.169156519395002,  0.149595988816576,
                         0.124628971255533,  0.095158511682492,  0.062253523938647,  0.027152459411754};
   
  Double_t xg20[20] = { -0.993128599185095, -0.963971927277914, -0.912234428251326, -0.839116971822219, -0.746331906460151, 
			-0.636053680726515, -0.510867001950827, -0.373706088715420, -0.227785851141645, -0.076526521133497, 
                        0.076526521133497,  0.227785851141645,  0.373706088715420,  0.510867001950827,  0.636053680726515,
                        0.746331906460151,  0.839116971822219,  0.912234428251326,  0.963971927277914,  0.993128599185095};

  Double_t wg20[20] = {  0.017614007139152,  0.040601429800387,  0.062672048334109,  0.083276741576705,  0.101930119817240, 
			 0.118194531961518,  0.131688638449177,  0.142096109318382,  0.149172986472603,  0.152753387130726, 
			 0.152753387130726,  0.149172986472604,  0.142096109318382,  0.131688638449177,  0.118194531961518,
			 0.101930119817240,  0.083276741576705,  0.062672048334109,  0.040601429800387,  0.017614007139152};

  Double_t xg32[32] = { -0.997263861849482, -0.985611511545268, -0.964762255587506, -0.934906075937740, 
			-0.896321155766052, -0.849367613732570, -0.794483795967942, -0.732182118740290, 
			-0.663044266930215, -0.587715757240762, -0.506899908932229, -0.421351276130635, 
			-0.331868602282128, -0.239287362252137, -0.144471961582796, -0.048307665687738,
                        0.048307665687738,  0.144471961582796,  0.239287362282127,  0.331868602282128,  
                        0.421351276130635,  0.506899908932229,  0.587715757240762,  0.663044266930215,  
                        0.732182118740290,  0.794483795967942,  0.849367613732570,  0.896321155766052,
                        0.934906075937740,  0.964762255587506,  0.985611511545268,  0.997263861849482};

  Double_t wg32[32] = {  0.007018610009470,  0.016274394730906,  0.025392065309262,  0.034273862913021,
                         0.042835898022227,  0.050998059262376,  0.058684093478536,  0.065822222776362,
                         0.072345794108849,  0.078193895787070,  0.083311924226947,  0.087652093004404,
                         0.091173878695764,  0.093844399080805,  0.095638720079275,  0.096540088514728,
			 0.096540088514728,  0.095638720079275,  0.093844399080805,  0.091173878695764,
                         0.087652093004404,  0.083311924226947,  0.078193895787070,  0.072345794108849,
                         0.065822222776362,  0.058684093478536,  0.050998059262376,  0.042835898022227,
                         0.034273862913021,  0.025392065309262,  0.016274394730906,  0.007018610009470};

  const Int_t poly = fOrder;

  Double_t xg[poly];
  Double_t wg[poly];

  for(Int_t a=0; a<poly; a++){

    if(poly==4) { 
      xg[a]=xg4[a];  
      wg[a]=wg4[a]; 
    }else if(poly==8) { 
      xg[a]=xg8[a];  
      wg[a]=wg8[a]; 
    }else if(poly==16){ 
      xg[a]=xg16[a]; 
      wg[a]=wg16[a];
    }else if(poly==20){ 
      xg[a]=xg20[a]; 
      wg[a]=wg20[a];
    }else if(poly==32){ 
      xg[a]=xg32[a]; 
      wg[a]=wg32[a];
    }

  }

  Double_t fac_all = 0.;
  Double_t fac_piplus = 0.;
  Double_t fac_piminus = 0.;
  Double_t fac_kplus = 0.;
  Double_t fac_kminus = 0.;
  Double_t fac_kzero = 0.;
  Double_t fac_akzero = 0.;
  Double_t fac_la = 0.;
  Double_t fac_ala = 0.;
  Double_t fac_sigmaplus = 0.;
  Double_t fac_asigmaplus = 0.;
  Double_t fac_sigmaminus = 0.;
  Double_t fac_asigmaminus = 0.;
  Double_t fac_proton = 0.;
  Double_t fac_aproton = 0.;
  Double_t fac_neutron = 0.;
  Double_t fac_aneutron = 0.;
  Double_t fac_deltaminus = 0.;
  Double_t fac_adeltaminus = 0.;
  Double_t fac_deltaplusplus = 0.;
  Double_t fac_adeltaplusplus = 0.;
  Double_t fac_ksiminus = 0.;
  Double_t fac_aksiminus = 0.;
  Double_t fac_ksi0 = 0.;
  Double_t fac_aksi0 = 0.;
  Double_t fac_Omega = 0.;
  Double_t fac_aOmega = 0.;

  // *********************************************************************** //

  Double_t afac_all = 0.;
  Double_t afac_piplus = 0.;
  Double_t afac_piminus = 0.;
  Double_t afac_kplus = 0.;
  Double_t afac_kminus = 0.;
  Double_t afac_kzero = 0.;
  Double_t afac_akzero = 0.;
  Double_t afac_la = 0.;
  Double_t afac_ala = 0.;
  Double_t afac_sigmaplus = 0.;
  Double_t afac_asigmaplus = 0.;
  Double_t afac_sigmaminus = 0.;
  Double_t afac_asigmaminus = 0.;
  Double_t afac_proton = 0.;
  Double_t afac_aproton = 0.;
  Double_t afac_neutron = 0.;
  Double_t afac_aneutron = 0.;
  Double_t afac_deltaminus = 0.;
  Double_t afac_adeltaminus = 0.;
  Double_t afac_deltaplusplus = 0.;
  Double_t afac_adeltaplusplus = 0.;
  Double_t afac_ksiminus = 0.;
  Double_t afac_aksiminus = 0.;
  Double_t afac_ksi0 = 0.;
  Double_t afac_aksi0 = 0.;
  Double_t afac_Omega = 0.;
  Double_t afac_aOmega = 0.;

  TF2 *fPartition = new TF2("Partition",Partition,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fPartition->SetParameters(par);

  TF2 *fPiplus = new TF2("Piplus",CorrPiplus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fPiplus->SetParameters(par);

  TF2 *fPiminus = new TF2("Piminus",CorrPiminus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fPiminus->SetParameters(par);

  TF2 *fKminus = new TF2("Kminus",CorrKminus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fKminus->SetParameters(par);

  TF2 *fKplus = new TF2("Kplus",CorrKplus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fKplus->SetParameters(par);

  TF2 *fKzero = new TF2("Kzero",CorrKzero,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fKzero->SetParameters(par);

  TF2 *fAKzero = new TF2("AKzero",CorrAKzero,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fAKzero->SetParameters(par);

  TF2 *fProton = new TF2("Proton",CorrProton,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fProton->SetParameters(par);

  TF2 *fAProton = new TF2("AProton",CorrAProton,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fAProton->SetParameters(par);

  TF2 *fNeutron = new TF2("Neutron",CorrNeutron,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fNeutron->SetParameters(par);

  TF2 *fANeutron = new TF2("ANeutron",CorrANeutron,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fANeutron->SetParameters(par);

  TF2 *fLa = new TF2("La",CorrLa,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fLa->SetParameters(par);

  TF2 *fALa = new TF2("ALa",CorrALa,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fALa->SetParameters(par);

  TF2 *fSigmaplus = new TF2("Sigmaplus",CorrSigmaplus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fSigmaplus->SetParameters(par);

  TF2 *fASigmaplus = new TF2("ASigmaplus",CorrASigmaplus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fASigmaplus->SetParameters(par);

  TF2 *fSigmaminus = new TF2("Sigmaminus",CorrSigmaminus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fSigmaminus->SetParameters(par);

  TF2 *fASigmaminus = new TF2("ASigmaminus",CorrASigmaminus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fASigmaminus->SetParameters(par);

  TF2 *fDeltaminus = new TF2("Deltaminus",CorrDeltaminus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fDeltaminus->SetParameters(par);

  TF2 *fADeltaminus = new TF2("ADeltaminus",CorrADeltaminus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fADeltaminus->SetParameters(par);

  TF2 *fDeltaplusplus = new TF2("Deltaplusplus",CorrDeltaplusplus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fDeltaplusplus->SetParameters(par);

  TF2 *fADeltaplusplus = new TF2("ADeltaplusplus",CorrADeltaplusplus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fADeltaplusplus->SetParameters(par);

  TF2 *fKsiminus = new TF2("Ksiminus",CorrKsiminus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fKsiminus->SetParameters(par);

  TF2 *fAKsiminus = new TF2("AKsiminus",CorrAKsiminus,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fAKsiminus->SetParameters(par);

  TF2 *fKsi0 = new TF2("Ksi0",CorrKsi0,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fKsi0->SetParameters(par);

  TF2 *fAKsi0 = new TF2("AKsi0",CorrAKsi0,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fAKsi0->SetParameters(par);

  TF2 *fOmega = new TF2("Omega",CorrOmega,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fOmega->SetParameters(par);

  TF2 *fAOmega = new TF2("AOmega",CorrAOmega,-2.*TMath::Pi(),2.*TMath::Pi(),-2.*TMath::Pi(),2.*TMath::Pi(),16);
  fAOmega->SetParameters(par);

  // *********************************************************************** //

  Double_t mpar[3];
  for(Int_t i = 0 ; i < 3 ; i++){
    mpar[i] = par[3+i];
  }

  Double_t bpar[10];
  for(Int_t j = 0 ; j < 10 ; j++){
    bpar[j] = par[6+j];
  }

  Int_t Btot = fParm->GetB();
  Double_t Stot = fParm->GetS();
  Double_t Qtot = fParm->GetQ();

  Double_t maxlB = 2.*(TMath::Abs(Btot)+1);
  Double_t maxlS = 2.*(TMath::Abs(Stot)+1);
  Double_t maxlQ = 2.*(TMath::Abs(Qtot)+1); 

  Int_t min = (Int_t) ceil( fStepScale * sqrt( pow(maxlB,2.) + pow(maxlS,2.) + pow(maxlQ,2.) ) );
  Int_t max = TMath::Max(min,fMin);

  Double_t step = TMath::Pi() / max;
  Double_t aQ, bQ, aS, bS;   		  // limits of integration
  Double_t phiQ, phiS;        		  // variables of integration
  
  Int_t count = 0;
  Int_t discard = 0;
  
  Bool_t STOP = false;
  Double_t TESTZ = 0.;

  for (Int_t lS = 0; lS < 2 * max; lS++) { // start of phiS integral

    if(lS<max){ 
      aS = lS * step;
    }else{ 
      aS = TMath::Pi() - (lS+1) * step;
    }

    bS = aS + step;
    
    STOP = false;
   
    for (Int_t lQ = 0; lQ < max; lQ++) { // start of phiQ integral

      aQ = lQ * step;

      bQ = aQ + step;

      Double_t radius = sqrt( pow(aS,2.) + pow(aQ,2.) );
     
      if( (radius<=fRadius) || (STOP==false) ) {

	count++;
        Double_t dpart = 0.;
        Double_t dpipl = 0.;
	Double_t dpim = 0.;
	Double_t dkm = 0.;
	Double_t dkp = 0.;
	Double_t dk0 = 0.;
	Double_t dak0 = 0.;
	Double_t dp = 0.;
	Double_t dap = 0.;
	Double_t dn = 0.;
	Double_t dan = 0.;
	Double_t dla = 0.;
	Double_t dala = 0.;
	Double_t dsp = 0.;
	Double_t dasp = 0.;
	Double_t dsm = 0.;
	Double_t dasm = 0.;
	Double_t ddm = 0.;
	Double_t dadm = 0.;
	Double_t ddpp = 0.;
	Double_t dadpp = 0.;
	Double_t dksim = 0.;
	Double_t daksim = 0.;
	Double_t dksi0 = 0.;
	Double_t daksi0 = 0.;
	Double_t dom = 0.;
	Double_t daom = 0.;

        Double_t xint[2];

	for (Int_t k = 0; k < poly; k++) {

	  phiS = (bS + aS) / 2. + step / 2. * xg[k];
	  
          xint[0] = phiS;

	  for (Int_t j = 0; j < poly; j++) {

	    phiQ = (bQ + aQ) / 2. + step / 2. * xg[j];

	    xint[1] = phiQ;

	    Double_t AbsOmega, ArgOmega;
	    
	    AbsOmega = FcnAbsOmega(xint,bpar);
	    ArgOmega = FcnArgOmega(xint,bpar);

	    Double_t intfactor = 2. *
	      wg[k] * wg[j] * step / 2. * step / 2.;

	    TESTZ = intfactor * Partition(xint,par);

            if(!TMath::IsNaN(TESTZ)){
              dpart += TESTZ;	
            }	    

	    Double_t term;

            term = intfactor * CorrPiplus(xint,par);

            if(!TMath::IsNaN(term)){
              dpipl += term;
            }

	    term = intfactor * CorrPiminus(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dpim += term;
            }
	    
	    term = intfactor * CorrKminus(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dkm += term;
            }

	    term = intfactor * CorrKplus(xint,par); 
	    
	    if(!TMath::IsNaN(term)){	    
              dkp += term;
            }
	    
	    term = intfactor * CorrKzero(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dk0 += term;
            }
	    
	    term = intfactor * CorrAKzero(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dak0 += term;
            }
	    
	    term = intfactor * CorrProton(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dp += term;
            }
	    
	    term = intfactor * CorrAProton(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dap += term;
            }
	    
	    term = intfactor * CorrNeutron(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dn += term;
            }
	   
	    term = intfactor * CorrANeutron(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dan += term;
            }
	   
	    term = intfactor * CorrLa(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dla += term;
            }

	    term = intfactor * CorrALa(xint,par);
	    
	    if(!TMath::IsNaN(term)){	    
              dala += term;
            }
	    
	    term = intfactor * CorrSigmaplus(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dsp += term;
            }
	    
	    term = intfactor * CorrASigmaplus(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dasp += term;
            }
	   
	    term = intfactor * CorrSigmaminus(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dsm += term;
            }
	   
	    term = intfactor * CorrASigmaminus(xint,par); 

	    if(!TMath::IsNaN(term)){	    
              dasm += term;
            }
	   
	    term = intfactor * CorrDeltaminus(xint,par);

	    if(!TMath::IsNaN(term)){	    
              ddm += term;
            }
	    
	    term = intfactor * CorrADeltaminus(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dadm += term;
            }
	    
	    term = intfactor * CorrDeltaplusplus(xint,par);

	    if(!TMath::IsNaN(term)){	    
              ddpp += term;
            }
	    
	    term = intfactor * CorrADeltaplusplus(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dadpp += term;
            }
	    
	    term = intfactor * CorrKsiminus(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dksim += term;
            }
	    
	    term = intfactor * CorrAKsiminus(xint,par);

	    if(!TMath::IsNaN(term)){	    
              daksim += term;
            }
	    
	    term = intfactor * CorrKsi0(xint,par);

	    if(!TMath::IsNaN(term)){	    
              dksi0 += term;
            }
	    
	    term = intfactor * CorrAKsi0(xint,par);

	    if(!TMath::IsNaN(term)){	    
              daksi0 += term;
            }

	    term = intfactor * CorrOmega(xint,par);
	    
	    if(!TMath::IsNaN(term)){	    
              dom += term;
            }

	    term = intfactor * CorrAOmega(xint,par);

	    if(!TMath::IsNaN(term)){ 
              daom += term;
            }	   
 
	  }

	}

        fac_all += dpart;
        fac_piplus += dpipl;
        fac_piminus += dpim;
        fac_kminus += dkm;
        fac_kplus += dkp;
        fac_kzero += dk0;
        fac_akzero += dak0;
        fac_proton += dp;
        fac_aproton += dap;
        fac_neutron += dn;
        fac_aneutron += dan;
        fac_la += dla;
        fac_ala += dala;
        fac_sigmaplus += dsp;
        fac_asigmaplus += dasp;
        fac_sigmaminus += dsm;
        fac_asigmaminus += dasm;
        fac_deltaminus += ddm;
        fac_adeltaminus += dadm;
        fac_deltaplusplus += ddpp;
        fac_adeltaplusplus += dadpp;
        fac_ksiminus += dksim;
        fac_aksiminus += daksim;
        fac_ksi0 += dksi0;
        fac_aksi0 += daksi0;
        fac_Omega += dom;
        fac_aOmega += daom;

        // *********************************************************************** //

        if(fTF2IntCorrFactorCalc){

          Double_t eps = 1e-6;

          afac_all += 2.*fPartition->Integral(aS,bS,aQ,bQ,dpart/2.*eps);
          afac_piplus += 2.*fPiplus->Integral(aS,bS,aQ,bQ,dpipl/2.*eps);
          afac_piminus += 2.*fPiminus->Integral(aS,bS,aQ,bQ,dpim/2.*eps);
          afac_kplus += 2.*fKplus->Integral(aS,bS,aQ,bQ,dkp/2.*eps);
          afac_kminus += 2.*fKminus->Integral(aS,bS,aQ,bQ,dkm/2.*eps);
          afac_kzero += 2.*fKzero->Integral(aS,bS,aQ,bQ,dk0/2.*eps);
          afac_akzero += 2.*fAKzero->Integral(aS,bS,aQ,bQ,dak0/2.*eps);
          afac_proton += 2.*fProton->Integral(aS,bS,aQ,bQ,dp/2.*eps);
          afac_aproton += 2.*fAProton->Integral(aS,bS,aQ,bQ,dap/2.*eps);
          afac_neutron += 2.*fNeutron->Integral(aS,bS,aQ,bQ,dn/2.*eps);
          afac_aneutron += 2.*fANeutron->Integral(aS,bS,aQ,bQ,dan/2.*eps);
          afac_la += 2.*fLa->Integral(aS,bS,aQ,bQ,dla/2.*eps);
          afac_ala += 2.*fALa->Integral(aS,bS,aQ,bQ,dala/2.*eps);
          afac_sigmaplus += 2.*fSigmaplus->Integral(aS,bS,aQ,bQ,dsp/2.*eps);
          afac_asigmaplus += 2.*fASigmaplus->Integral(aS,bS,aQ,bQ,dasp/2.*eps);
          afac_sigmaminus += 2.*fSigmaminus->Integral(aS,bS,aQ,bQ,dsm/2.*eps);
          afac_asigmaminus += 2.*fASigmaminus->Integral(aS,bS,aQ,bQ,dasm/2.*eps);
          afac_deltaminus += 2.*fDeltaminus->Integral(aS,bS,aQ,bQ,ddm/2.*eps);
          afac_adeltaminus += 2.*fADeltaminus->Integral(aS,bS,aQ,bQ,dadm/2.*eps);
          afac_deltaplusplus += 2.*fDeltaplusplus->Integral(aS,bS,aQ,bQ,ddpp/2.*eps);
          afac_adeltaplusplus += 2.*fADeltaplusplus->Integral(aS,bS,aQ,bQ,dadpp/2.*eps);
          afac_ksiminus += 2.*fKsiminus->Integral(aS,bS,aQ,bQ,dksim/2.*eps);
          afac_aksiminus += 2.*fAKsiminus->Integral(aS,bS,aQ,bQ,daksim/2.*eps);
          afac_ksi0 += 2.*fKsi0->Integral(aS,bS,aQ,bQ,dksi0/2.*eps);
          afac_aksi0 += 2.*fAKsi0->Integral(aS,bS,aQ,bQ,daksi0/2.*eps);
          afac_Omega += 2.*fOmega->Integral(aS,bS,aQ,bQ,dom/2.*eps);
          afac_aOmega += 2.*fAOmega->Integral(aS,bS,aQ,bQ,daom/2.*eps);

        }

        // *********************************************************************** //

	if(TMath::Abs( TESTZ/fac_all ) < fTol){ 
          STOP = true; 
        }else{ 
          STOP = false; 
        }

      }else{

        discard++;

      }

    }

  }

  if(fListCorrFactorOutput){

    cout<<"************************ Gauss-Legendre Method ************************"<<endl;
    cout<<"fac_all: "<<fac_all<<endl;
    cout<<"fac_piplus: "<<fac_piplus<<endl;
    cout<<"fac_piminus: "<<fac_piminus<<endl;
    cout<<"fac_kplus: "<<fac_kplus<<endl;
    cout<<"fac_kminus: "<<fac_kminus<<endl;
    cout<<"fac_kzero: "<<fac_kzero<<endl;
    cout<<"fac_akzero: "<<fac_akzero<<endl;
    cout<<"fac_la: "<<fac_la<<endl;
    cout<<"fac_ala: "<<fac_ala<<endl;
    cout<<"fac_sigmaplus: "<<fac_sigmaplus<<endl;
    cout<<"fac_asigmaplus: "<<fac_asigmaplus<<endl;
    cout<<"fac_sigmaminus: "<<fac_sigmaminus<<endl;
    cout<<"fac_asigmaminus: "<<fac_asigmaminus<<endl;
    cout<<"fac_proton: "<<fac_proton<<endl;
    cout<<"fac_aproton: "<<fac_aproton<<endl;
    cout<<"fac_neutron: "<<fac_neutron<<endl;
    cout<<"fac_aneutron: "<<fac_aneutron<<endl;
    cout<<"fac_deltaminus: "<<fac_deltaminus<<endl;
    cout<<"fac_adeltaminus: "<<fac_adeltaminus<<endl;
    cout<<"fac_deltaplusplus: "<<fac_deltaplusplus<<endl;
    cout<<"fac_adeltaplusplus: "<<fac_adeltaplusplus<<endl;
    cout<<"fac_ksiminus: "<<fac_ksiminus<<endl;
    cout<<"fac_aksiminus: "<<fac_aksiminus<<endl;
    cout<<"fac_ksi0: "<<fac_ksi0<<endl;
    cout<<"fac_aksi0: "<<fac_aksi0<<endl;
    cout<<"fac_Omega: "<<fac_Omega<<endl;
    cout<<"fac_aOmega: "<<fac_aOmega<<endl;

    if(fTF2IntCorrFactorCalc){

      cout<<"************************ TF2 Method ************************"<<endl;
      cout<<"afac_all: "<<afac_all<<endl;
      cout<<"afac_piplus: "<<afac_piplus<<endl;
      cout<<"afac_piminus: "<<afac_piminus<<endl;
      cout<<"afac_kplus: "<<afac_kplus<<endl;
      cout<<"afac_kminus: "<<afac_kminus<<endl;
      cout<<"afac_kzero: "<<afac_kzero<<endl;
      cout<<"afac_akzero: "<<afac_akzero<<endl;
      cout<<"afac_la: "<<afac_la<<endl;
      cout<<"afac_ala: "<<afac_ala<<endl;
      cout<<"afac_sigmaplus: "<<afac_sigmaplus<<endl;
      cout<<"afac_asigmaplus: "<<afac_asigmaplus<<endl;
      cout<<"afac_sigmaminus: "<<afac_sigmaminus<<endl;
      cout<<"afac_asigmaminus: "<<afac_asigmaminus<<endl;
      cout<<"afac_proton: "<<afac_proton<<endl;
      cout<<"afac_aproton: "<<afac_aproton<<endl;
      cout<<"afac_neutron: "<<afac_neutron<<endl;
      cout<<"afac_aneutron: "<<afac_aneutron<<endl;
      cout<<"afac_deltaminus: "<<afac_deltaminus<<endl;
      cout<<"afac_adeltaminus: "<<afac_adeltaminus<<endl;
      cout<<"afac_deltaplusplus: "<<afac_deltaplusplus<<endl;
      cout<<"afac_adeltaplusplus: "<<afac_adeltaplusplus<<endl;
      cout<<"afac_ksiminus: "<<afac_ksiminus<<endl;
      cout<<"afac_aksiminus: "<<afac_aksiminus<<endl;
      cout<<"afac_ksi0: "<<afac_ksi0<<endl;
      cout<<"afac_aksi0: "<<afac_aksi0<<endl;
      cout<<"afac_Omega: "<<afac_Omega<<endl;
      cout<<"afac_aOmega: "<<afac_aOmega<<endl;
      cout<<"****************************************************** "<<endl;

    }

  }
 
  Int_t sign_check = 0;  

  if(fac_all <= 0.){
    cout<<"fac_all <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_piplus <= 0.){
    cout<<"fac_piplus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_piminus <= 0.){
    cout<<"fac_piminus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_kplus <= 0.){
    cout<<"fac_kplus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_kminus <= 0.){
    cout<<"fac_kminus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_kzero <= 0.){
    cout<<"fac_kzero <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_akzero <= 0.){
    cout<<"fac_akzero <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_la <= 0.){
    cout<<"fac_la <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_ala <= 0.){
    cout<<"fac_ala <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_sigmaplus <= 0.){
    cout<<"fac_sigmaplus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_asigmaplus <= 0.){
    cout<<"fac_asigmaplus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_sigmaminus <= 0.){
    cout<<"fac_sigmaminus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_asigmaminus <= 0.){
    cout<<"fac_asigmaminus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_proton <= 0.){
    cout<<"fac_proton <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_aproton <= 0.){
    cout<<"fac_aproton <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_neutron <= 0.){
    cout<<"fac_neutron <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_aneutron <= 0.){
    cout<<"fac_aneutron <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_deltaminus <= 0.){
    cout<<"fac_deltaminus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_adeltaminus <= 0.){
    cout<<"fac_adeltaminus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_deltaplusplus <= 0.){
    cout<<"fac_deltaplusplus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_adeltaplusplus <= 0.){
    cout<<"fac_adeltaplusplus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_ksiminus <= 0.){
    cout<<"fac_ksiminus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_aksiminus <= 0.){
    cout<<"fac_aksiminus <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_ksi0 <= 0.){
    cout<<"fac_ksi0 <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_aksi0 <= 0.){
    cout<<"fac_aksi0 <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_Omega <= 0.){
    cout<<"fac_Omega <= 0"<<endl;
    sign_check = 1;
  } 
  if(fac_aOmega <= 0.){
    cout<<"fac_aOmega <= 0"<<endl;
    sign_check = 1;
  } 

  if(sign_check){
    cout<<"Problems: Negative Correction Factors!"<<endl;
    return 1;
  }

  // *********************************************************************** //

  Int_t asign_check = 0;  

  if(fTF2IntCorrFactorCalc){

    if(afac_all <= 0.){
      cout<<"afac_all <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_piplus <= 0.){
      cout<<"afac_piplus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_piminus <= 0.){
      cout<<"afac_piminus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_kplus <= 0.){
      cout<<"afac_kplus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_kminus <= 0.){
      cout<<"afac_kminus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_kzero <= 0.){
      cout<<"afac_kzero <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_akzero <= 0.){
      cout<<"afac_akzero <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_la <= 0.){
      cout<<"afac_la <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_ala <= 0.){
      cout<<"afac_ala <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_sigmaplus <= 0.){
      cout<<"afac_sigmaplus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_asigmaplus <= 0.){
      cout<<"afac_asigmaplus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_sigmaminus <= 0.){
      cout<<"afac_sigmaminus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_asigmaminus <= 0.){
      cout<<"afac_asigmaminus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_proton <= 0.){
      cout<<"afac_proton <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_aproton <= 0.){
      cout<<"afac_aproton <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_neutron <= 0.){
      cout<<"afac_neutron <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_aneutron <= 0.){
      cout<<"afac_aneutron <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_deltaminus <= 0.){
      cout<<"afac_deltaminus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_adeltaminus <= 0.){
      cout<<"afac_adeltaminus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_deltaplusplus <= 0.){
      cout<<"afac_deltaplusplus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_adeltaplusplus <= 0.){
      cout<<"afac_adeltaplusplus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_ksiminus <= 0.){
      cout<<"afac_ksiminus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_aksiminus <= 0.){
      cout<<"afac_aksiminus <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_ksi0 <= 0.){
      cout<<"afac_ksi0 <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_aksi0 <= 0.){
      cout<<"afac_aksi0 <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_Omega <= 0.){
      cout<<"afac_Omega <= 0"<<endl;
      asign_check = 1;
    } 
    if(afac_aOmega <= 0.){
      cout<<"afac_aOmega <= 0"<<endl;
      asign_check = 1;
    } 

  }

  if(asign_check){
    cout<<"Problems with TF2 integrals: Negative Correction Factors!"<<endl;
    return 1;
  }

  // *********************************************************************** //

  Double_t lfac_all;
  Double_t lfac_piplus;
  Double_t lfac_piminus;
  Double_t lfac_kplus;
  Double_t lfac_kminus;
  Double_t lfac_kzero;
  Double_t lfac_akzero;
  Double_t lfac_la;
  Double_t lfac_ala;
  Double_t lfac_sigmaplus;
  Double_t lfac_asigmaplus;
  Double_t lfac_sigmaminus;
  Double_t lfac_asigmaminus;
  Double_t lfac_proton;
  Double_t lfac_aproton;
  Double_t lfac_neutron;
  Double_t lfac_aneutron;
  Double_t lfac_deltaminus;
  Double_t lfac_adeltaminus; 
  Double_t lfac_deltaplusplus;
  Double_t lfac_adeltaplusplus;
  Double_t lfac_ksiminus;
  Double_t lfac_aksiminus;
  Double_t lfac_ksi0;
  Double_t lfac_aksi0;
  Double_t lfac_Omega;
  Double_t lfac_aOmega;

  if(!fTF2IntCorrFactorCalc){

    lfac_all = log(fac_all);
    lfac_piplus = log(fac_piplus);
    lfac_piminus = log(fac_piminus);
    lfac_kplus = log(fac_kplus);
    lfac_kminus = log(fac_kminus);
    lfac_kzero = log(fac_kzero);
    lfac_akzero = log(fac_akzero);
    lfac_la = log(fac_la);
    lfac_ala = log(fac_ala);
    lfac_sigmaplus = log(fac_sigmaplus);
    lfac_asigmaplus = log(fac_asigmaplus);
    lfac_sigmaminus = log(fac_sigmaminus);
    lfac_asigmaminus = log(fac_asigmaminus);
    lfac_proton = log(fac_proton);
    lfac_aproton = log(fac_aproton);
    lfac_neutron = log(fac_neutron);
    lfac_aneutron = log(fac_aneutron);
    lfac_deltaminus = log(fac_deltaminus);
    lfac_adeltaminus = log(fac_adeltaminus); 
    lfac_deltaplusplus = log(fac_deltaplusplus);
    lfac_adeltaplusplus = log(fac_adeltaplusplus);
    lfac_ksiminus = log(fac_ksiminus);
    lfac_aksiminus = log(fac_aksiminus);
    lfac_ksi0 = log(fac_ksi0);
    lfac_aksi0 = log(fac_aksi0);
    lfac_Omega = log(fac_Omega);
    lfac_aOmega = log(fac_aOmega);

  }else{

    lfac_all = log(afac_all);
    lfac_piplus = log(afac_piplus);
    lfac_piminus = log(afac_piminus);
    lfac_kplus = log(afac_kplus);
    lfac_kminus = log(afac_kminus);
    lfac_kzero = log(afac_kzero);
    lfac_akzero = log(afac_akzero);
    lfac_la = log(afac_la);
    lfac_ala = log(afac_ala);
    lfac_sigmaplus = log(afac_sigmaplus);
    lfac_asigmaplus = log(afac_asigmaplus);
    lfac_sigmaminus = log(afac_sigmaminus);
    lfac_asigmaminus = log(afac_asigmaminus);
    lfac_proton = log(afac_proton);
    lfac_aproton = log(afac_aproton);
    lfac_neutron = log(afac_neutron);
    lfac_aneutron = log(afac_aneutron);
    lfac_deltaminus = log(afac_deltaminus);
    lfac_adeltaminus = log(afac_adeltaminus); 
    lfac_deltaplusplus = log(afac_deltaplusplus);
    lfac_adeltaplusplus = log(afac_adeltaplusplus);
    lfac_ksiminus = log(afac_ksiminus);
    lfac_aksiminus = log(afac_aksiminus);
    lfac_ksi0 = log(afac_ksi0);
    lfac_aksi0 = log(afac_aksi0);
    lfac_Omega = log(afac_Omega);
    lfac_aOmega = log(afac_aOmega);

  }

  // *********************************************************************** //

  fCorrpip = exp(lfac_piplus - lfac_all);
  fCorrpim = exp(lfac_piminus - lfac_all);
  fCorrkm = exp(lfac_kminus - lfac_all);
  fCorrkp = exp(lfac_kplus - lfac_all);
  fCorrk0 = exp(lfac_kzero - lfac_all);
  fCorrak0 = exp(lfac_akzero - lfac_all);
  fCorrproton = exp(lfac_proton - lfac_all);
  fCorraproton = exp(lfac_aproton - lfac_all);
  fCorrneutron = exp(lfac_neutron - lfac_all);
  fCorraneutron = exp(lfac_aneutron - lfac_all);
  fCorrlambda = exp(lfac_la - lfac_all);
  fCorralambda = exp(lfac_ala - lfac_all);
  fCorrsigmap = exp(lfac_sigmaplus - lfac_all);
  fCorrasigmap = exp(lfac_asigmaplus - lfac_all);
  fCorrsigmam = exp(lfac_sigmaminus - lfac_all);
  fCorrasigmam = exp(lfac_asigmaminus - lfac_all);
  fCorrdeltam = exp(lfac_deltaminus - lfac_all);
  fCorradeltam = exp(lfac_adeltaminus - lfac_all);
  fCorrdeltapp = exp(lfac_deltaplusplus - lfac_all);
  fCorradeltapp = exp(lfac_adeltaplusplus - lfac_all);
  fCorrksim = exp(lfac_ksiminus - lfac_all);
  fCorraksim = exp(lfac_aksiminus - lfac_all);
  fCorrksi0 = exp(lfac_ksi0 - lfac_all);
  fCorraksi0 = exp(lfac_aksi0 - lfac_all);
  fCorromega = exp(lfac_Omega - lfac_all);
  fCorraomega = exp(lfac_aOmega - lfac_all);

  Double_t Z_pip = par[3];
  Double_t Z_kp = par[4];
  Double_t Z_k0 = par[5];
  Double_t Z_pi0 = par[16];

  flnZtot = lfac_all + 2. * (Z_pip + Z_kp + Z_k0) + Z_pi0;

  fMuB = fParm->GetT() * log(fCorrneutron);
  fMuS = fParm->GetT() * log(fCorrk0);
  fMuQ = fParm->GetT() * log(fCorrpip);
  
  if(fListCorrFactorOutput){
    cout<<"lnZtot : "<<flnZtot<<endl;
    cout<<"fCorrpip: "<< fCorrpip<<endl;
    cout<<"fCorrpim: "<< fCorrpim<<endl;
    cout<<"fCorrkm: "<< fCorrkm<<endl;
    cout<<"fCorrkp: "<< fCorrkp<<endl;
    cout<<"fCorrk0: "<< fCorrk0<<endl;
    cout<<"fCorrak0: "<< fCorrak0<<endl;
    cout<<"fCorrproton: "<< fCorrproton<<endl;
    cout<<"fCorraproton: "<< fCorraproton<<endl;
    cout<<"fCorrneutron: "<< fCorrneutron<<endl;
    cout<<"fCorraneutron: "<< fCorraneutron<<endl;
    cout<<"fCorrlambda: "<< fCorrlambda<<endl;
    cout<<"fCorralambda: "<< fCorralambda<<endl;
    cout<<"fCorrsigmap: "<< fCorrsigmap<<endl;
    cout<<"fCorrasigmap: "<< fCorrasigmap<<endl;
    cout<<"fCorrsigmam: "<< fCorrsigmam<<endl;
    cout<<"fCorrasigmam: "<< fCorrasigmam<<endl;
    cout<<"fCorrdeltam: "<< fCorrdeltam<<endl;
    cout<<"fCorradeltam: "<< fCorradeltam<<endl;
    cout<<"fCorrdeltapp: "<< fCorrdeltapp<<endl;
    cout<<"fCorradeltapp: "<< fCorradeltapp<<endl;
    cout<<"fCorrksim: "<< fCorrksim<<endl;
    cout<<"fCorraksim: "<< fCorraksim<<endl;
    cout<<"fCorrksi0: "<< fCorrksi0<<endl;
    cout<<"fCorraksi0: "<< fCorraksi0<<endl;
    cout<<"fCorromega: "<< fCorromega<<endl;
    cout<<"fCorraomega: "<< fCorraomega<<endl;
  }
 
  delete fPartition;
  delete fPiplus;
  delete fPiminus;
  delete fKminus;
  delete fKplus;
  delete fKzero;
  delete fAKzero;
  delete fProton;
  delete fAProton;
  delete fNeutron;
  delete fANeutron;
  delete fLa;
  delete fALa;
  delete fSigmaplus;
  delete fASigmaplus;
  delete fSigmaminus;
  delete fASigmaminus;
  delete fDeltaminus;
  delete fADeltaminus;
  delete fDeltaplusplus;
  delete fADeltaplusplus;
  delete fKsiminus;
  delete fAKsiminus;
  delete fKsi0;
  delete fAKsi0;
  delete fOmega;
  delete fAOmega;

  next.Reset();
   
  while ((part = (TTMParticle *) next())) {
      
    Double_t CorrFactor = GetCorrFactor(part);	   
    TTMThermalParticleCanBSQ *ptr = new TTMThermalParticleCanBSQ(part, 
                                                                 fParm, CorrFactor);
    TTMDensObj *dens = new TTMDensObj(part->GetID());
      
    Double_t PartDens;

    if(part->GetWidth() == 0. || !fWidth){
      PartDens = ptr->DensityBoltzmannNoWidth();
    }else{
      PartDens = ptr->DensityBoltzmannWidth();
    }	      

    delete ptr;      
   
    dens->SetPrimDensity(PartDens);
    fDensTable->Add(dens);

    fDensity += PartDens;

    if(part->GetS() > 0.){
      fSplus += PartDens * (part->GetS());
    }else if(part->GetS() < 0.){
      fSminus += PartDens * (part->GetS()); 
    }

    if(part->GetB() > 0.){
      fBplus += PartDens * (part->GetB()); 
    }else if(part->GetB() < 0.){
      fBminus += PartDens * (part->GetB()); 
    }

    if(part->GetQ() > 0.){
      fQplus += PartDens * (part->GetQ());
    }else if(part->GetQ() < 0.){
      fQminus += PartDens * (part->GetQ());
    }

    if(part->GetCharm() > 0.){
      fCplus += PartDens * (part->GetCharm());
    }else if(part->GetCharm() < 0.){
      fCminus += PartDens * (part->GetCharm()); 
    }

    if(part->GetBeauty() > 0.){
      fbplus += PartDens * (part->GetBeauty());
    }else if(part->GetBeauty() < 0.){
      fbminus += PartDens * (part->GetBeauty()); 
    }
   // delete dens; // 09/02/2016 B.H. to be checked in order to avoid memory leaks
  }

  fStrange = fSplus + fSminus;
  fBaryon = fBplus + fBminus;
  fCharge = fQplus + fQminus;
  fCharm = fCplus + fCminus;
  fBeauty = fbplus + fbminus;

  if(fListCorrFactorOutput){

    cout<<"fSplus: "<<fSplus<<endl;
    cout<<"fSminus: "<<fSminus<<endl;
    cout<<"fBplus: "<<fBplus<<endl;
    cout<<"fBminus: "<<fBminus<<endl;
    cout<<"fQplus: "<<fQplus<<endl;
    cout<<"fQminus: "<<fQminus<<endl;
    cout<<"fCplus: "<<fCplus<<endl;
    cout<<"fCminus: "<<fCminus<<endl;
    cout<<"fbplus: "<<fbplus<<endl;
    cout<<"fbminus: "<<fbminus<<endl;

    cout<<"B:"<<fBaryon*fParm->GetVolume()<<endl;
    cout<<"S:"<<fStrange*fParm->GetVolume()<<endl;
    cout<<"Q:"<<fCharge*fParm->GetVolume()<<endl;
    cout<<"C:"<<fCharm*fParm->GetVolume()<<endl;
    cout<<"Beauty:"<<fBeauty*fParm->GetVolume()<<endl;

  }

  Double_t volume = fParm->GetVolume();

  Double_t errB;
  if(Btot!=0){
    errB = TMath::Abs(fBaryon*volume-Btot)/Btot;
  }else{
    errB = TMath::Abs(fBaryon*volume-Btot);
  }

  Double_t errS;
  if(Stot!=0){
    errS = TMath::Abs(fStrange*volume-Stot)/Stot;
  }else{
    errS = TMath::Abs(fStrange*volume-Stot);
  }

  Double_t errQ;
  if(Qtot!=0){
    errQ = TMath::Abs(fCharge*volume-Qtot)/Qtot;
  }else{
    errQ = TMath::Abs(fCharge*volume-Qtot);
  }

  if(errB>0.001||errS>0.001||errQ>0.001){
    cout<<"WARNING: Constraints not satisfied!!"<<endl;
    return 1;
  }else{
    return 0;
  }

}

//__________________________________________________________________________
void TTMThermalModelCanBSQ::PopulateZHistograms(TH2F *hZ[27])
{
  // Since the histograms are allocated off the heap they must be cleaned up 
  // afterwards
  //

  Int_t nphiS = 100;
  Int_t nphiQ = 100;

  hZ[0] = new TH2F("Partition","Partition",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[1] = new TH2F("Piplus","Piplus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[2] = new TH2F("Piminus","Piminus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[3] = new TH2F("Kminus","Kminus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[4] = new TH2F("Kplus","Kplus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[5] = new TH2F("Kzero","Kzero",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[6] = new TH2F("AKzero","AKzero",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[7] = new TH2F("Proton","Proton",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[8] = new TH2F("AProton","AProton",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[9] = new TH2F("Neutron","Neutron",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[10] = new TH2F("ANeutron","ANeutron",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[11] = new TH2F("La","La",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[12] = new TH2F("ALa","ALa",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[13] = new TH2F("Sigmaplus","Sigmaplus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[14] = new TH2F("ASigmaplus","ASigmaplus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[15] = new TH2F("Sigmaminus","Sigmaminus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[16] = new TH2F("ASigmaminus","ASigmaminus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[17] = new TH2F("Deltaminus","Deltaminus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[18] = new TH2F("ADeltaminus","ADeltaminus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[19] = new TH2F("Deltaplusplus","Deltaplusplus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[20] = new TH2F("ADeltaplusplus","ADeltaplusplus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[21] = new TH2F("Ksiminus","Ksiminus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[22] = new TH2F("AKsiminus","AKsiminus",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[23] = new TH2F("Ksi0","Ksi0",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[24] = new TH2F("AKsi0","AKsi0",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[25] = new TH2F("Omega","Omega",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());
  hZ[26] = new TH2F("AOmega","AOmega",nphiS,-TMath::Pi(),+TMath::Pi(),nphiQ,-TMath::Pi(),+TMath::Pi());

  for(Int_t i = 0 ; i < 27 ; i++){
    hZ[i]->SetStats(kFALSE);
    hZ[i]->GetXaxis()->SetTitleOffset(0.9);
    hZ[i]->GetXaxis()->SetTitleSize(0.05);
    hZ[i]->GetXaxis()->SetLabelSize(0.045);
    hZ[i]->GetYaxis()->SetTitleOffset(0.9);
    hZ[i]->GetYaxis()->SetTitleSize(0.05);
    hZ[i]->GetYaxis()->SetLabelSize(0.045);
    hZ[i]->GetXaxis()->SetTitle("#phi_{S}");
    hZ[i]->GetYaxis()->SetTitle("#phi_{Q}");
  }

  Double_t phi[2];

  for(Int_t bx = 1 ; bx <= nphiS ; bx++){

    for(Int_t by = 1 ; by <= nphiQ ; by++){

      phi[0] = hZ[0]->GetXaxis()->GetBinCenter(bx);
      phi[1] = hZ[0]->GetYaxis()->GetBinCenter(by);

      hZ[0]->SetBinContent(hZ[0]->GetBin(bx,by),Partition(phi,par));
      hZ[1]->SetBinContent(hZ[1]->GetBin(bx,by),CorrPiplus(phi,par));
      hZ[2]->SetBinContent(hZ[2]->GetBin(bx,by),CorrPiminus(phi,par));
      hZ[3]->SetBinContent(hZ[3]->GetBin(bx,by),CorrKminus(phi,par));
      hZ[4]->SetBinContent(hZ[4]->GetBin(bx,by),CorrKplus(phi,par));
      hZ[5]->SetBinContent(hZ[5]->GetBin(bx,by),CorrKzero(phi,par));
      hZ[6]->SetBinContent(hZ[6]->GetBin(bx,by),CorrAKzero(phi,par));
      hZ[7]->SetBinContent(hZ[7]->GetBin(bx,by),CorrProton(phi,par));
      hZ[8]->SetBinContent(hZ[8]->GetBin(bx,by),CorrAProton(phi,par));
      hZ[9]->SetBinContent(hZ[9]->GetBin(bx,by),CorrNeutron(phi,par));
      hZ[10]->SetBinContent(hZ[10]->GetBin(bx,by),CorrANeutron(phi,par));
      hZ[11]->SetBinContent(hZ[11]->GetBin(bx,by),CorrLa(phi,par));
      hZ[12]->SetBinContent(hZ[12]->GetBin(bx,by),CorrALa(phi,par));
      hZ[13]->SetBinContent(hZ[13]->GetBin(bx,by),CorrSigmaplus(phi,par));
      hZ[14]->SetBinContent(hZ[14]->GetBin(bx,by),CorrASigmaplus(phi,par));
      hZ[15]->SetBinContent(hZ[15]->GetBin(bx,by),CorrSigmaminus(phi,par));
      hZ[16]->SetBinContent(hZ[16]->GetBin(bx,by),CorrASigmaminus(phi,par));
      hZ[17]->SetBinContent(hZ[17]->GetBin(bx,by),CorrDeltaminus(phi,par));
      hZ[18]->SetBinContent(hZ[18]->GetBin(bx,by),CorrADeltaminus(phi,par));
      hZ[19]->SetBinContent(hZ[19]->GetBin(bx,by),CorrDeltaplusplus(phi,par));
      hZ[20]->SetBinContent(hZ[20]->GetBin(bx,by),CorrADeltaplusplus(phi,par));
      hZ[21]->SetBinContent(hZ[21]->GetBin(bx,by),CorrKsiminus(phi,par));
      hZ[22]->SetBinContent(hZ[22]->GetBin(bx,by),CorrAKsiminus(phi,par));
      hZ[23]->SetBinContent(hZ[23]->GetBin(bx,by),CorrKsi0(phi,par));
      hZ[24]->SetBinContent(hZ[24]->GetBin(bx,by),CorrAKsi0(phi,par));
      hZ[25]->SetBinContent(hZ[25]->GetBin(bx,by),CorrOmega(phi,par));
      hZ[26]->SetBinContent(hZ[26]->GetBin(bx,by),CorrAOmega(phi,par));

    }

  }

}

//__________________________________________________________________________
Int_t TTMThermalModelCanBSQ::GenerateParticleDens()
{
  // Calculates the Primordial particle densities and populates the density 
  // hash table. The Wroblewski factor and the decay contributions are 
  // also calculated provided the decays have been entered into the 
  // particle set through TTMParticleSet::InputDecays(). 
  //
   
  Int_t check = 1;

  check = PrimPartDens();
  
  if(!check){ 
    CalcWroblewski();
    GenerateDecayPartDens();
  }

  return check;
}

//__________________________________________________________________________
void TTMThermalModelCanBSQ::GenerateEnergyDens()
{
  // Iterates through the density hash table calculating the primordial
  // energy density of each particle in the hash table. Must first run 
  // GenerateParticleDens() to populate the hash table. If the parameters 
  // change then GenerateParticleDens() should be run again before this
  // function to calculate the canonical correction factors.
  //  

  TIter next(fDensTable);
  TTMDensObj *dens;

  fEnergy = 0.;

  while ((dens = (TTMDensObj *) next())) {
      
    TTMParticle *part = fPartSet->GetParticle(dens->GetID());
    Double_t CorrFactor = GetCorrFactor(part);

    TTMThermalParticleCanBSQ *ptr = new TTMThermalParticleCanBSQ(part, 
                                                                 fParm, CorrFactor);
    Double_t PartEnergy;
      
    if(part->GetWidth() == 0 || !fWidth){
      PartEnergy = ptr->EnergyBoltzmannNoWidth();
    }else{
      PartEnergy = ptr->EnergyBoltzmannWidth();
    }	      
    delete ptr;      
      
    dens->SetPrimEnergy(PartEnergy);
    fEnergy += PartEnergy;
  }
}

//__________________________________________________________________________
void TTMThermalModelCanBSQ::GenerateEntropyDens()
{
  // Iterates through the density hash table calculating the primordial
  // entropy density of each particle in the hash table. Must first run 
  // GenerateParticleDens() to populate the hash table. If the parameters 
  // change then GenerateParticleDens() should be run again before this
  // function to calculate the canonical correction factors. In the 
  // canonical approach the total entropy can't be split entirely into
  // the sum of particle entropies (i.e. there is a contribution to the 
  // total entropy not included in the hash table entries!).   
  //

  TIter next(fDensTable);
  TTMDensObj *dens;

  fEntropy = 0.;

  while ((dens = (TTMDensObj *) next())) {
      
    TTMParticle *part = fPartSet->GetParticle(dens->GetID());
    Double_t CorrFactor = GetCorrFactor(part);

    TTMThermalParticleCanBSQ *ptr = new TTMThermalParticleCanBSQ(part, 
                                                                 fParm, CorrFactor);
    Double_t PartEntropy;

    if(part->GetWidth() == 0 || !fWidth){
      PartEntropy = ptr->EnergyBoltzmannNoWidth()/fParm->GetT();
    }else{
      PartEntropy = ptr->EnergyBoltzmannWidth()/fParm->GetT();
    }	      
    delete ptr;      

    dens->SetPrimEntropy(PartEntropy);
    fEntropy += PartEntropy;
  }
   
  fEntropy += flnZtot/fParm->GetVolume();
}

//__________________________________________________________________________
void TTMThermalModelCanBSQ::GeneratePressure()
{
  // Iterates through the density hash table calculating the pressure
  // of each particle in the hash table. Must first run 
  // GenerateParticleDens() to populate the hash table. If the parameters 
  // change then GenerateParticleDens() should be run again before this
  // function to calculate the canonical correction factors.
  //  

  TIter next(fDensTable);
  TTMDensObj *dens;

  fPressure = 0.;

  while ((dens = (TTMDensObj *) next())) {
      
    TTMParticle *part = fPartSet->GetParticle(dens->GetID());
    Double_t CorrFactor = GetCorrFactor(part);

    TTMThermalParticleCanBSQ *ptr = new TTMThermalParticleCanBSQ(part, 
                                                                 fParm, CorrFactor);
    Double_t PartPressure;
      
    if(part->GetWidth() == 0 || !fWidth){
      PartPressure = ptr->PressureBoltzmannNoWidth();
    }else{
      PartPressure = ptr->PressureBoltzmannWidth();
    }	      
    delete ptr;      
      
    dens->SetPrimPressure(PartPressure);
    fPressure += PartPressure;
  }
}

//__________________________________________________________________________
void TTMThermalModelCanBSQ::ListInfo()
{
  // List model information
  //

  cout << "  ************************************************************"
       << "*****************" << endl;
  cout << "  ***************************** Thermal Model Info ***********"
       << "*****************" << endl << endl;

  cout << "\t Particle set: " << endl << "\t\t" << fPartSet->GetFilename()
       << endl << endl;

  cout << "\t Boltzmann Statistics " << endl;
   
  if (fWidth) {
    cout << "\t Resonance width included " << endl;
  } else {
    cout << "\t Resonance width excluded " << endl;
  }
  cout << endl << endl;

  fParm->List();
  cout << "  ***************************** Thermal Quantities ***********"
       << "****************** " << endl << endl;

  cout << "\t S/V        = ";
  cout.width(10);
  cout << fStrange;
  cout << endl << endl;
  cout << "\t B/2Q     = ";
  cout.width(10);
  cout << fBaryon / 2. / fCharge;
  cout << "\t";
  cout << endl << endl;
  cout << "\t lambda_s = ";
  cout.width(8);
  cout << fWroblewski << endl;
  cout << endl << endl;
  cout << "\t Primordial Densities: " << endl;
  if (fDensity != 0.) {
    cout << "\t\t\t\t n = " << fDensity << endl;
  }
  if (fEnergy != 0.) {
    cout << "\t\t\t\t e = " << fEnergy << endl;
  }
  if (fEntropy != 0.) {
    cout << "\t\t\t\t s = " << fEntropy << endl << endl;
  }
  if (fPressure != 0.) {
    cout << "\t\t\t\t P = " << fPressure << endl << endl;
  }
  ListStableParticles();
  cout << endl << endl;
  cout << "  ***********************************************************"
       << "*******************" << endl;
  cout << "  ***********************************************************"
       << "*******************" << endl << endl;
}

//__________________________________________________________________________
Double_t TTMThermalModelCanBSQ::GetCorrFactor(TTMParticle * part) const
{
  // Returns the canonical correction factor. Must first run 
  // GeneratePartDens() to calculate these!
  //

  Double_t B = part->GetB();
  Double_t S = part->GetS();
  Double_t Q = part->GetQ();

  Double_t CorrFactor = 1.;

  if (S == 0 && B == 0 && Q == 0) {
    CorrFactor = 1.;
  } else if (S == 0 && B == 0 && Q == +1) {
    CorrFactor = fCorrpip;
  } else if (S == 0 && B == 0 && Q == -1) {
    CorrFactor = fCorrpim;
  } else if (S == -1 && B == 0 && Q == -1) {
    CorrFactor = fCorrkm;
  } else if (S == +1 && B == 0 && Q == +1) {
    CorrFactor = fCorrkp;
  } else if (S == 1 && B == 0 && Q == 0) {
    CorrFactor = fCorrk0;
  } else if (S == -1 && B == 0 && Q == 0) {
    CorrFactor = fCorrak0;
  } else if (S == 0 && B == 1 && Q == +1) {
    CorrFactor = fCorrproton;
  } else if (S == 0 && B == -1 && Q == -1) {
    CorrFactor = fCorraproton;
  } else if (S == 0 && B == 1 && Q == 0) {
    CorrFactor = fCorrneutron;
  } else if (S == 0 && B == -1 && Q == 0) {
    CorrFactor = fCorraneutron;
  } else if (S == -1 && B == 1 && Q == 0) {
    CorrFactor = fCorrlambda;
  } else if (S == 1 && B == -1 && Q == 0) {
    CorrFactor = fCorralambda;
  } else if (S == -1 && B == 1 && Q == +1) {
    CorrFactor = fCorrsigmap;
  } else if (S == 1 && B == -1 && Q == -1) {
    CorrFactor = fCorrasigmap;
  } else if (S == -1 && B == 1 && Q == -1) {
    CorrFactor = fCorrsigmam;
  } else if (S == 1 && B == -1 && Q == +1) {
    CorrFactor = fCorrasigmam;
  } else if (S == 0 && B == 1 && Q == -1) {
    CorrFactor = fCorrdeltam;
  } else if (S == 0 && B == -1 && Q == +1) {
    CorrFactor = fCorradeltam;
  } else if (S == 0 && B == 1 && Q == +2) {
    CorrFactor = fCorrdeltapp;
  } else if (S == 0 && B == -1 && Q == -2) {
    CorrFactor = fCorradeltapp;
  } else if (S == -2 && B == 1 && Q == -1) {
    CorrFactor = fCorrksim;
  } else if (S == 2 && B == -1 && Q == +1) {
    CorrFactor = fCorraksim;
  } else if (S == -2 && B == 1 && Q == 0) {
    CorrFactor = fCorrksi0;
  } else if (S == 2 && B == -1 && Q == 0) {
    CorrFactor = fCorraksi0;
  } else if (S == -3 && B == 1 && Q == -1) {
    CorrFactor = fCorromega;
  } else if (S == 3 && B == -1 && Q == +1) {
    CorrFactor = fCorraomega;
  }

  return CorrFactor;
}
