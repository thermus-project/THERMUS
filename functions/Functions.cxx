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

#include <TROOT.h>
#include <TObject.h>
#include <TMath.h>

using namespace std;

Double_t FcnDistr(Double_t *x, Double_t *par)
{
  // x[0] : p/T

  Double_t mu = par[0];
  Double_t mass = par[1];
  Double_t G = par[2];
  Double_t Stat = par[3];

  Double_t E = TMath::Sqrt(x[0]*x[0] + pow(mass,2));

  return  G*exp(-E + mu)/(1 + Stat*G*exp(-E + mu));
}


Double_t FcnDens(Double_t *x, Double_t *par)
{
  // x[0] : p/T
  // This function returns 1/T^3*dn/d(p/T)- it is calculated based on 
  // n_i = -d(Omega_GC)/d(mu_i)
  //

  Double_t mu = par[0];
  Double_t mass = par[1];
  Double_t G = par[2];
  Double_t Stat = par[3];
  Double_t Deg = par[4];

  Double_t E = TMath::Sqrt(x[0]*x[0] + pow(mass,2));

  return 1. / (2. * pow(TMath::Pi(), 2))*Deg*x[0]*x[0]*G*exp(-E + mu)/(1 + Stat*G*exp(-E + mu));
}

Double_t FcnEnergyDens(Double_t *x, Double_t *par)
{
  // x[0] : p/T
  // This function returns 1/T^4*de/d(p/T)- it follows from the definition 
  // of average energy
  //

  Double_t mu = par[0];
  Double_t mass = par[1];
  Double_t G = par[2];
  Double_t Stat = par[3];
  Double_t Deg = par[4];

  Double_t E = TMath::Sqrt(x[0]*x[0] + pow(mass,2));

  return 1. / (2. * pow(TMath::Pi(), 2))*Deg*x[0]*x[0]*G*E*exp(-E + mu)/(1 + Stat*G*exp(-E + mu));
}

Double_t FcnEntropyDens(Double_t *x, Double_t *par)
{
  // x[0] : p/T
  // This function returns 1/T^3*ds/d(p/T)- it is calculated based on
  // s = -d(Omega_GC)/dT
  //
 
  Double_t mu = par[0];
  Double_t mass = par[1];
  Double_t G = par[2];
  Double_t Stat = par[3];
  Double_t Deg = par[4];

  Double_t E = TMath::Sqrt(x[0]*x[0] + pow(mass,2));

  return 1. / (2. * pow(TMath::Pi(), 2))*Deg*x[0]*x[0]*G*exp(-E + mu)/(1 + Stat*G*exp(-E + mu))*(x[0]*x[0]/3./E + E -mu);
}

Double_t FcnPressure(Double_t *x, Double_t *par)
{
  // x[0] : p/T
  // This function returns 1/T^4*dP/d(p/T)- it is calculated based on
  // P = -d(Omega_GC)/dV
  //

  Double_t mu = par[0];
  Double_t mass = par[1];
  Double_t G = par[2];
  Double_t Stat = par[3];
  Double_t Deg = par[4];

  Double_t E = TMath::Sqrt(x[0]*x[0] + pow(mass,2));

  return 1. / (2. * pow(TMath::Pi(), 2))*Deg*x[0]*x[0]*G*exp(-E + mu)/(1 + Stat*G*exp(-E + mu))*x[0]*x[0]/3./E;
}

Double_t FcnDensWidth(Double_t *x, Double_t *par)
{
  // x[0] : p/T
  // x[1] : m/T
  // This function returns 1/T^3 * d^2n/d(p/T)d(m/T) un-normalised
  // consistent with the no-width density
  //

  Double_t mu = par[0];		// mu/T
  Double_t mass = par[1];	// mass/T
  Double_t G = par[2];
  Double_t Stat = par[3];
  Double_t Deg = par[4];
  Double_t width = par[5];	// width/T

  Double_t E = TMath::Sqrt(x[0]*x[0] + x[1]*x[1]);

  return 1. / (2. * pow(TMath::Pi(), 2)) * Deg * x[0] * x[0] * G * exp(-E+mu) / (1. + G * exp(-E+mu) * Stat) / TMath::Pi() * mass * width * 2. * x[1] / ( pow((x[1] * x[1] - mass * mass), 2) + mass * mass * width * width);

}

Double_t FcnEnergyDensWidth(Double_t *x, Double_t *par)
{
  // x[0] : p/T
  // x[1] : m/T
  // This function returns 1/T^4 * d^2e/d(p/T)d(m/T) un-normalised
  // consistent with the no-width energy density
  //

  Double_t mu = par[0];		// mu/T
  Double_t mass = par[1];	// mass/T
  Double_t G = par[2];
  Double_t Stat = par[3];
  Double_t Deg = par[4];
  Double_t width = par[5];	// width/T

  Double_t E = TMath::Sqrt(x[0]*x[0] + x[1]*x[1]);

  return 1. / (2. * pow(TMath::Pi(), 2)) * Deg * x[0] * x[0] * E * G * exp(-E+mu) / (1. + G * exp(-E+mu) * Stat) / TMath::Pi() * mass * width * 2. * x[1] / ( pow((x[1] * x[1] - mass * mass), 2) + mass * mass * width * width);

}

Double_t FcnEntropyDensWidth(Double_t *x, Double_t *par)
{
  // x[0] : p/T
  // x[1] : m/T
  // This function returns 1/T^3 * d^2s/d(p/T)d(m/T) un-normalised
  // consistent with the no-width entropy density
  //

  Double_t mu = par[0];		// mu/T
  Double_t mass = par[1];	// mass/T
  Double_t G = par[2];
  Double_t Stat = par[3];
  Double_t Deg = par[4];
  Double_t width = par[5];	// width/T

  Double_t E = TMath::Sqrt(x[0]*x[0] + x[1]*x[1]);

  return 1. / (2. * pow(TMath::Pi(), 2)) * Deg * x[0] * x[0] * G * exp(-E+mu) / (1. + G * exp(-E+mu) * Stat) * (x[0] * x[0] / 3. / E + E - mu) / TMath::Pi() * mass * width * 2. * x[1] / ( pow((x[1] * x[1] - mass * mass), 2) + mass * mass * width * width);

}

Double_t FcnPressureWidth(Double_t *x, Double_t *par)
{
  // x[0] : p/T
  // x[1] : m/T
  // This function returns 1/T^4 * d^2P/d(p/T)d(m/T) un-normalised
  // consistent with the no-width pressure
  //

  Double_t mu = par[0];		// mu/T
  Double_t mass = par[1];	// mass/T
  Double_t G = par[2];
  Double_t Stat = par[3];
  Double_t Deg = par[4];
  Double_t width = par[5];	// width/T

  Double_t E = TMath::Sqrt(x[0]*x[0] + x[1]*x[1]);

  return 1. / (2. * pow(TMath::Pi(), 2)) * Deg * x[0] * x[0] * G * exp(-E+mu) / (1. + G * exp(-E+mu) * Stat) * x[0] * x[0] / 3. / E  / TMath::Pi() * mass * width * 2. * x[1] / ( pow((x[1] * x[1] - mass * mass), 2) + mass * mass * width * width);

}

Double_t FcnDensNormWidth(Double_t *x, Double_t *par)
{
  // x[0] : m/T
  // This function when integrated gives the normalisation of the 
  // Breit-Wigner distribution
  //

  Double_t mass = par[0];	// mass/T
  Double_t width = par[1];	// width/T

  return 1. / TMath::Pi() * mass * width * 2. * x[0] / ( pow((x[0] * x[0] - mass * mass), 2) + mass * mass * width * width);

}

Double_t FcnDensBoltzmannWidth(Double_t *x, Double_t *par)
{
  // x[0] : m/T
  // This function returns 1/T^3 * dn/d(m/T)

  Double_t mu = par[0];		// mu/T
  Double_t mass = par[1];	// mass/T
  Double_t G = par[2];
  Double_t Deg = par[3];
  Double_t width = par[4];	// width/T

  return 1. / (2. * pow(TMath::Pi(), 2)) * G * Deg * exp(mu) * x[0] * x[0] * TMath::BesselK( 2, x[0]) / 
    TMath::Pi() * mass * width * 2. * x[0] / ( pow((x[0] * x[0] - mass * mass), 2) + mass * mass * width * width);

}

Double_t FcnEnergyBoltzmannWidth(Double_t *x, Double_t *par)
{
  // x[0] : m/T
  // This function returns 1/T^4 * de/d(m/T)

  Double_t mu = par[0];		// mu/T
  Double_t mass = par[1];	// mass/T
  Double_t G = par[2];
  Double_t Deg = par[3];
  Double_t width = par[4];	// width/T

  return 1. / (2. * pow(TMath::Pi(), 2)) * G * Deg * exp(mu) * x[0] * x[0] * x[0] * (TMath::BesselK( 1, x[0]) + 3. / x[0] * TMath::BesselK( 2, x[0])) / TMath::Pi() * mass * width * 2. * x[0] / ( pow((x[0] * x[0] - mass * mass), 2) + mass * mass * width * width);

}

Double_t FcnEntropyBoltzmannWidth(Double_t *x, Double_t *par)
{
  // x[0] : m/T
  // This function returns 1/T^3 * ds/d(m/T)

  Double_t mu = par[0];		// mu/T
  Double_t mass = par[1];	// mass/T
  Double_t G = par[2];
  Double_t Deg = par[3];
  Double_t width = par[4];	// width/T

  return 1. / (2. * pow(TMath::Pi(), 2)) * G * Deg * exp(mu) * x[0] * x[0] * (x[0] * TMath::BesselK( 1, x[0]) + (4. - mu) * TMath::BesselK( 2, x[0])) / TMath::Pi() * mass * width * 2. * x[0] / ( pow((x[0] * x[0] - mass * mass), 2) + mass * mass * width * width);

}

Double_t FcnEnergyFlucBoltzmannWidth(Double_t *x, Double_t *par)
{
  // x[0] : m/T
  // This function returns 1/T^4 * integrand

  Double_t mu = par[0];		// mu/T
  Double_t mass = par[1];	// mass/T
  Double_t G = par[2];
  Double_t Deg = par[3];
  Double_t width = par[4];	// width/T

  return 1. / (2. * pow(TMath::Pi(), 2)) * G * Deg * exp(mu) * x[0] * x[0] * x[0] * 
      ( TMath::BesselK(0,x[0]) + 5. / x[0] * TMath::BesselK(1,x[0]) + 12. * pow(1./x[0],2.) * TMath::BesselK(2,x[0]) ) 
      / TMath::Pi() * mass * width * 2. * x[0] / ( pow((x[0] * x[0] - mass * mass), 2) + mass * mass * width * width); 

}

Double_t GenericQStatDeriv(Double_t *x, Double_t *par)
{
    // this function returns
    // 1/T^3 * g/2/Pi^2 (e^{-E/T+mu/T})^a / (1+-e^{-E/T+mu/T})^b * E^c 

    // x[0] : p/T

    Double_t mu = par[0];   // mu/T
    Double_t mass = par[1]; // mass/T
    Double_t G = par[2];    // gamma
    Double_t Stat = par[3];
    Double_t Deg = par[4];
    Double_t a = par[5];
    Double_t b = par[6];
    Double_t c = par[7];

    Double_t E = TMath::Sqrt(x[0]*x[0] + pow(mass,2));

    return  1. / (2.* pow(TMath::Pi(),2))*pow(x[0],2)*Deg*pow(G,a)*pow(E,c)*
	pow(Stat,a+1.)*exp(-a*E + a*mu)/pow(1 + Stat*G*exp(-E + mu),b);

}

Double_t GenericQStatDerivWidth(Double_t *x, Double_t *par)
{
    // this function returns
    // 1/T^3 * g/2/Pi^2 (e^{-E/T+mu/T})^a / (1+-e^{-E/T+mu/T})^b * E^c * BW Factor

    // x[0] : p/T
    // x[1] : m/T

    Double_t mu = par[0];   // mu/T
    Double_t mass = par[1]; // mass/T
    Double_t G = par[2];    // gamma
    Double_t Stat = par[3];
    Double_t Deg = par[4];
    Double_t a = par[5];
    Double_t b = par[6];
    Double_t c = par[7];
    Double_t width = par[8]; // in fact width/T

    Double_t E = TMath::Sqrt(x[0]*x[0] + x[1]*x[1]);  // in fact E/T

    return  1. / (2.* pow(TMath::Pi(),2))*pow(x[0],2)*Deg*pow(G,a)*pow(E,c)*
	pow(Stat,a+1.)*exp(-a*E + a*mu)/pow(1 + Stat*G*exp(-E + mu),b)/ 
	TMath::Pi()*mass*width*2.*x[1]/(pow((x[1]*x[1]-mass*mass), 2)+mass*mass*width*width);

}
