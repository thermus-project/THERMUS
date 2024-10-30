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
// Adapted for GSL: Yves Schutz October 2017    //

#include "FncsConstrain.h"

#include "TTMDensObj.h"
#include "TThermalModelBSQ.h"

//__________________________________________________________________________
Int_t  BSQfuncSQPercolation(const gsl_vector* x, void* p, gsl_vector* f)
{
    Int_t rv = 0;
    TTMThermalModelBSQ* model = ((PARAMETERSS *)p)->p0;
    
    (model->GetParameterSet())->GetParameter(2)->SetValue(gsl_vector_get(x, 0));
    (model->GetParameterSet())->GetParameter(3)->SetValue(gsl_vector_get(x, 1));
    (model->GetParameterSet())->GetParameter(1)->SetValue(gsl_vector_get(x, 2));
    
    Double_t vh = 4. / 3. *TMath::Pi() * TMath::Power(0.8, 3.);
    
    Int_t check = model->PrimPartDens();
    
    if (!check) {
        model->GenerateEnergyDens();
        Double_t p1 = ((PARAMETERSS *)p)->p1;
        if (p1 != 0.)
            gsl_vector_set(f, 0, (model->GetStrange() - p1) / p1);
        else
            gsl_vector_set(f, 0, (model->GetStrange() - p1) / (TMath::Abs(model->GetSplus()) + TMath::Abs(model->GetSminus())));
        Double_t p2 = ((PARAMETERSS *)p)->p2;
        gsl_vector_set(f, 1, (model->GetBaryon() / 2. / model->GetCharge() - p2) / p2);
        gsl_vector_set(f, 2, (model->GetDensity() - 1.24 / vh * (1. - model->GetBaryon() / model->GetDensity()) - 0.34 / vh * model->GetBaryon() / model->GetDensity()) / (model->GetDensity()));
    } else {
        cout << "Primary particles density problems!" << endl;
        rv = 1;
    }
    return rv;
}
