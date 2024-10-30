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

// Author: Spencer Wheaton 22 September 2009    //
// Adapted for GSL: Yves Schutz October 2017    //

#include "FncsConstrain.h"

#include "TThermalModelBSQ.h"

//__________________________________________________________________________
Int_t BSQfuncBSQ(const gsl_vector* x, void* p, gsl_vector* f)
{
    Int_t rv = 0;
    TTMThermalModelBSQ* model = ((PARAMETERSS *)p)->p0;
    (model->GetParameterSet())->GetParameter(1)->SetValue(gsl_vector_get(x, 0));
    (model->GetParameterSet())->GetParameter(2)->SetValue(gsl_vector_get(x, 1));
    (model->GetParameterSet())->GetParameter(3)->SetValue(gsl_vector_get(x, 2));
    
    Int_t check = model->PrimPartDens();
    
    if (!check) {
        
        Double_t vol = model->GetParameterSet()->GetVolume();
        Double_t y0 = ((PARAMETERSS *)p)->p1;
        Double_t y1 = ((PARAMETERSS *)p)->p2;
        Double_t y2 = ((PARAMETERSS *)p)->p3;
        
        Double_t scale = 1.;
        
        if (y0 != 0.)
            gsl_vector_set(f, 0, (model->GetBaryon() * vol - y0) / y0 * scale);
        else
            gsl_vector_set(f, 0, (model->GetBaryon() * vol - y0) / ((TMath::Abs(model->GetBplus()) + TMath::Abs(model->GetBminus())) * vol) * scale);
        
        if (y1 != 0.)
            gsl_vector_set(f, 1, (model->GetStrange() * vol - y1) / y1 * scale);
        else
            gsl_vector_set(f, 1, (model->GetStrange() * vol - y1) / ((TMath::Abs(model->GetSplus()) + TMath::Abs(model->GetSminus())) * vol) * scale);
        
        if (y2 != 0.)
            gsl_vector_set(f, 2, (model->GetCharge() * vol - y2) / y2 * scale);
        else
            gsl_vector_set(f, 2, (model->GetCharge() * vol - y2) / ((TMath::Abs(model->GetQplus()) + TMath::Abs(model->GetQminus())) * vol) * scale);
    } else {
        cout << "Primary particles density problems!" << endl;
        rv = 1;
    }
    return rv;
}
