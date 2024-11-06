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

#include "TThermalModelBSQ.h"

//__________________________________________________________________________
Int_t BSQfuncSQST3(const gsl_vector* x, void* p, gsl_vector* f)
{
    Int_t rv = 0;
    TTMThermalModelBSQ* model = ((PARAMETERSS *)p)->p0;
    
    (model->GetParameterSet())->GetParameter(2)->SetValue(gsl_vector_get(x, 0));
    (model->GetParameterSet())->GetParameter(3)->SetValue(gsl_vector_get(x, 1));
    (model->GetParameterSet())->GetParameter(1)->SetValue(gsl_vector_get(x, 2));
    
    Int_t check = model->PrimPartDens();
    
    if (!check) {
        model->GenerateEntropyDens();
        
        Double_t ss0 = ((PARAMETERSS *)p)->p1;
        
        if (ss0 != 0.)
            gsl_vector_set(f, 0, (model->GetStrange() - ss0) / ss0);
        else
            gsl_vector_set(f, 0, (model->GetStrange() - ss0 )/ (TMath::Abs(model->GetSplus()) + TMath::Abs(model->GetSminus())));
        
        Double_t ss1 = ((PARAMETERSS *)p)->p2;
        gsl_vector_set(f, 1, (model->GetBaryon()/2./ model->GetCharge() - ss1) / ss1);
        
        Double_t ss2 = ((PARAMETERSS *)p)->p3;
        gsl_vector_set(f, 2, (model->GetEntropy()/ TMath::Power(model->GetParameterSet()->GetT(), 3.) * TMath::Power(0.197, 3.) - ss2) / ss2);
        
    } else {
        cout << "Primary particles density problems!" << endl;
        rv = 1;
    }
    return rv;
}
