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
// Adapted for GSL: Yves Schutz September 2017  //

#include "FncsConstrain.h"

#include "TTMParameterSetBQ.h"
#include "TThermalModelBQ.h"

//__________________________________________________________________________
Int_t BQConstrainBDens(TTMThermalModelBQ *model, Double_t nb)
{
    Int_t rv = 0;
    model->GetParameterSet()->GetParameter(1)->SetStatus("");
    if (nb == 0.) {
        cout << "Cannot constrain nb to zero" << endl;
        rv = 1;
    } else {
        const size_t ndim = 1;
        gsl_vector *x = gsl_vector_alloc(ndim);
        gsl_vector_set(x, 0, model->GetParameterSet()->GetMuB());
        Int_t  check = 0;
        PARAMETERS p;
        p.p0 = model;
        p.p1 = nb;
        broyden(x, ndim, check, p, BQfuncBDens);
        if(check) {
            cout << gsl_strerror(check) << endl;
            model->GetParameterSet()->SetConstraintInfo("Unable to Constrain nb");
            model->GetParameterSet()->GetParameter(1)->SetStatus("(Unable to constrain)");
            model->GetParameterSet()->GetParameter(1)->SetValue(0.);
            rv = check;
        } else {
            model->GetParameterSet()->GetParameter(1)->SetValue(gsl_vector_get(x, 0));
            model->GetParameterSet()->SetConstraintInfo("nb Successfully Constrained");
            model->GetParameterSet()->GetParameter(1)->SetStatus("(*CONSTRAINED*)");
        }
        gsl_vector_free(x);
    }
    return rv;
}
