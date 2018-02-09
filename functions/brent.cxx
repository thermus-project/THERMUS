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

// Author: Yves Schutz 22 Janvier 2018          //

#include <gsl/gsl_roots.h>

#include "FncsConstrain.h"

//__________________________________________________________________________
Double_t brent(Double_t xLow, Double_t xHigh, Double_t step, Int_t& status, PARAMETERSS p, Double_t (*myfunction)(Double_t x, void* p))
{
    Int_t iter = 0, max_iter = 100;
    const gsl_root_fsolver_type *T;
    gsl_root_fsolver *s;
    Double_t result = 0;
    
    gsl_function F;
    F.function = myfunction;
    F.params = &p;
    
    T = gsl_root_fsolver_brent;
    s = gsl_root_fsolver_alloc (T);
    gsl_root_fsolver_set (s, &F, xLow, xHigh);
    do
    {
        iter++;
        status = gsl_root_fsolver_iterate (s);
        result = gsl_root_fsolver_root (s);
        xLow = gsl_root_fsolver_x_lower (s);
        xHigh = gsl_root_fsolver_x_upper (s);
        status = gsl_root_test_interval (xLow, xHigh, 0, step);
    }
    while (status == GSL_CONTINUE && iter < max_iter);
    
    gsl_root_fsolver_free (s);
    
    return result;
}
