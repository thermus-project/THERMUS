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

#include <gsl/gsl_multiroots.h>

#include "FncsConstrain.h"

//__________________________________________________________________________
void broyden(gsl_vector* x, size_t ndim, Int_t& status, PARAMETERS p, Int_t (*myfunction)(const gsl_vector* x, void* p, gsl_vector* f))
{
    // the broyden function to solve mutidimentional equations
    
    gsl_multiroot_function function;
    function.f = myfunction;
    function.n = ndim;
    function.params = &p;
    
    const gsl_multiroot_fsolver_type *T;
    T = gsl_multiroot_fsolver_hybrid;
    
    
    
    gsl_multiroot_fsolver *s;
    s = gsl_multiroot_fsolver_alloc(T, ndim);
    gsl_multiroot_fsolver_set (s, &function, x);
    
    size_t iter = 0;
    do
    {
        iter++;
        status = gsl_multiroot_fsolver_iterate (s);
        if (status)   /* check if solver is stuck */
            break;
        status = gsl_multiroot_test_residual (s->f, 1e-7);
        for (size_t i = 0; i < ndim; i++)
            gsl_vector_set(x, i, gsl_vector_get(s->x, i));
    } while (status == GSL_CONTINUE && iter < 10);
    
    gsl_multiroot_fsolver_free (s);
}
