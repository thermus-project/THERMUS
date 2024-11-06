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

#include <gsl/gsl_vector.h>

#include <TThermalModelBSQ.h>
#include <TThermalModelBQ.h>

struct PARAMETERS  { TTMThermalModelBQ* p0; Double_t p1;};
struct PARAMETERSS { TTMThermalModelBSQ* p0; Double_t p1; Double_t p2; Double_t p3;};

//__________________________________________________________________________
Int_t BSQConstrainSQ            (TTMThermalModelBSQ *model);
Int_t BSQConstrainSQC           (TTMThermalModelBSQ *model);
Int_t BSQConstrainSQCb          (TTMThermalModelBSQ *model);
Int_t BSQConstrainBSQ           (TTMThermalModelBSQ *model, Double_t B, Double_t S, Double_t Q);
Int_t BSQConstrainQQ            (TTMThermalModelBSQ *model, Double_t Q);
Int_t BSQConstrainS             (TTMThermalModelBSQ *model);
Int_t BSQConstrainQ             (TTMThermalModelBSQ *model);
Int_t BSQConstrainSQEN          (TTMThermalModelBSQ *model, Double_t eovern);
Int_t BSQConstrainSQPercolation (TTMThermalModelBSQ *model);
Int_t BSQConstrainSPercolation  (TTMThermalModelBSQ *model);
Int_t BSQConstrainSEN           (TTMThermalModelBSQ *model, Double_t eovern);
Int_t BSQConstrainSQST3         (TTMThermalModelBSQ *model, Double_t sovert3);
Int_t BSQConstrainSST3          (TTMThermalModelBSQ *model, Double_t sovert3);
Int_t BQConstrainQ              (TTMThermalModelBQ *model);
Int_t BQConstrainQEN            (TTMThermalModelBQ *model, Double_t eovern);
Int_t BQConstrainQPercolation   (TTMThermalModelBQ *model);
Int_t BQConstrainST3            (TTMThermalModelBQ *model, Double_t SoverT3);
Int_t BQConstrainBDens          (TTMThermalModelBQ *model, Double_t nb);
Int_t BQConstrainQNetBDens      (TTMThermalModelBQ *model, Double_t nb);
Int_t BQConstrainQST3           (TTMThermalModelBQ *model, Double_t SoverT3);
Int_t BQConstrainQBDens         (TTMThermalModelBQ *model, Double_t nb);
Int_t BQConstrainEN             (TTMThermalModelBQ *model, Double_t eovern);
Int_t BSQConstrainSBDens        (TTMThermalModelBSQ *model, Double_t nb);
Int_t BSQConstrainSQBDens       (TTMThermalModelBSQ *model, Double_t nb);
Int_t BSQConstrainSQNetBDens    (TTMThermalModelBSQ *model, Double_t nb);

Int_t BQfuncBDens       (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BQfuncEN          (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BQfuncQ           (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BQfuncQBDens      (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BQfuncQBDens      (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BQfuncQEN         (const gsl_vector* x, void *p, gsl_vector *f);
Int_t BQfuncQNetBDens   (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BQfuncQPercolation(const gsl_vector* x, void* p, gsl_vector* f);
Int_t BQfuncQST3        (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BQfuncST3         (const gsl_vector* x, void* p, gsl_vector* f);

Int_t BSQfuncBSQ           (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncQ             (const gsl_vector *x, void* p, gsl_vector* f);
Int_t BSQfuncQQ            (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncS             (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncSBDens        (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncSEN           (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncSPercolation  (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncSQ            (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncSQBDens       (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncSQC           (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncSQCb          (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncSQEN          (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncSQNetBDens    (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncSQPercolation (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncSQST3         (const gsl_vector* x, void* p, gsl_vector* f);
Int_t BSQfuncSST3          (const gsl_vector* x, void* p, gsl_vector* f);

void broyden(gsl_vector* x, size_t n, Int_t& status, PARAMETERS p, Int_t (*f)(const gsl_vector* x, void* p, gsl_vector* f));
void broyden(gsl_vector* x, size_t n, Int_t& status, PARAMETERSS p, Int_t (*f)(const gsl_vector* x, void* p, gsl_vector* f));


Double_t ExclVolPressureFunc(Double_t x, void * p);

Double_t brent(Double_t xLow, Double_t xHigh, Double_t step, Int_t& status, PARAMETERSS p, Double_t (*f)(Double_t x, void* p));


Int_t funcTest(const gsl_vector* x, void* p, gsl_vector* f);
