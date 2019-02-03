/*
**  MaxSA
**
**  Purpose:
**    Implement Simulated Annealing in Ox following 
**      Goffe, William L., Gary D. Ferrier, and John Rogers (1994).  
**      Global Optimization of Statistical Functions with Simulated Annealing.  
**      Journal of Econometrics, 60(1/2):65 99.
**      http://emlab.berkeley.edu/Software/abstracts/goffe895.html
**
**  Date:
**    2/10/2002
**
**  Author:
**    Charles Bos
**
*/
#ifndef MAXSA_INCLUDED
#define MAXSA_INCLUDED

enum
{
    MAXSA_CONV, MAXSA_MAXEV, MAXSA_FUNC_FAIL, MAXSA_TEMP
};
MaxSAConvergenceMsg(const iCode);
MaxSAControl(const mxEval, const iPrint);
MaxSAControlEps(const dEps, const iNEps);
MaxSAControlStep(const iNS, const iNT, const dRT, const vM, const vC);
GetMaxSAControl();
GetMaxSAControlEps();
GetMaxSAControlStep();
MaxSA(const func, const avP, const adFunc, const adT);

#include <packages/maxsa/maxsa.ox>

#endif /* MAXSA_INCLUDED */
