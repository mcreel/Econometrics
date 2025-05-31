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
  #include "maxsa.h"
#else
  #include <oxfloat.h>

/*========================= function maximization ==========================*/

static decl s_mxEval = 1e5;
static decl s_dEps = 1e-6, s_iNEps= 4;        /* convergence criteria */
static decl s_iPrint = 0;                     /* print results */
static decl s_iNS= 20, s_iNT= M_NAN, s_vC = 2, s_vM= 1, s_dRT= 0.85;
static decl s_bInit= FALSE;

MaxSAConvergenceMsg(const iCode)
{
    if (iCode == MAXSA_CONV)
        return "Strong convergence";
    else if (iCode == MAXSA_MAXEV)
        return "No convergence (maximum no of function evaluations reached)";
    else if (iCode == MAXSA_FUNC_FAIL)
        return "No convergence (function evaluation failed)";
    else if (iCode == MAXSA_TEMP)
        return "No convergence (initial temperature negative)";
    else
        return "No convergence (for unknown reason)";
}

/*
**  MaxSAControl(const mxEval, const iPrint)
**   mxEval - The maximum number of function evaluations. If it is
**            exceeded, IER = 1. (INT)
**   iPrint - controls printing inside SA. (INT)
**            Values: 0 - Nothing printed.
**                    1 - Function value for the starting value and
**                        summary results before each temperature
**                        reduction. This includes the optimal
**                        function value found so far, the total
**                        number of moves (broken up into uphill,
**                        downhill, accepted and rejected), the
**                        number of out of bounds trials, the
**                        number of new optima found at this
**                        temperature, the current optimal X and
**                        the step length VM. Note that there are
**                        N*NS*NT function evalutations before each
**                        temperature reduction. Finally, notice is
**                        is also given upon achieveing the termination
**                        criteria.
**                    2 - Each new step length (VM), the current optimal
**                        X (XOPT) and the current trial X (X). This
**                        gives the user some idea about how far X
**                        strays from XOPT as well as how VM is adapting
**                        to the function.
**            Suggested value: 1
**            Note: For a given value of IPRINT, the lower valued
**                  options (other than 0) are utilized.
*/
MaxSAControl(const mxEval, const iPrint)
{
    if (mxEval >= 0) ::s_mxEval = mxEval;
    if (iPrint >= 0) ::s_iPrint = iPrint;
}
/*
**  MaxSAControlEps(const dEps, const iNEps)
**
**  Purpose:
**    Control the precision of the optimization
**
**   EPS -  Error tolerance for termination. If the final function
**          values from the last neps temperatures differ from the
**          corresponding value at the current temperature by less than
**          EPS and the final function value at the current temperature
**          differs from the current optimal function value by less than
**          EPS, execution terminates and IER = 0 is returned. (EP)
**   NEPS - Number of final function values used to decide upon termi-
**          nation. See EPS. Suggested value is 4. (INT)
*/
MaxSAControlEps(const dEps, const iNEps)
{
    if (dEps > 0) ::s_dEps = dEps;
    if ((iNEps > 0) && isint(iNEps)) ::s_iNEps = iNEps;
}

/*
**  MaxSAControlStep(const iNS, const iNT, const vC)
**
**  Purpose:
**    Control the size of the steps
**
**   NS - Number of cycles. After NS*N function evaluations, each
**        element of VM is adjusted so that approximately half of
**        all function evaluations are accepted. The suggested value
**        is 20. (INT)
**   NT - Number of iterations before temperature reduction. After
**        NT*NS*N function evaluations, temperature (T) is changed
**        by the factor RT. Value suggested by Corana et al. is
**        MAX(100, 5*N). See Goffe et al. for further advice. (INT)
**        Value of M_NAN leads to computation of MAX(100, 5*N)
**   RT - The temperature reduction factor. The value suggested by
**        Corana et al. is .85. See Goffe et al. for more advice. (DP)
**   VM - The step length vector. On input it should encompass the
**        region of interest given the starting value X. For point
**        X(I), the next trial point is selected is from X(I) - VM(I)
**        to  X(I) + VM(I). Since VM is adjusted so that about half
**        of all points are accepted, the input value is not very
**        important (i.e. is the value is off, SA adjusts VM to the
**        correct value). (DP(N))
**   C - Vector that controls the step length adjustment. The suggested
**       value for all elements is 2.0. (DP(N))
*/
MaxSAControlStep(const iNS, const iNT, const dRT, const vM, const vC)
{
  if (iNS > 0) :: s_iNS= iNS;
  if ((iNT > 0) || isnan(iNT)) ::s_iNT= iNT;
  if (dRT > 0) ::s_dRT= dRT;
  if (vM > 0) ::s_vM= vM;
  if (vC > 0) ::s_vC= vC;
  
  if (s_dRT >= 1)
    println ("Warning: Temperature reduction factor >= 1, no temperature reduction");
  s_bInit= TRUE;    
}
GetMaxSAControl()
{
	return { s_mxEval, s_iPrint };
}
GetMaxSAControlEps()
{
	return { s_dEps, s_iNEps, s_vC };
}
GetMaxSAControlStep()
{
	return { s_iNS, s_iNT, s_dRT, s_vM, s_vC };
}


/*
**  MaxSA
**
** Input Parameters:
**
** Input/Output Parameters:
**   T - On input, the initial temperature. See Goffe et al. for advice.
**       On output, the final temperature. (DP)
*/
MaxSA(const func, const avP, const adFunc, const adT)
{
  decl ir, inAcc, inObds, inFnEv, j, h, m, dFOpt, vXOpt, vFStar, 
       vnAcp, nUp, nRej, nNew, nDown, vX, vXP, dF, dFP, dU, dP, iN,
       vRatio, bRep, iNT, vM;

  inAcc= inObds= inFnEv= 0;
  vX= vXOpt= avP[0];
  iN= sizerc(avP[0]);
  
  iNT= isnan(s_iNT) ? max(100, 5*iN) : s_iNT;
  vM= ones(iN, 1) .* s_vM;
  vFStar= constant(M_INF_NEG, s_iNEps, 1);
  
  ir= func(vX, &dF, 0, 0);
  ++inFnEv;
  
  if (adT[0] <= 0)
    return MAXSA_TEMP;
  if (!ir)
    return MAXSA_FUNC_FAIL;

  if (!s_bInit)
    println ("Warning: Use MaxSAControlStep to init step-parameters;\n",
             "         default values may be suboptimal.");  

  if (s_iPrint > 0)
    print("Initial result ", double(dF), "%r", {"at parameters"}, vX');

  dFOpt= vFStar[0]= dF;
  bRep= TRUE;

  while (bRep)
    {  
      nUp= nRej= nNew= nDown= 0;
      for (m= 0; m < iNT; ++m)
        {
          vnAcp= zeros(vX);
          for (j= 0; j < s_iNS; ++j)
            for (h= 0; h < iN; ++h)
              {
                vXP= vX;
                ir= 0;
                while (!ir) 
                  {
                    vXP[h]= vX[h] + (2*ranu(1,1)-1)*vM[h];

                    ir= func(vXP, &dFP, 0, 0);
                    ++inFnEv;
                    if (inFnEv > s_mxEval)
                      {
                        println ("Error: Too many function evaluations");
                        avP[0]= vXOpt;
                        adFunc[0]= dFOpt;
                        return MAXSA_MAXEV;
                      }  
                  }  

                if (dFP >= dF)
                  {
                    vX= vXP;
                    dF= dFP;
                    ++inAcc;
                    ++vnAcp[h];
                    ++nUp;

                    if (dFP > dFOpt)
                      {
                        vXOpt= vXP;
                        dFOpt= dFP;
                        ++nNew;
                      }      
                  }
                else
                  {
                    dP= exp((dFP-dF)/adT[0]);
                    dU= ranu(1,1);
                    if (dU < dP)
                      { // Accept
                        vX= vXP;
                        dF= dFP;
                        ++inAcc;
                        ++vnAcp[h];
                        ++nDown;
                      }
                    else
                      ++nRej;
                  }
              }
          vRatio= vnAcp / s_iNS;
          vM= vRatio .> .6 .? vM.*(1 + s_vC .* (vRatio-.6) ./ .4) .:
              vRatio .< .4 .? vM./(1 + s_vC .* (.4-vRatio) ./ .4) .:
              vM;
        
          if (s_iPrint > 1)
            print ("\nIntermediate results after step length adjustment",
                   "%c", {"vM", "vX", "vXOpt"}, vM~vX~vXOpt);
        }      

      vFStar[0]= dF;
      bRep= (dFOpt - vFStar[0] > s_dEps) || 
               !(fabs(dF - vFStar) <= s_dEps);
      
      vX= vXOpt;
      dF= dFOpt;
      vFStar= lag0(vFStar, 1, M_NAN);
  
      // Adapt temperature
      if (bRep)
        {
          if (s_iPrint > 0)
            print ("\nIntermediate results before next temperature reduction",
                   "%r", {"Current temperature", "Min func", "Total moves", 
                          "Downhill", "Accepted uphill", "Rejected uphill",
                          "New minima"},
                   adT[0]|dFOpt|(nUp+nDown+nRej)|nUp|nDown|nRej|nNew);
          adT[0]*= s_dRT;
        }  
    }
  
  if (s_iPrint > 0)
    println ("Function evaluations: ", double(inFnEv));

  avP[0]= vXOpt;
  adFunc[0]= dFOpt;
  
  return MAXSA_CONV;  
}
#endif
