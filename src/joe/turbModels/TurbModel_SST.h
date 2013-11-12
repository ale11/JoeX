#ifndef RANSTURBMODEL_SST_H
#define RANSTURBMODEL_SST_H

#include "UgpWithCvCompFlow.h"

class RansTurbSST : virtual public UgpWithCvCompFlow
{
public:

  RansTurbSST()
  {
    if (mpi_rank == 0)
      cout << "RansTurbSST()" << endl;

    turbModel = SST;

    sigma_k1  = getDoubleParam("sigma_k1" , "0.85"   );
    sigma_k2  = getDoubleParam("sigma_k2" , "1.0"    );
    sigma_om1 = getDoubleParam("sigma_om1", "0.5"    );
    sigma_om2 = getDoubleParam("sigma_om2", "0.856"  );
    beta_1    = getDoubleParam("beta_1"   , "0.075"  );
    beta_2    = getDoubleParam("beta_2"   , "0.0828" );
    betaStar  = getDoubleParam("betaStar" , "0.09"   );
    a1        = getDoubleParam("a1"       , "0.31"   );

    sst_form = getStringParam("SST_FORM", "STANDARD");

    if (sst_form == "STANDARD")
    {
      gamma_1 = beta_1/betaStar - sigma_om1*pow(0.41, 2.0)/sqrt(betaStar);
      gamma_2 = beta_2/betaStar - sigma_om2*pow(0.41, 2.0)/sqrt(betaStar);
      boundPk = 20.0;
      boundPw = 1.0e+12; // no effective limiter
      boundCrossDiff = 1.0e-20;
    }
    else if (sst_form == "2003")
    {
      gamma_1 = getDoubleParam("gamma_1"  , "0.55555");
      gamma_2 = getDoubleParam("gamma_2"  , "0.44"   );
      boundPk = 10.0;
      boundPw = 10.0;
      boundCrossDiff = 1.0e-10;
    }
    else
    {
      cerr << "ERROR: SST model form specified not implemented." << endl;
      throw(-1);
    }

    ScalarTranspEq *eq;
    eq = registerScalarTransport("kine", CV_DATA);
    eq->relax = getDoubleParam("RELAX_kine", "0.6");
    eq->phiZero = getDoubleParam("ZERO_kine", "1.0e-8");
    eq->phiZeroRel = getDoubleParam("ZERO_REL_kine", "1.0e-2");
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0e10;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    eq = registerScalarTransport("omega", CV_DATA);
    eq->relax = getDoubleParam("RELAX_omega", "0.4");
    eq->phiZero = getDoubleParam("ZERO_omega", "1.0e-8");
    eq->phiZeroRel = getDoubleParam("ZERO_REL_omega", "1.0e-2");
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-4;
    eq->upperBound = 1.0e15;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    strMag      = NULL;   registerScalar(strMag     , "strMag"     , CV_DATA);
    vortMag     = NULL;   registerScalar(vortMag    , "vortMag", CV_DATA);
    diverg      = NULL;   registerScalar(diverg     , "diverg"     , CV_DATA);
    muT         = NULL;   registerScalar(muT        , "muT"        , CV_DATA);
    crossDiff   = NULL;   registerScalar(crossDiff  , "crossDiff"  , CV_DATA);
    blendFuncF1 = NULL;   registerScalar(blendFuncF1, "blendFuncF1", CV_DATA);
    blendFuncF2 = NULL;   registerScalar(blendFuncF2, "blendFuncF2", CV_DATA);
    wallDist    = NULL;   registerScalar(wallDist   , "wallDist"   , CV_DATA);
    wallConn    = NULL;   // array of integers
  }

  virtual ~RansTurbSST() {}

public:

  // model variables
  int *wallConn;           ///< index of closest wall face

  double *omega;           ///< specific dissipation
  double (*grad_kine)[3];  ///< gradient of tke
  double (*grad_omega)[3]; ///< gradient of omega
  double *kine_bfa;        ///< tke at the boundaries
  double *omega_bfa;       ///< omega at the boundaries
  double *muT;             ///< turbulent viscosity at cell center for output
  double *wallDist;        ///< distance to closest wall face
  double *crossDiff;       ///< cross-diffusion
  double *blendFuncF1;     ///< first blending function
  double *blendFuncF2;     ///< second blending function
  double *limiterFunc;     ///< vort. vs. strain visc. limiter

  string sst_form;         ///< various model forms: standard, 2003, etc.
  
  // model constants
  double sigma_k1, sigma_k2, sigma_om1, sigma_om2;
  double gamma_1, gamma_2, beta_1, beta_2;
  double betaStar, a1, boundPk, boundPw, boundCrossDiff;

public:

  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0) 
      cout << "initialHookScalarRansTurbModel()" << endl;

    wallConn = new int[ncv];
    calcWallDistance(wallConn, wallDist);

    if      (sst_form == "STANDARD") limiterFunc = vortMag;
    else if (sst_form == "2003")     limiterFunc = strMag;

    // connect pointers
    ScalarTranspEq *eq;
    eq = getScalarTransportData("kine");      kine = eq->phi;    kine_bfa = eq->phi_bfa;   grad_kine = eq->grad_phi;
    eq = getScalarTransportData("omega");     omega = eq->phi;   omega_bfa = eq->phi_bfa;  grad_omega = eq->grad_phi;
  }
  
  virtual void calcRansTurbViscMuet()
  {
    // update velocity gradients, cross-diffusion, blending functions
    calcGradVel();
    calcStrainRateAndDivergence();
    calcVorticity();
    calcMenterBlendingFunctions();

    // eddy-viscosity at cell center
    for (int icv=0; icv<ncv; icv++)
    {
      double zeta = min(1.0/omega[icv], a1/(limiterFunc[icv]*blendFuncF2[icv]));
      muT[icv] = min(max(rho[icv]*kine[icv]*zeta, 0.0), 1.0);
    }
    updateCvData(muT, REPLACE_DATA);

    // internal faces
    for (int ifa=nfa_b; ifa<nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];

      double dx0[3], dx1[3];
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));

      mut_fa[ifa] = (w1*muT[icv0] + w0*muT[icv1])/(w0+w1);
    }

    // boundary faces
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            mut_fa[ifa] = 0.0;                      // set mut zero at walls
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            mut_fa[ifa] = muT[icv0];
          }
      }
  }

  virtual void diffusivityHookScalarRansTurb(const string &name)
  {
    ScalarTranspEq *eq;

    if (name == "kine")
    {
      eq = getScalarTransportData(name);
      
      // internal faces
      for (int ifa=nfa_b; ifa<nfa; ifa++)
      {
        int icv0 = cvofa[ifa][0];
        int icv1 = cvofa[ifa][1];

        double dx0[3], dx1[3];
        vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
        vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
        double w0 = sqrt(vecDotVec3d(dx0, dx0));
        double w1 = sqrt(vecDotVec3d(dx1, dx1));

        double f1_fa = (w1*blendFuncF1[icv0] + w0*blendFuncF1[icv1])/(w0+w1);

        double coeff = sigma_k1*f1_fa + (1. - f1_fa)*sigma_k2;
        eq->diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
      }

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = mul_fa[ifa];
          }
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            double coeff = sigma_k1*blendFuncF1[icv0] + (1.0 - blendFuncF1[icv0])*sigma_k2;
            eq->diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
          }
      }
    }

    if (name == "omega")
    {
      eq = getScalarTransportData(name);

      // internal faces
      for (int ifa=nfa_b; ifa<nfa; ifa++)
      {
        int icv0 = cvofa[ifa][0];
        int icv1 = cvofa[ifa][1];

        double dx0[3], dx1[3];
        vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
        vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
        double w0 = sqrt(vecDotVec3d(dx0, dx0));
        double w1 = sqrt(vecDotVec3d(dx1, dx1));

        double f1_fa = (w1*blendFuncF1[icv0] + w0*blendFuncF1[icv1])/(w0+w1);

        double coeff = sigma_om1*f1_fa + (1.0 - f1_fa)*sigma_om2;
        eq->diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
      }

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = mul_fa[ifa];
          }
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            double coeff = sigma_om1*blendFuncF1[icv0] + (1. - blendFuncF1[icv0])*sigma_om2;
            eq->diff[ifa] = mul_fa[ifa] + coeff*mut_fa[ifa];
          }
      }
    }
  }

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")
      for (int icv=0; icv<ncv; icv++)
      {
        double Pk = muT[icv]*strMag[icv]*strMag[icv] - 2.0/3.0*rho[icv]*kine[icv]*diverg[icv];

        Pk = min(Pk, boundPk*betaStar*rho[icv]*kine[icv]*omega[icv]);
        Pk = max(Pk, 0.0);

        double src = Pk - betaStar*rho[icv]*omega[icv]*kine[icv];
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = - betaStar*omega[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }

    if (name == "omega")
    {
      for (int icv=0; icv<ncv; icv++)
      {
        double F1 = blendFuncF1[icv];
        double alfa = F1*gamma_1 + (1.0 - F1)*gamma_2;
        double beta = F1*beta_1 + (1.0 - F1)*beta_2;

        double zeta = max(omega[icv], limiterFunc[icv]*blendFuncF2[icv]/a1);
        double Pw = strMag[icv]*strMag[icv] - 2.0/3.0*zeta*diverg[icv];

        Pw = min(Pw, boundPw*betaStar*omega[icv]*zeta);
        Pw = max(Pw, 0.0);

        double src = alfa*rho[icv]*Pw + (1.0 - F1)*crossDiff[icv] - beta*rho[icv]*omega[icv]*omega[icv];
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi =  - 2.0*beta*omega[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    }
  }

  virtual void sourceHookRansTurbCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
    int kine_Index = getScalarTransportIndex("kine");
    int omega_Index = getScalarTransportIndex("omega");
    
    for (int icv = 0; icv < ncv; icv++)
    {
      double F1 = blendFuncF1[icv];
      double alfa = F1 * gamma_1 + (1.0 - F1) * gamma_2;
      double beta = F1 * beta_1 + (1.0 - F1) * beta_2;

      double Pk = muT[icv]*strMag[icv]*strMag[icv] - 2.0/3.0*rho[icv]*kine[icv]*diverg[icv];

      Pk = min(Pk, boundPk*betaStar*rho[icv]*kine[icv]*omega[icv]);
      Pk = max(Pk, 0.0);

      double src = Pk - betaStar*rho[icv]*omega[icv]*kine[icv];
      rhs[icv][5+kine_Index] += src*cv_volume[icv];

      double zeta = max(omega[icv], limiterFunc[icv]*blendFuncF2[icv]/a1);
      double Pw = strMag[icv]*strMag[icv] - 2.0/3.0*zeta*diverg[icv];

      Pw = min(Pw, boundPw*betaStar*omega[icv]*zeta);
      Pw = max(Pw, 0.0);

      src = alfa*rho[icv]*Pw + (1.0 - F1)*crossDiff[icv] - beta*rho[icv]*omega[icv]*omega[icv];
      rhs[icv][5+omega_Index] += src*cv_volume[icv];

      if (flagImplicit)
      {
        int noc00 = nbocv_i[icv];
        A[noc00][5+kine_Index][5+kine_Index]   -= - betaStar*omega[icv]*cv_volume[icv];
        A[noc00][5+omega_Index][5+omega_Index] -= - 2.0*beta*omega[icv]*cv_volume[icv];
     }
    }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_fa, FaZone *zone, const string &name)
  {
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      Param *param;
      if (getParam(param, zone->getName()))
        if ((param->getString() == "WALL") && (name == "omega"))
        {
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv = cvofa[ifa][0];
            double muLamCV = calcMuLam(icv);
            phi_fa[ifa] = 60.0*muLamCV/(rho[icv]*beta_1*wallDist[icv]*wallDist[icv]);
          }
        }
    }
  }

  virtual void calcMenterBlendingFunctions()
  {
    calcCvScalarGrad(grad_kine,  kine,  kine_bfa,  gradreconstruction, limiterNavierS, kine,  epsilonSDWLS);
    calcCvScalarGrad(grad_omega, omega, omega_bfa, gradreconstruction, limiterNavierS, omega, epsilonSDWLS);

    for (int icv=0; icv<ncv; icv++)
    {
      double d = wallDist[icv];
      double mue = InterpolateAtCellCenterFromFaceValues(mul_fa, icv);

      crossDiff[icv] = max(2.0*rho[icv]*sigma_om2/omega[icv]*vecDotVec3d(grad_kine[icv], grad_omega[icv]), boundCrossDiff);

      double gamma1 = 500.0*mue/(pow(d, 2.0)*rho[icv]*omega[icv]);
      double gamma2 = 4.0*sigma_om2*rho[icv]*kine[icv]/(d*d*crossDiff[icv]);
      double gamma3 = sqrt(kine[icv])/(betaStar*omega[icv]*d);

      double gamma = min(max(gamma1, gamma3), gamma2);
      blendFuncF1[icv] = tanh(pow(gamma,4.0));
      gamma = max(2.0*gamma3, gamma1);
      blendFuncF2[icv] = tanh(pow(gamma,2.0));
    }

    updateCvData(crossDiff, REPLACE_DATA);
    updateCvData(blendFuncF1, REPLACE_DATA);
    updateCvData(blendFuncF2, REPLACE_DATA);
  }

  virtual void finalHookScalarRansTurbModel()
  {
    if (wallConn != NULL) {delete [] wallConn; wallConn = NULL;}

    if (checkParam("CALC_RS_BOUSS"))
    {
      if (mpi_rank == 0)
        cout << "calculating Rs Bouss" << endl;
      for (int icv = 0; icv < ncv; icv++)
      {
        double term1 = (2.0/3.0) * rho[icv] * kine[icv];
        rij_diag[icv][0] = -term1 + muT[icv] * 2.0 * (grad_u[icv][0][0] - 1.0/3.0*diverg[icv]);
        rij_diag[icv][1] = -term1 + muT[icv] * 2.0 * (grad_u[icv][1][1] - 1.0/3.0*diverg[icv]);
        rij_diag[icv][2] = -term1 + muT[icv] * 2.0 * (grad_u[icv][2][2] - 1.0/3.0*diverg[icv]);
        rij_offdiag[icv][0] =     + muT[icv] * (grad_u[icv][0][1] + grad_u[icv][1][0]);
        rij_offdiag[icv][1] =     + muT[icv] * (grad_u[icv][0][2] + grad_u[icv][2][0]);
        rij_offdiag[icv][2] =     + muT[icv] * (grad_u[icv][1][2] + grad_u[icv][2][1]);
      }
    }
  }

};

#endif


