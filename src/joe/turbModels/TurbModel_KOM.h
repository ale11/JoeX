#ifndef RANSTURBMODEL_KOM_H
#define RANSTURBMODEL_KOM_H

#include "UgpWithCvCompFlow.h"

//######################################################//
//                                                      //
// Wilcox k-omega Two-Equation Model (Wilcox2006)       //
//                                                      //
//######################################################//

class RansTurbKOm : virtual public UgpWithCvCompFlow
{
public:   // constructors

  RansTurbKOm()
  {
    if (mpi_rank == 0)
      cout << "RansTurbKOm()" << endl;
    
    turbModel = KOM;

    betaStar = getDoubleParam("betaStar", "0.09");
    sigma_k  = getDoubleParam("sigma_k",  "0.6");
    sigma_om = getDoubleParam("sigma_om", "0.5");
    alfa     = getDoubleParam("alfa",     "0.52");
    sigmad0  = getDoubleParam("sigmad0",  "0.125");
    beta0    = getDoubleParam("beta0",    "0.0708");
    cLim     = getDoubleParam("cLim",     "0.875");

    ScalarTranspEq *eq;

    eq = registerScalarTransport("kine", CV_DATA);
    eq->relax = getDoubleParam("RELAX_kine", "0.4");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0e10;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    eq = registerScalarTransport("omega", CV_DATA);
    eq->relax = getDoubleParam("RELAX_omega", "0.4");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-4;
    eq->upperBound = 1.0e15;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    strMag   = NULL;       registerScalar(strMag, "strMag", CV_DATA);
    vortMag  = NULL;       registerScalar(vortMag, "vortMag", CV_DATA);
    diverg   = NULL;       registerScalar(diverg, "diverg", CV_DATA);
    muT      = NULL;       registerScalar(muT, "muT", CV_DATA);
    wallDist = NULL;       registerScalar(wallDist, "wallDist",  CV_DATA);
  }

  virtual ~RansTurbKOm() {}

public:

  // model variables
  double *omega;           ///< specific dissipation
  double (*grad_kine)[3];  ///< gradient of tke
  double (*grad_omega)[3]; ///< gradient of omega
  double *kine_bfa;        ///< tke at the boundaries
  double *omega_bfa;       ///< omega at the boundaries
  double *muT;             ///< turbulent viscosity at cell center
  double *wallDist;        ///< wall distance
  
  // model constants
  double betaStar, sigma_om, sigma_k, alfa, sigmad0, beta0, cLim;

public:
  
  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0) 
      cout << "initialHook WILCOX KOM model" << endl;

    if (!checkDataFlag(wallDist))
    {
      for (int icv=0; icv<ncv; icv++)
        wallDist[icv] = 0.0;
      calcWallDistance(NULL, wallDist);
    }

    // connect pointers
    ScalarTranspEq *eq;
    eq = getScalarTransportData("kine");    kine = eq->phi;     kine_bfa = eq->phi_bfa;   grad_kine = eq->grad_phi;
    eq = getScalarTransportData("omega");   omega = eq->phi;    omega_bfa = eq->phi_bfa;  grad_omega = eq->grad_phi;
  }

  virtual void calcRansTurbViscMuet()
  {
    // update velocity gradients
    calcGradVel();
    calcStrainRateAndDivergence();
    calcVorticity();

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

      double rho_fa = (w1*rho[icv0] + w0*rho[icv1])/(w0+w1);
      double kine_fa = (w1*kine[icv0] + w0*kine[icv1])/(w0+w1);
      double om_fa = (w1*omega[icv0] + w0*omega[icv1])/(w0+w1);
      double strMag_fa = (w1*strMag[icv0] + w0*strMag[icv1])/(w0+w1);

      double omega_tilde = max(om_fa, cLim*strMag_fa/sqrt(betaStar));
      mut_fa[ifa] = min(rho_fa*kine_fa/omega_tilde, 100.0);
    }

    // boundary faces
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    if (zone->getKind() == FA_ZONE_BOUNDARY)
    {
      if (zoneIsWall(zone->getName()))
        for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          mut_fa[ifa] = 0.0;                                     // set mut zero at walls
      else
        for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
        {
          int icv0 = cvofa[ifa][0];
          
          double omega_tilde = max(omega[icv0], cLim*strMag[icv0]/sqrt(betaStar));
          mut_fa[ifa] = min(rho[icv0]*kine[icv0]/omega_tilde, 100.0);    // zero order extrapolation for others
        }
    }

    // just calculate mut in cell center by averaging over all face values
    for (int icv=0; icv<ncv; icv++)
      muT[icv] = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);
  }

  virtual void diffusivityHookScalarRansTurb(const string &name)
  {
    ScalarTranspEq *eq;

    if ((name == "kine") || (name == "omega"))
    {
      double sigma;
      if (name == "kine")       sigma = sigma_k;
      if (name == "omega")      sigma = sigma_om;

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

        double rho_fa = (w1*rho[icv0] + w0*rho[icv1])/(w0+w1);
        double kine_fa = (w1*kine[icv0] + w0*kine[icv1])/(w0+w1);
        double om_fa = (w1*omega[icv0] + w0*omega[icv1])/(w0+w1);

        eq->diff[ifa] = mul_fa[ifa] + sigma*rho_fa*kine_fa/om_fa;
      }

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
        {
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = mul_fa[ifa];            // set mut zero at walls
          }
        }
        else
        {
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            eq->diff[ifa] = mul_fa[ifa] + sigma*rho_bfa[ifa]*kine_bfa[ifa]/omega_bfa[ifa];
          }
        }
      }
    }
  }

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")
      for (int icv=0; icv<ncv; icv++)
      {
        double src = calcTurbProd(icv) - betaStar*rho[icv]*omega[icv]*kine[icv];
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
      double OM[3][3], STR_hat[3][3];

      calcCvScalarGrad(grad_kine,  kine,  kine_bfa,  gradreconstruction, limiterNavierS, kine,  epsilonSDWLS);
      calcCvScalarGrad(grad_omega, omega, omega_bfa, gradreconstruction, limiterNavierS, omega, epsilonSDWLS);

      for (int icv=0; icv<ncv; icv++)
      {
        // cross-diffusion
        double sigmad;
        double crossDiff = vecDotVec3d(grad_kine[icv], grad_omega[icv]);
        if (crossDiff <= 0.0) sigmad = 0.0;
        else                  sigmad = sigmad0;

        // rate of strain and rotation tensors
        for (int i=0; i<3; i++)
          for (int j=0; j<3; j++)
          {
            OM[i][j] = 0.5*(grad_u[icv][i][j] - grad_u[icv][j][i]);
            // STR_hat subtracts 0.5 divergence instead of 1/3!!!
            if (i==j)   STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i] - diverg[icv]);
            else        STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]);
          }

        // chi_omega
        double chiOm = 0.0;
        for (int i=0; i<3; i++)
          for (int j=0; j<3; j++)
            for (int k=0; k<3; k++)
              chiOm += OM[i][j]*OM[j][k]*STR_hat[k][i];
        chiOm = fabs(chiOm/pow(betaStar*omega[icv], 3.0));

        double fbeta = (1.0 + 85.0*chiOm)/(1.0 + 100.0*chiOm);
        double beta = beta0*fbeta;

        double src =  alfa*omega[icv]/kine[icv]*calcTurbProd(icv)
                    - beta*rho[icv]*omega[icv]*omega[icv]
                    + sigmad*rho[icv]/omega[icv]*crossDiff;

        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dsrcdphi = -2.0*beta*omega[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    }
  }

  virtual void sourceHookRansTurbCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
    int kine_Index = getScalarTransportIndex("kine");
    int omega_Index = getScalarTransportIndex("omega");
    
    double OM[3][3], STR_hat[3][3];

    calcCvScalarGrad(grad_kine,  kine,  kine_bfa,  gradreconstruction, limiterNavierS, kine,  epsilonSDWLS);
    calcCvScalarGrad(grad_omega, omega, omega_bfa, gradreconstruction, limiterNavierS, omega, epsilonSDWLS);

    for (int icv = 0; icv < ncv; icv++)
    {
      // cross-diffusion
      double sigmad;
      double crossDiff = vecDotVec3d(grad_kine[icv], grad_omega[icv]);
      if (crossDiff <= 0.0) sigmad = 0.0;
      else                  sigmad = sigmad0;

      // rate of strain and rotation tensors
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
        {
          OM[i][j] = 0.5*(grad_u[icv][i][j] - grad_u[icv][j][i]);

          // STR_hat subtracts 0.5 divergence instead 1/3, no bug!!!
          if (i==j)   STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i] - diverg[icv]);
          else        STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]);
        }

      // chi_omega
      double chiOm = 0.0;
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
          for (int k=0; k<3; k++)
            chiOm += OM[i][j]*OM[j][k]*STR_hat[k][i];
      chiOm = fabs(chiOm/pow(betaStar*omega[icv], 3.0));

      double fbeta = (1.0 + 85.0*chiOm)/(1.0 + 100.0*chiOm);
      double beta = beta0*fbeta;

      double src =  alfa*omega[icv]/kine[icv]*calcTurbProd(icv)
                  - beta*rho[icv]*omega[icv]*omega[icv]
                  + sigmad*rho[icv]/omega[icv]*crossDiff;

      rhs[icv][5+kine_Index]  += (calcTurbProd(icv) - betaStar*rho[icv]*omega[icv]*kine[icv])*cv_volume[icv];
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
            phi_fa[ifa] = 6.0*muLamCV/(rho[icv]*beta0*wallDist[icv]*wallDist[icv]);
          }
        }
    }
  }

  virtual double calcTurbProd(int icv)
  {
    double omega_tilde = max(omega[icv], cLim*strMag[icv]/sqrt(betaStar));
    double mu_t = rho[icv]*kine[icv]/omega_tilde;
    return max(mu_t*strMag[icv]*strMag[icv] - 2.0/3.0*rho[icv]*kine[icv]*diverg[icv], 0.0);
  }

};



#endif
