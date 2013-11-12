#ifndef RANSTURBMODEL_EASM_H
#define RANSTURBMODEL_EASM_H

#include "UgpWithCvCompFlow.h"

class RansTurbEASM : virtual public UgpWithCvCompFlow
{
public:
  RansTurbEASM()
  {
    cmus   = NULL;       registerScalar(cmus,   "cmus",   CV_DATA);
  }

  virtual ~RansTurbEASM() {}

public:
  double *cmus; ///< variable CMU

public:
  virtual void calcRsCenterEASM()
  {
    double ceps1 = 1.44;
    double ceps2 = 1.83;
    double c10   = 3.4;
    double c11   = 1.8;
    double c2    = 0.36;
    double c3    = 1.25;
    double c4    = 0.4;

    double gamma0 = 0.5*c11;
    double gamma1 = 0.5*c10 + (ceps2 - ceps1)/(ceps1 - 1.0);
    double a1     = 0.5*(4./3. - c2);
    double a2     = 0.5*(2.0 - c4);
    double a3     = 0.5*(2.0 - c3);

    double S11, S12, S13, S21, S22, S23, S31, S32, S33;
    double W11, W12, W13, W21, W22, W23, W31, W32, W33;
    double SijSij, WijWij;
    double S1kWk1, W1kSk1, S1kSk1;
    double S1kWk2, W1kSk2, S1kSk2;
    double S2kWk2, W2kSk2, S2kSk2;
    double S3kWk3, W3kSk3, S3kSk3;
    double S1kWk3, W1kSk3, S1kSk3;
    double S2kWk3, W2kSk3, S2kSk3;

    double tau, e2t2, e2t2G0, R2e2t2, p, q, r, a4;
    double cmus_temp, cmustau;
    double a, b, d, t1, t2, t3, theta;
    double PI = 3.141592653589793238462;

    double r11, r22, r33, r12, r13, r23;

    // Compute cmus and Reynolds stresses
    for (int icv = 0; icv < ncv; icv++)
    {
      // Rate of strain tensor
      S11 = grad_u[icv][0][0];
      S12 = 0.5*(grad_u[icv][0][1] + grad_u[icv][1][0]);
      S13 = 0.5*(grad_u[icv][0][2] + grad_u[icv][2][0]);

      S21 = S12;
      S22 = grad_u[icv][1][1];
      S23 = 0.5*(grad_u[icv][1][2] + grad_u[icv][2][1]);

      S31 = S13;
      S32 = S23;
      S33 = grad_u[icv][2][2];

      // Rate of mean rotation tensor
      W11 = 0.0;
      W12 = 0.5*(grad_u[icv][0][1] - grad_u[icv][1][0]);
      W13 = 0.5*(grad_u[icv][0][2] - grad_u[icv][2][0]);

      W21 = -W12;
      W22 = 0.0;
      W23 = 0.5*(grad_u[icv][1][2] - grad_u[icv][2][1]);

      W31 = -W13;
      W32 = -W23;
      W33 = 0.0;

      SijSij = S11*S11 + S22*S22 + S33*S33 + 2.0*S12*S12 + 2.0*S13*S13 + 2.0*S23*S23;
      WijWij = 2.*W12*W12 + 2.*W13*W13 + 2.*W23*W23;

      S1kWk1 = S11*W11 + S12*W21 + S13*W31;
      W1kSk1 = W11*S11 + W12*S21 + W13*S31;
      S1kSk1 = S11*S11 + S12*S21 + S13*S31;

      S2kWk2 = S21*W12 + S22*W22 + S23*W32;
      W2kSk2 = W21*S12 + W22*S22 + W23*S32;
      S2kSk2 = S21*S12 + S22*S22 + S23*S32;

      S3kWk3 = S31*W13 + S32*W23 + S33*W33;
      W3kSk3 = W31*S13 + W32*S23 + W33*W33;
      S3kSk3 = S31*S13 + S32*S23 + S33*S33;

      S1kWk2 = S11*W12 + S12*W22 + S13*W32;
      W1kSk2 = W11*S12 + W12*S22 + W13*S32;
      S1kSk2 = S11*S12 + S12*S22 + S13*S32;

      S1kWk3 = S11*W13 + S12*W23 + S13*W33;
      W1kSk3 = W11*S13 + W12*S23 + W13*S33;
      S1kSk3 = S11*S13 + S12*S23 + S13*S33;

      S2kWk3 = S21*W13 + S22*W23 + S23*W33;
      W2kSk3 = W21*S13 + W22*S23 + W23*S33;
      S2kSk3 = S21*S13 + S22*S23 + S23*S33;

      tau    = turbTS[icv];
      e2t2   = SijSij*tau*tau;
      R2e2t2 = WijWij*tau*tau;
      e2t2G0 = e2t2*gamma0;

      p = -gamma1/(e2t2G0);
      q = (gamma1*gamma1 - 2.0*e2t2G0*a1 - 2.0/3.0*e2t2*a3*a3 + 2.0*R2e2t2*a2*a2);
      q /= 4.0*e2t2G0*e2t2G0;
      r = gamma1*a1/(4.0*e2t2G0*e2t2G0);

      // Compute Cmustar
      if (e2t2 < 1.0e-6)
        // degenerate case
        cmus_temp = gamma1*a1/(gamma1*gamma1 + 2.*R2e2t2*a2*a2);
      else
      {
        // find the roots
        a = q - p*p/3.;
        b = 1.0/27.0*(2*p*p*p - 9.0*p*q + 27.0*r);
        d = (b*b)/4.0 + (a*a*a)/27.0;
        if (d > 0)
        {
          t1 = -0.5*b + sqrt(d);
          if (t1 > 0.0)
            t1 = pow(fabs(t1), 1.0/3.0);
          else
            t1 = -pow(fabs(t1), 1.0/3.0);

          t2 = -0.5*b - sqrt(d);
          if (t2 > 0.0)
            t2 = pow(fabs(t2), 1.0/3.0);
          else
            t2 = -pow(fabs(t2), 1.0/3.0);

          cmus_temp = -min(-p/3.0 + t1 + t2, -p/3.0 - 0.5*t1 - 0.5*t2);
        }
        else
        {
          theta = acos(-0.5*b/sqrt(-a*a*a/27.0));
          t1 = -p/3.0 + 2.0*sqrt(-a/3.0)*cos(theta/3.);
          t2 = -p/3.0 + 2.0*sqrt(-a/3.0)*cos(2.0/3.0*PI + theta/3.0);
          t3 = -p/3.0 + 2.0*sqrt(-a/3.0)*cos(4.0/3.0*PI + theta/3.0);
          cmus_temp = -min(min(t1,t2),t3);
        }
      }

      cmus_temp = max(cmus_temp,0.0005);
      cmus[icv] = cmus_temp;

      cmustau = cmus_temp*tau;
      a4 = tau/(gamma1 + 2.0*gamma0*cmus_temp*e2t2);

      // Compute the Reynolds stresses
      r11 = 1.0/3.0 - cmustau*(S11 - 1.0/3.0*diverg[icv] + a2*a4*(S1kWk1 - W1kSk1) - 2.0*a3*a4*(S1kSk1 - SijSij/3.0));
      r22 = 1.0/3.0 - cmustau*(S22 - 1.0/3.0*diverg[icv] + a2*a4*(S2kWk2 - W2kSk2) - 2.0*a3*a4*(S2kSk2 - SijSij/3.0));
      r33 = 1.0/3.0 - cmustau*(S33 - 1.0/3.0*diverg[icv] + a2*a4*(S3kWk3 - W3kSk3) - 2.0*a3*a4*(S3kSk3 - SijSij/3.0));
      r12 = -cmustau*(S12 + a2*a4*(S1kWk2 - W1kSk2) - 2.0*a3*a4*S1kSk2);
      r13 = -cmustau*(S13 + a2*a4*(S1kWk3 - W1kSk3) - 2.0*a3*a4*S1kSk3);
      r23 = -cmustau*(S23 + a2*a4*(S2kWk3 - W2kSk3) - 2.0*a3*a4*S2kSk3);

      rij_diag[icv][0] = -max(r11,0.)*2*kine[icv]*rho[icv];
      rij_diag[icv][1] = -max(r22,0.)*2*kine[icv]*rho[icv];
      rij_diag[icv][2] = -max(r33,0.)*2*kine[icv]*rho[icv];

      rij_offdiag[icv][0] = -r12*2*kine[icv]*rho[icv];
      rij_offdiag[icv][1] = -r13*2*kine[icv]*rho[icv];
      rij_offdiag[icv][2] = -r23*2*kine[icv]*rho[icv];
    }

    updateCvData(cmus, REPLACE_DATA);
    updateCvData(rij_diag, REPLACE_ROTATE_DATA);
    updateCvData(rij_offdiag, REPLACE_ROTATE_DATA);

    // Interpolate cell centered stresses to cell faces
    interpolateStressToFace();
  }

  virtual void interpolateStressToFace()
  {
    // ====================================================================
    // Compute internal face Reynolds stresses
    // ====================================================================
    for (int ifa = nfa_b; ifa < nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];
      assert( icv0 >= 0 );
      assert( icv1 >= 0 );

      double dx0[3], dx1[3];
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));

      // Update Reynolds stresses
      rij_diag_fa[ifa][0] = (w1*rij_diag[icv0][0] + w0*rij_diag[icv1][0])/(w0 + w1);
      rij_diag_fa[ifa][1] = (w1*rij_diag[icv0][1] + w0*rij_diag[icv1][1])/(w0 + w1);
      rij_diag_fa[ifa][2] = (w1*rij_diag[icv0][2] + w0*rij_diag[icv1][2])/(w0 + w1);

      rij_offdiag_fa[ifa][0] = (w1*rij_offdiag[icv0][0] + w0*rij_offdiag[icv1][0])/(w0 + w1);
      rij_offdiag_fa[ifa][1] = (w1*rij_offdiag[icv0][1] + w0*rij_offdiag[icv1][1])/(w0 + w1);
      rij_offdiag_fa[ifa][2] = (w1*rij_offdiag[icv0][2] + w0*rij_offdiag[icv1][2])/(w0 + w1);
    }

    // ====================================================================
    // Compute boundary face Reynolds stresses
    // ====================================================================
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    {
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
        {
          // .............................................................................................
          // SYMMETRY BOUNDARY CONDITION
          // .............................................................................................
          if (param->getString() == "SYMMETRY")
          {
            // No viscous flux in this case (or yes!)

            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              int icv0 = cvofa[ifa][0];
              assert( icv0 >= 0 );

              // define Reynolds stresses
              rij_diag_fa[ifa][0] = rij_diag[icv0][0];
              rij_diag_fa[ifa][1] = rij_diag[icv0][1];
              rij_diag_fa[ifa][2] = rij_diag[icv0][2];

              rij_offdiag_fa[ifa][0] = rij_offdiag[icv0][0];
              rij_offdiag_fa[ifa][1] = rij_offdiag[icv0][1];
              rij_offdiag_fa[ifa][2] = rij_offdiag[icv0][2];

            }
          }
          // .............................................................................................
          // WALL BOUNDARY CONDITION
          // .............................................................................................
          else if (param->getString() == "WALL")
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              rij_diag_fa[ifa][0] = 0.0;
              rij_diag_fa[ifa][1] = 0.0;
              rij_diag_fa[ifa][2] = 0.0;

              rij_offdiag_fa[ifa][0] = 0.0;
              rij_offdiag_fa[ifa][1] = 0.0;
              rij_offdiag_fa[ifa][2] = 0.0;
            }
          }
          // .............................................................................................
          // OTHER BOUNDARY CONDITIONS
          // .............................................................................................
          else
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              int icv0 = cvofa[ifa][0];
              assert( icv0 >= 0 );

              // define Reynolds stresses
              rij_diag_fa[ifa][0] = rij_diag[icv0][0];
              rij_diag_fa[ifa][1] = rij_diag[icv0][1];
              rij_diag_fa[ifa][2] = rij_diag[icv0][2];

              rij_offdiag_fa[ifa][0] = rij_offdiag[icv0][0];
              rij_offdiag_fa[ifa][1] = rij_offdiag[icv0][1];
              rij_offdiag_fa[ifa][2] = rij_offdiag[icv0][2];
            }
          }
        }
      }
    }
  }

};

/*
 * Nonlinear EASM k-omega (2003) Model (EASMko2003)
 */

// differences with standard EASM:
// muT is being limited between zero and one.
// -for production of w, I divide by a limited k/w.
// -the production of w (instead of just k) is also being limited.

class RansTurbEASMkom : virtual public RansTurbEASM
{
public:
  RansTurbEASMkom()
  {
    if (mpi_rank == 0)
      cout << "RansTurbEASMkom()" << endl;
    
    turbModel = EASMkom;

    C_MU = getDoubleParam("C_MU", "0.0895");
    BETA = getDoubleParam("BETA", "0.83");
    G_OM = getDoubleParam("G_OM", "0.53");

    SIG_K = getDoubleParam("SIG_K"  , "1.0");
    SIG_O = 0.41*0.41/(sqrt(C_MU)*(BETA - G_OM));

    LIMIT_PK   = getIntParam("LIMIT_PK"  , "1");
    REALIZABLE = getIntParam("REALIZABLE", "0");

    ScalarTranspEq *eq;
    eq = registerScalarTransport("kine",  CV_DATA);
    eq->relax = getDoubleParam("RELAX_kine", "0.7");
    eq->phiZero = 1.0e-8;
    eq->phiMaxiter = 500;
    eq->lowerBound = 1.0e-6;
    eq->upperBound = 1.0e10;
    eq->turbSchmidtNumber = SIG_K;

    eq = registerScalarTransport("omega", CV_DATA);
    eq->relax = getDoubleParam("RELAX_omega", "0.4");
    eq->phiZero = getDoubleParam("ZERO_omega", "1.0e-8");
    eq->phiZeroRel = getDoubleParam("ZERO_REL_omega", "1.0e-2");
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-4;
    eq->upperBound = 1.0e15;
    eq->turbSchmidtNumber = SIG_O;

    strMag = NULL;       registerScalar(strMag, "strMag", CV_DATA);
    diverg = NULL;       registerScalar(diverg, "diverg", CV_DATA);
    muT    = NULL;       registerScalar(muT,    "muT",    CV_DATA);
    fbstar = NULL;       registerScalar(fbstar, "fbstar", CV_DATA);
    turbTS = NULL;       registerScalar(turbTS, "turbTS", CV_DATA);
  }

  virtual ~RansTurbEASMkom() {}

public:

  double *omega;                            ///< introduced to have access to variables
  double *kine_bfa, *omega_bfa;             ///< turbulent scalars at the boundary
  double (*grad_kine)[3], (*grad_omega)[3]; ///< gradients for turb variables
  double *muT;                              ///< turbulent viscosity at cell center for output
  double *fbstar;
  
  double C_MU, SIG_K, SIG_O, BETA, G_OM; ///< model constants

  int LIMIT_PK;   ///< limiter for tke production
  int REALIZABLE; ///< realizability limits for turbTS, turbLS

public:

  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0)
      cout << "initialHookScalarRansTurbModel()" << endl;

    ScalarTranspEq *eq;
    eq = getScalarTransportData("kine");   kine = eq->phi;    kine_bfa = eq->phi_bfa;     grad_kine = eq->grad_phi;
    eq = getScalarTransportData("omega");  omega = eq->phi;   omega_bfa = eq->phi_bfa;    grad_omega = eq->grad_phi;

    for (int ifa = 0; ifa < nfa; ifa++)
      nonLinear[ifa] = 1.0;
    if (mpi_rank == 0)
      cout << "NON-LINEAR DOMAIN ACTIVATED." << endl;
  }

  virtual void calcRansTurbViscMuet()
  {
    // update velocity gradients, fbstar, timescale
    calcGradVel();
    calcStrainRateAndDivergence();
    calcFbstarFunction();
    calcTurbTimeScale();

    // compute Reynolds stresses, cmus
    calcRsCenterEASM();

    // Calculate turbulent viscosity used for diffusion only
    // internal faces
    for (int ifa = nfa_b; ifa < nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];

      double dx0[3], dx1[3];
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));
      double invwtot = 1.0/(w0+w1);

      double rho_fa    = (w1*rho[icv0] + w0*rho[icv1])*invwtot;
      double kine_fa   = (w1*kine[icv0] + w0*kine[icv1])*invwtot;
      double turbTS_fa = (w1*turbTS[icv0] + w0*turbTS[icv1])*invwtot;
      double cmus_fa   = (w1*cmus[icv0] + w0*cmus[icv1])*invwtot;

      mut_fa[ifa] = min(max(cmus_fa*rho_fa*kine_fa*turbTS_fa, 0.0), 1.0);
    }

    //boundary faces
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    {
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))                             // if wall ensure nu_t = 0.0
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            mut_fa[ifa] = 0.0;
        }
        else
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)     // otherwise make first order extrapolation to face
          {
            int icv0 = cvofa[ifa][0];
            mut_fa[ifa] = min(max(cmus[icv0]*rho[icv0]*kine[icv0]*turbTS[icv0], 0.0), 1.0);
          }
      }
    }

    for (int icv=0; icv<ncv; icv++)
      muT[icv] = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);
  }

  /**
   * calculate diffusivity scalars
   */
  virtual void diffusivityHookScalarRansTurb(const string &name)
  {
    ScalarTranspEq *scal;

    if (name == "kine")
    {
      scal = getScalarTransportData(name);
      for (int ifa = 0; ifa < nfa; ifa++)
        scal->diff[ifa] = mul_fa[ifa] + mut_fa[ifa]/scal->turbSchmidtNumber;
    }

    if (name == "omega")
    {
      scal = getScalarTransportData(name);
      for (int ifa = 0; ifa < nfa; ifa++)
        scal->diff[ifa] = mul_fa[ifa] + mut_fa[ifa]/scal->turbSchmidtNumber;
    }
  }

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")
    {
      for (int icv = 0; icv < ncv; icv++)
      {
        double src  = getTurbProd(icv) - fbstar[icv]*rho[icv]*kine[icv]*omega[icv];
        rhs[icv] += src*cv_volume[icv];

        //d(rho*kine*eps/kine)/d(rho*kine)
        if (flagImplicit)
        {
          double dsrcdphi = -fbstar[icv]*omega[icv];
          int noc00 = nbocv_i[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    }

    if (name == "omega")
      for (int icv = 0; icv < ncv; icv++)
      {
        double TS = sqrt(kine[icv]*kine[icv]/(omega[icv]*omega[icv]) + 1.e-12);

        double src = G_OM*getTurbProd(icv)/TS - BETA*rho[icv]*omega[icv]*omega[icv];
        rhs[icv]  += src*cv_volume[icv];

        // d(ce2*rho*eps/TS)/d(rho*eps)
        if (flagImplicit)
        {
          double dsrcdphi = -2.0*BETA*omega[icv];
          int noc00 = nbocv_i[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
  }

  virtual void sourceHookRansTurbCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
    int kine_Index, omega_Index;
    double dsrcdphi;

    kine_Index  = 5 + getScalarTransportIndex("kine");
    omega_Index = 5 + getScalarTransportIndex("omega");

    for (int icv = 0; icv < ncv; icv++)
    {
      double TS = sqrt(kine[icv]*kine[icv]/(omega[icv]*omega[icv]) + 1.e-12);

      double kine_src  = getTurbProd(icv) - fbstar[icv]*rho[icv]*kine[icv]*omega[icv];
      double omega_src = G_OM*getTurbProd(icv)/TS - BETA*rho[icv]*omega[icv]*omega[icv];

      rhs[icv][kine_Index]  += kine_src*cv_volume[icv];
      rhs[icv][omega_Index] += omega_src*cv_volume[icv];

      if (flagImplicit)
      {
        int noc00 = nbocv_i[icv];

        dsrcdphi = -fbstar[icv]*omega[icv];
        A[noc00][kine_Index][kine_Index] -= dsrcdphi*cv_volume[icv];
        dsrcdphi = -fbstar[icv]*kine[icv];
        A[noc00][kine_Index][omega_Index] -= dsrcdphi*cv_volume[icv];

        dsrcdphi = -BETA*omega[icv];
        A[noc00][omega_Index][omega_Index] -= dsrcdphi*cv_volume[icv];
      }
    }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_fa, FaZone *zone, const string &name)
  {
    if (name == "omega")
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int icv0 = cvofa[ifa][0];

        double nVec[3], s_half[3];
        normVec3d(nVec, fa_normal[ifa]);
        vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
        double wallDist = fabs(vecDotVec3d(s_half, nVec));
        double nuLamCV = calcMuLam(icv0)/rho[icv0];
        double omegaWall = 60.0*nuLamCV/(BETA*wallDist*wallDist);

        phi_fa[ifa] = omegaWall;
      }
    }
  }

  void calcTurbTimeScale()
  {
    double nu, eps, timeScale, realScale;
    for (int icv = 0; icv < ncv; icv++)
    {
      timeScale = 1.0/omega[icv];
      //timeScale = sqrt(timeScale*timeScale + 1.0e-12);//1.0e-6 limiter

      if (REALIZABLE)
      {
        realScale = 0.6/max(sqrt(3.0)*strMag[icv],1.0e-14);
        timeScale = min(timeScale, realScale);
      }
      turbTS[icv] = timeScale;
    }
    updateCvData(turbTS, REPLACE_DATA);
  }

  double getTurbProd(int icv)
  {
    double mu_t = min(max(rho[icv]*cmus[icv]*kine[icv]/omega[icv], 0.0), 1.0);
    double Pk = mu_t*strMag[icv]*strMag[icv] - 2.0/3.0*rho[icv]*kine[icv]*diverg[icv];

    if (LIMIT_PK == 1)
      Pk = min(Pk, 20.0*fbstar[icv]*rho[icv]*kine[icv]*omega[icv]);

    return Pk;
  }

  virtual void calcFbstarFunction()
  {
    calcCvScalarGrad(grad_kine,  kine,  kine_bfa,  gradreconstruction, limiterNavierS, kine,  epsilonSDWLS);
    calcCvScalarGrad(grad_omega, omega, omega_bfa, gradreconstruction, limiterNavierS, omega, epsilonSDWLS);

    for (int icv = 0; icv < ncv; icv++)
    {
      double chi = C_MU*C_MU/(omega[icv]*omega[icv]*omega[icv])*vecDotVec3d(grad_kine[icv], grad_omega[icv]);

      if (chi > 0.0)
        fbstar[icv] = (1.0 + 680.0*chi*chi)/(1.0 + 400*chi*chi);
      else
        fbstar[icv] = 1.0;
    }
  }
};

/*
 * Nonlinear EASM k-epsilon (2003) Model (EASMke2003)
 */

// differences with standard EASM:
// muT is being limited between zero and one.
// -the time scale, used for cmus and eps sources, has a Kolmogorov lower limit
// -the production of w (instead of just k) is also being limited.

class RansTurbEASMkeps : virtual public RansTurbEASM
{
public:
  RansTurbEASMkeps()
  {
    if (mpi_rank == 0)
      cout << "RansTurbEASMkeps()" << endl;

    turbModel = EASMkeps;

    C_MU  = getDoubleParam("C_MU", "0.0885");
    CEPS1 = getDoubleParam("CEPS1"  , "1.44");
    CEPS2 = getDoubleParam("CEPS2"  , "1.83");

    SIG_K = getDoubleParam("SIG_K"  , "1.0");
    SIG_D = 0.41*0.41/(sqrt(C_MU)*(CEPS2 - CEPS1));

    LIMIT_PK   = getIntParam("LIMIT_PK", "1");
    REALIZABLE = getIntParam("REALIZABLE", "0");

    ScalarTranspEq *eq;
    eq = registerScalarTransport("kine",  CV_DATA);
    eq->relax = getDoubleParam("RELAX_kine", "0.7");
    eq->phiZero = 1.0e-8;
    eq->phiMaxiter = 500;
    eq->lowerBound = 1.0e-6;
    eq->upperBound = 1.0e10;
    eq->turbSchmidtNumber = SIG_K;

    eq = registerScalarTransport("eps", CV_DATA);
    eq->relax = getDoubleParam("RELAX_eps", "0.7");
    eq->phiZero = 1.0e-8;
    eq->phiMaxiter = 500;
    eq->lowerBound = 1.0e-6;
    eq->upperBound = 1.0e10;
    eq->turbSchmidtNumber = SIG_D;

    strMag = NULL;       registerScalar(strMag, "strMag", CV_DATA);
    diverg = NULL;       registerScalar(diverg, "diverg", CV_DATA);
    muT    = NULL;       registerScalar(muT,    "muT",    CV_DATA);
    turbTS = NULL;       registerScalar(turbTS, "turbTS", CV_DATA);

    wallDist = NULL;       registerScalar(wallDist, "wallDist",  CV_DATA);
    wallConn = NULL;       // array of integers
  }

  virtual ~RansTurbEASMkeps() {}

public:

  double *eps;                ///< introduced to have access to variables
  double *kine_bfa, *eps_bfa; ///< turbulent scalars at the boundary
  double *muT;                ///< turbulent viscosity at cell center for output

  double C_MU, SIG_K, SIG_D, CEPS1, CEPS2; ///< model constants

  double *wallDist; ///< distance to closest wall face
  int *wallConn;    ///< index of closest wall face

  int LIMIT_PK;   ///< limiter for tke production
  int REALIZABLE; ///< realizability limits for turbTS, turbLS

public:

  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0)
      cout << "initialHookScalarRansTurbModel()" << endl;

    ScalarTranspEq *eq;
    eq = getScalarTransportData("kine");   kine = eq->phi;    kine_bfa = eq->phi_bfa;
    eq = getScalarTransportData("eps");    eps = eq->phi;     eps_bfa = eq->phi_bfa;

    wallConn = new int[ncv];
    calcWallDistance(wallConn, wallDist);

    for (int ifa = 0; ifa < nfa; ifa++)
      nonLinear[ifa] = 1.0;
    if (mpi_rank == 0)
      cout << "NON-LINEAR DOMAIN ACTIVATED." << endl;
  }

  virtual void calcRansTurbViscMuet()
  {
    // update velocity gradients, timescale
    calcGradVel();
    calcStrainRateAndDivergence();
    calcTurbTimeScale();

    // compute Reynolds stresses, cmus
    calcRsCenterEASM();

    // Calculate turbulent viscosity used for diffusion only
    // internal faces
    for (int ifa = nfa_b; ifa < nfa; ifa++)
    {
      int icv0 = cvofa[ifa][0];
      int icv1 = cvofa[ifa][1];

      double dx0[3], dx1[3];
      vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
      vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
      double w0 = sqrt(vecDotVec3d(dx0, dx0));
      double w1 = sqrt(vecDotVec3d(dx1, dx1));
      double invwtot = 1.0/(w0+w1);

      double rho_fa    = (w1*rho[icv0] + w0*rho[icv1])*invwtot;
      double kine_fa   = (w1*kine[icv0] + w0*kine[icv1])*invwtot;
      double turbTS_fa = (w1*turbTS[icv0] + w0*turbTS[icv1])*invwtot;

      mut_fa[ifa] = min(max(C_MU*rho_fa*kine_fa*turbTS_fa, 0.0), 1.0);
    }

    //boundary faces
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    {
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))                             // if wall ensure nu_t = 0.0
        {
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            mut_fa[ifa] = 0.0;
        }
        else
          for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)     // otherwise make first order extrapolation to face
          {
            int icv0 = cvofa[ifa][0];
            mut_fa[ifa] = min(max(C_MU*rho[icv0]*kine[icv0]*turbTS[icv0], 0.0), 1.0);
          }
      }
    }

    for (int icv=0; icv<ncv; icv++)
      muT[icv] = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);
  }

  /**
   * calculate diffusivity scalars
   */
  virtual void diffusivityHookScalarRansTurb(const string &name)
  {
    ScalarTranspEq *scal;

    if (name == "kine")
    {
      scal = getScalarTransportData(name);
      for (int ifa = 0; ifa < nfa; ifa++)
        scal->diff[ifa] = mul_fa[ifa] + mut_fa[ifa]/scal->turbSchmidtNumber;
    }

    if (name == "eps")
    {
      scal = getScalarTransportData(name);
      for (int ifa = 0; ifa < nfa; ifa++)
        scal->diff[ifa] = mul_fa[ifa] + mut_fa[ifa]/scal->turbSchmidtNumber;
    }
  }

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    if (name == "kine")
    {
      for (int icv = 0; icv < ncv; icv++)
      {
        double src  = getTurbProd(icv) - rho[icv]*eps[icv];
        rhs[icv] += src*cv_volume[icv];

        //d(rho*kine*eps/kine)/d(rho*kine)
        if (flagImplicit)
        {
          double dsrcdphi = -eps[icv]/kine[icv];
          int noc00 = nbocv_i[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    }

    if (name == "eps")
      for (int icv = 0; icv < ncv; icv++)
      {
        double TS = turbTS[icv];

        double muLamCV = calcMuLam(icv);
        double Rek = rho[icv]*sqrt(kine[icv])*wallDist[icv]/muLamCV;
        double fe = 1.0 - exp(-Rek/10.8);

        double src = (CEPS1*getTurbProd(icv) - CEPS2*fe*rho[icv]*eps[icv])/TS;
        rhs[icv]  += src*cv_volume[icv];

        // d(ce2*rho*eps/TS)/d(rho*eps)
        if (flagImplicit)
        {
          double dsrcdphi = -CEPS2*fe/TS;
          int noc00 = nbocv_i[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
  }

  virtual void sourceHookRansTurbCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
    int kine_Index, omega_Index, eps_Index;
    double dsrcdphi;

    kine_Index = 5 + getScalarTransportIndex("kine");
    eps_Index  = 5 + getScalarTransportIndex("eps");

    for (int icv = 0; icv < ncv; icv++)
    {
      double TS = turbTS[icv];

      double muLamCV = calcMuLam(icv);
      double Rek = rho[icv]*sqrt(kine[icv])*wallDist[icv]/muLamCV;
      double fe = 1.0 - exp(-Rek/10.8);

      double kine_src = getTurbProd(icv) - rho[icv]*eps[icv];
      double eps_src  = (CEPS1*getTurbProd(icv) - CEPS2*fe*rho[icv]*eps[icv])/TS;

      rhs[icv][kine_Index] += kine_src*cv_volume[icv];
      rhs[icv][eps_Index]  += eps_src*cv_volume[icv];

      if (flagImplicit)
      {
        int noc00 = nbocv_i[icv];

        dsrcdphi = -eps[icv]*kine[icv];
        A[noc00][kine_Index][kine_Index] -= dsrcdphi*cv_volume[icv];
        dsrcdphi = -1.0;
        A[noc00][kine_Index][eps_Index] -= dsrcdphi*cv_volume[icv];

        dsrcdphi = -CEPS2*fe/TS;
        A[noc00][omega_Index][eps_Index] -= dsrcdphi*cv_volume[icv];
      }
    }
  }

  virtual void boundaryHookScalarRansTurb(double *phi_fa, FaZone *zone, const string &name)
  {
    if (name == "eps")
    {
      for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      {
        int icv0 = cvofa[ifa][0];

        double nVec[3], s_half[3];
        normVec3d(nVec, fa_normal[ifa]);
        vecMinVec3d(s_half, x_fa[ifa], x_cv[icv0]);
        double wallDist = fabs(vecDotVec3d(s_half, nVec));
        double nuLamCV = calcMuLam(icv0)/rho[icv0];
        double epsWall = 2.0*nuLamCV*kine[icv0]/(wallDist*wallDist);

        phi_fa[ifa] = epsWall;
      }
    }
  }

  void calcTurbTimeScale()
  {
    double nu, timeScale, realScale;
    for (int icv = 0; icv < ncv; icv++)
    {
      nu = calcMuLam(icv)/rho[icv];
      timeScale = max(kine[icv]/eps[icv], 6.0*sqrt(nu/eps[icv])); //Kolmogorov limiter

      if (REALIZABLE)
      {
        realScale = 0.6/max(sqrt(3.0)*strMag[icv],1.0e-14);
        timeScale = min(timeScale, realScale);
      }
      turbTS[icv] = timeScale;
    }
    updateCvData(turbTS, REPLACE_DATA);
  }

  double getTurbProd(int icv)
  {
    double mu_t = min(max(rho[icv]*C_MU*kine[icv]*kine[icv]/eps[icv], 0.0), 1.0);
    double Pk = mu_t*strMag[icv]*strMag[icv] - 2.0/3.0*rho[icv]*kine[icv]*diverg[icv];

    if (LIMIT_PK == 1)
      Pk = min(Pk, 20.0*rho[icv]*eps[icv]);

    return Pk;
  }

  virtual void finalHookScalarRansTurbModel()
  {
    if (wallConn != NULL) {delete [] wallConn; wallConn = NULL;}
  }
};



#endif
