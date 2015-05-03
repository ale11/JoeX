#ifndef RANSTURBMODEL_EASM_H
#define RANSTURBMODEL_EASM_H

#include "UgpWithCvCompFlow.h"

class RansTurbEASM : virtual public UgpWithCvCompFlow
{
public: // constructor, destructor
  RansTurbEASM()
  {
    cmus   = NULL;       registerScalar(cmus,   "cmus",   CV_DATA);

    rij_diag_nd   = NULL;     registerVector(rij_diag_nd,    "rij_diag_nd",  CV_DATA);
    rij_offdiag_nd = NULL;    registerVector(rij_offdiag_nd, "rij_offdiag_nd", CV_DATA);

    ver2 = getIntParam("VERSION2", "0");

    debug1 = NULL;     registerVector(debug1,         "debug1",         CV_DATA);
    debug2 = NULL;     registerVector(debug2,         "debug2",         CV_DATA);
    debug3 = NULL;     registerVector(debug3,         "debug3",         CV_DATA);
  }

  virtual ~RansTurbEASM() {}

public: // member variables
  double *cmus; ///< non-constant muT coefficient

  double (*rij_diag_nd)[3];
  double (*rij_offdiag_nd)[3];

  double (*debug1)[3];
  double (*debug2)[3];
  double (*debug3)[3];

  int ver2;

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

    double c2r = 1.2;
    double cs = 0.84;
    double ca = 4.0;
    double cb = -10.0;

    double gamma0 = 0.5*c11;
    double gamma1 = 0.5*c10 + (ceps2 - ceps1)/(ceps1 - 1.0);

    double a1     = 0.5*(4.0/3.0 - c2);
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

      double f = 0.5*(1.0 + tanh(0.5*ca*strMag[icv]*tau + cb));
      if (ver2) a1 = 2.0/3.0 - 0.5*c2r + 0.5*(c2r - c2)*f;

      p = -gamma1/(e2t2G0);
      q = (gamma1*gamma1 - 2.0*e2t2G0*a1 - 2.0/3.0*e2t2*a3*a3 + 2.0*R2e2t2*a2*a2);
      q /= 4.0*e2t2G0*e2t2G0;
      r = gamma1*a1/(4.0*e2t2G0*e2t2G0);

      if (ver2)
      {
      	p += cs*(1.0 - f)/(2.0*gamma0);
      	q -= gamma1*(1.0 - f)*cs*e2t2/(4.0*e2t2G0*e2t2G0);
      }

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

      rij_diag[icv][0] = -max(r11,0.)*2.0*kine[icv]*rho[icv];
      rij_diag[icv][1] = -max(r22,0.)*2.0*kine[icv]*rho[icv];
      rij_diag[icv][2] = -max(r33,0.)*2.0*kine[icv]*rho[icv];

      rij_offdiag[icv][0] = -r12*2.0*kine[icv]*rho[icv];
      rij_offdiag[icv][1] = -r13*2.0*kine[icv]*rho[icv];
      rij_offdiag[icv][2] = -r23*2.0*kine[icv]*rho[icv];

      rij_diag_nd[icv][0] = max(r11,0.);
      rij_diag_nd[icv][1] = max(r22,0.);
      rij_diag_nd[icv][2] = max(r33,0.);

      rij_offdiag_nd[icv][0] = r12;
      rij_offdiag_nd[icv][1] = r13;
      rij_offdiag_nd[icv][2] = r23;

      double b1kSk1 = (r11 - 1.0/3.0)*S11 + r12*S21 + r13*S31;
      double b2kSk2 = r12*S12 + (r22 - 1.0/3.0)*S22 + r23*S32;
      double b3kSk3 = r13*S13 + r23*S23 + (r33 - 1.0/3.0)*S33;
      double b1kSk2 = (r11 - 1.0/3.0)*S12 + r12*S22 + r13*S32;
      double b1kSk3 = (r11 - 1.0/3.0)*S13 + r12*S23 + r13*S33;
      double b2kSk3 = r12*S13 + (r22 - 1.0/3.0)*S23 + r23*S33;
      double b2kSk1 = r12*S11 + (r22 - 1.0/3.0)*S21 + r23*S31;
      double b3kSk1 = r13*S11 + r23*S21 + (r33 - 1.0/3.0)*S31;
      double b3kSk2 = r13*S12 + r23*S22 + (r33 - 1.0/3.0)*S32;

      double bnmSnm = b1kSk1 + b2kSk2 + b3kSk3;

      double b1kWk1 = (r11 - 1.0/3.0)*W11 + r12*W21 + r13*W31;
      double b2kWk2 = r12*W12 + (r22 - 1.0/3.0)*W22 + r23*W32;
      double b3kWk3 = r13*W13 + r23*W23 + (r33 - 1.0/3.0)*W33;
      double b1kWk2 = (r11 - 1.0/3.0)*W12 + r12*W22 + r13*W32;
      double b1kWk3 = (r11 - 1.0/3.0)*W13 + r12*W23 + r13*W33;
      double b2kWk3 = r12*W13 + (r22 - 1.0/3.0)*W23 + r23*W33;
      double b2kWk1 = r12*W11 + (r22 - 1.0/3.0)*W21 + r23*W31;
      double b3kWk1 = r13*W11 + r23*W21 + (r33 - 1.0/3.0)*W31;
      double b3kWk2 = r13*W12 + r23*W22 + (r33 - 1.0/3.0)*W32;

      double Peps = -2.0*tau*(r11*S11 + r22*S22 + r33*S33 + 2.0*(r12*S12 + r13*S13 + r23*S23));
      double c2star = c2r + (c2 - c2r)*f + cs*Peps*(f - 1.0);

      //debug1[icv][0] = Peps;
      //debug1[icv][1] = f;
      //debug1[icv][2] = c2star;

      if (ver2) c2 = c2star;

      //debug2[icv][0] = -(c10 + c11*Peps)*(r11 - 1.0/3.0) + c2*S11*tau
      //		             +c3*(b1kSk1 + b1kSk1 - 2.0/3.0*bnmSnm)*tau - c4*(b1kWk1 + b1kWk1)*tau;
      //debug2[icv][1] = -(c10 + c11*Peps)*(r22 - 1.0/3.0) + c2*S22*tau
      //		             +c3*(b2kSk2 + b2kSk2 - 2.0/3.0*bnmSnm)*tau - c4*(b2kWk2 + b2kWk2)*tau;
      //debug2[icv][2] = -(c10 + c11*Peps)*(r33 - 1.0/3.0) + c2*S33*tau
      //		             +c3*(b3kSk3 + b3kSk3 - 2.0/3.0*bnmSnm)*tau - c4*(b3kWk3 + b3kWk3)*tau;

      //debug3[icv][0] = -(c10 + c11*Peps)*r12 + c2*S12*tau
      //		             +c3*(b1kSk2 + b2kSk1)*tau - c4*(b1kWk2 + b2kWk1)*tau;
      //debug3[icv][1] = -(c10 + c11*Peps)*r13 + c2*S13*tau
      //		             +c3*(b1kSk3 + b3kSk1)*tau - c4*(b1kWk3 + b3kWk1)*tau;
      //debug3[icv][2] = -(c10 + c11*Peps)*r23 + c2*S23*tau
      //		             +c3*(b2kSk3 + b3kSk2)*tau - c4*(b2kWk3 + b3kWk2)*tau;
    }

    updateCvData(cmus, REPLACE_DATA);
    updateCvData(rij_diag, REPLACE_ROTATE_DATA);
    updateCvData(rij_offdiag, REPLACE_ROTATE_DATA);

    // Interpolate cell centered stresses to cell faces
    interpolateReStressToFace();
  }
};

//######################################################//
//                                                      //
// Nonlinear EASM k-omega Model (EASMko2003)            //
//                                                      //
//######################################################//

// differences with standard EASM:
// -muT is being limited between zero and one.
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

    if (ver2) SIG_O = 1.0/0.65;

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
  int REALIZABLE; ///< realizability limits for turbTS

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
    	calcCvScalarGrad(grad_kine,  kine,  kine_bfa,  gradreconstruction, limiterNavierS, kine,  epsilonSDWLS);
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

        debug1[icv][0] = getTurbProd(icv)/rho[icv]; //*cv_volume[icv];
        debug1[icv][1] = -fbstar[icv]*kine[icv]*omega[icv]; //*cv_volume[icv];
        debug1[icv][2] = vel[icv][0]*grad_kine[icv][0] + vel[icv][1]*grad_kine[icv][1] + vel[icv][2]*grad_kine[icv][2];
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

//######################################################//
//                                                      //
// Nonlinear EASM k-epsilon Model (EASMko2003)          //
//                                                      //
//######################################################//

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
