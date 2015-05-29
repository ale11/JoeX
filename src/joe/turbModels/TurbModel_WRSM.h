#ifndef RANSTURBMODEL_WRSM_H
#define RANSTURBMODEL_WRSM_H

#include "UgpWithCvCompFlow.h"

//######################################################//
//                                                      //
// Wilcox Stress-omega Full Reynolds Stress Model       //
//                                                      //
//######################################################//

class RansTurbWRSM : virtual public UgpWithCvCompFlow
{
public:   // constructors

  RansTurbWRSM()
  {
    if (mpi_rank == 0)
      cout << "RansTurbWRSM()" << endl;

    turbModel = WRSM;

    C1      = getDoubleParam("C1",      "1.8"  );
    C2      = getDoubleParam("C2",      "0.52632");
    alpha   = getDoubleParam("alpha",   "0.52");
    beta_s  = getDoubleParam("beta_s",  "0.09"     );
    sigma   = getDoubleParam("sigma",   "0.5"      );
    sigma_s = getDoubleParam("sigma_s", "0.6"      );
    beta_0  = getDoubleParam("beta_0",  "0.0708"   );

    alpha_h = (8.0 + C2)/11.0;
    beta_h  = (8.0*C2 - 2.0)/11.0;
    gamma_h = (60.0*C2 - 4.0)/55.0;

    ScalarTranspEq *eq;

    eq = registerScalarTransport("R11", CV_DATA);
    eq->relax = getDoubleParam("RELAX_R11", "0.4");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0e10;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    eq = registerScalarTransport("R22", CV_DATA);
    eq->relax = getDoubleParam("RELAX_R22", "0.4");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0e10;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    eq = registerScalarTransport("R33", CV_DATA);
    eq->relax = getDoubleParam("RELAX_R33", "0.4");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 1000;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0e10;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    eq = registerScalarTransport("R12", CV_DATA);
    eq->relax = getDoubleParam("RELAX_R12", "0.4");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 1000;
    eq->lowerBound = -1.0e10;
    eq->upperBound = 1.0e10;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    eq = registerScalarTransport("R13", CV_DATA);
    eq->relax = getDoubleParam("RELAX_R13", "0.4");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 1000;
    eq->lowerBound = -1.0e10;
    eq->upperBound = 1.0e10;
    eq->reconstruction = getStringParam("SCALAR_TURB_RECONSTRUCTION", "STANDARD");

    eq = registerScalarTransport("R23", CV_DATA);
    eq->relax = getDoubleParam("RELAX_R23", "0.4");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 1000;
    eq->lowerBound = -1.0e10;
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

    kine     = NULL;       registerScalar(kine,   "kine",   CV_DATA);
    strMag   = NULL;       registerScalar(strMag, "strMag", CV_DATA);
    diverg   = NULL;       registerScalar(diverg, "diverg", CV_DATA);
    muT      = NULL;       registerScalar(muT,    "muT",    CV_DATA);
    debug1   = NULL;       registerScalar(debug1, "debug1", CV_DATA);
    debug2   = NULL;       registerScalar(debug2, "debug2", CV_DATA);
    debug3   = NULL;       registerScalar(debug3, "debug3", CV_DATA);
  }

  virtual ~RansTurbWRSM() {}

public:

  // model variables
  double *R11, *R22, *R33; ///< normal Reynolds stresses
  double *R12, *R13, *R23; ///< shear Reynolds stresses
  double *omega;           ///< specific dissipation
  double (*grad_R11)[3];   ///< gradient of R11
  double (*grad_R22)[3];   ///< gradient of R22
  double (*grad_R33)[3];   ///< gradient of R33
  double (*grad_R12)[3];   ///< gradient of R12
  double (*grad_R13)[3];   ///< gradient of R13
  double (*grad_R23)[3];   ///< gradient of R23
  double (*grad_omega)[3]; ///< gradient of omega
  double *R11_bfa;         ///< R11 at the boundaries
  double *R22_bfa;         ///< R22 at the boundaries
  double *R33_bfa;         ///< R33 at the boundaries
  double *R12_bfa;         ///< R12 at the boundaries
  double *R13_bfa;         ///< R13 at the boundaries
  double *R23_bfa;         ///< R23 at the boundaries
  double *omega_bfa;       ///< omega at the boundaries
  double *muT;             ///< turbulent viscosity at cell center

  double *debug1, *debug2, *debug3;

  // model constants
  double C1, C2, alpha, beta_s, sigma, sigma_s, beta_0, alpha_h, beta_h, gamma_h;

public:

  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0)
      cout << "initialHook WILCOX Stress-omega model" << endl;

    // connect pointers
    ScalarTranspEq *eq;
    eq = getScalarTransportData("R11");     R11 = eq->phi;      R11_bfa = eq->phi_bfa;    grad_R11 = eq->grad_phi;
    eq = getScalarTransportData("R22");     R22 = eq->phi;      R22_bfa = eq->phi_bfa;    grad_R22 = eq->grad_phi;
    eq = getScalarTransportData("R33");     R33 = eq->phi;      R33_bfa = eq->phi_bfa;    grad_R33 = eq->grad_phi;
    eq = getScalarTransportData("R12");     R12 = eq->phi;      R12_bfa = eq->phi_bfa;    grad_R12 = eq->grad_phi;
    eq = getScalarTransportData("R13");     R13 = eq->phi;      R13_bfa = eq->phi_bfa;    grad_R13 = eq->grad_phi;
    eq = getScalarTransportData("R23");     R23 = eq->phi;      R23_bfa = eq->phi_bfa;    grad_R23 = eq->grad_phi;
    eq = getScalarTransportData("omega");   omega = eq->phi;    omega_bfa = eq->phi_bfa;  grad_omega = eq->grad_phi;

    for (int ifa = 0; ifa < nfa; ifa++)
      nonLinear[ifa] = 1.0;
    if (mpi_rank == 0)
      cout << "NON-LINEAR DOMAIN ACTIVATED." << endl;
  }

  virtual void calcRansTurbViscMuet()
  {
    // update velocity gradients
    calcGradVel();
    calcStrainRateAndDivergence();

    for (int icv = 0; icv < ncv; icv++)
    {
      // Enforce realizability
    	if (R12[icv]*R12[icv] > R11[icv]*R22[icv])
    		R12[icv] = sign(R12[icv])*sqrt(R11[icv]*R22[icv]);

    	if (R13[icv]*R13[icv] > R11[icv]*R33[icv])
    		R13[icv] = sign(R13[icv])*sqrt(R11[icv]*R33[icv]);

    	if (R23[icv]*R23[icv] > R22[icv]*R33[icv])
    		R23[icv] = sign(R23[icv])*sqrt(R22[icv]*R33[icv]);

    	// Store stresses
    	rij_diag[icv][0] = -rho[icv]*R11[icv];
    	rij_diag[icv][1] = -rho[icv]*R22[icv];
    	rij_diag[icv][2] = -rho[icv]*R33[icv];

    	rij_offdiag[icv][0] = -rho[icv]*R12[icv];
    	rij_offdiag[icv][1] = -rho[icv]*R13[icv];
    	rij_offdiag[icv][2] = -rho[icv]*R23[icv];

    	// TKE and eddy viscosity
    	kine[icv] = 0.5*(R11[icv] + R22[icv] + R33[icv]);
      muT[icv] = min(rho[icv]*kine[icv]/omega[icv], 1000.0);
    }

    updateCvData(R11, REPLACE_DATA);
    updateCvData(R22, REPLACE_DATA);
    updateCvData(R33, REPLACE_DATA);
    updateCvData(R12, REPLACE_DATA);
    updateCvData(R13, REPLACE_DATA);
    updateCvData(R23, REPLACE_DATA);
    updateCvData(rij_diag, REPLACE_ROTATE_DATA);
    updateCvData(rij_offdiag, REPLACE_ROTATE_DATA);
    updateCvData(kine, REPLACE_DATA);
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
            mut_fa[ifa] = 0.0;        // set mut zero at walls
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            mut_fa[ifa] = muT[icv0];  // zero order extrapolation for others
          }
      }

    // Interpolate cell-centered stresses to faces
    interpolateReStressToFace();
  }

  virtual void diffusivityHookScalarRansTurb(const string &name)
  {
    ScalarTranspEq *eq;

    if ((name == "R11") || (name == "R22") || (name == "R33") ||
    		(name == "R12") || (name == "R13") || (name == "R23") || (name == "omega"))
    {
      double sigma_muT;
      if (name != "omega") sigma_muT = sigma_s;
      else                 sigma_muT = sigma;

      eq = getScalarTransportData(name);

      for (int ifa = 0; ifa < nfa; ifa++)
        eq->diff[ifa] = mul_fa[ifa] + sigma_muT*mut_fa[ifa];
    }
  }

  virtual void sourceHookScalarRansTurb_new(double *rhs, double *A, const string &name, int flagImplicit)
  {
    double S[3][3], W[3][3], b[3][3];

    if (name != "omega")
    {
    	int i,j;
    	if      (name == "R11") {i = 0; j = 0;}
    	else if (name == "R22") {i = 1; j = 1;}
    	else if (name == "R33") {i = 2; j = 2;}
    	else if (name == "R12") {i = 0; j = 1;}
    	else if (name == "R13") {i = 0; j = 2;}
    	else if (name == "R23")	{i = 1; j = 2;}
    	else{
    		if (mpi_rank == 0)
    			cout << "WARNING: wrong name for sourceHookScalarRansTurb_new!" << endl;
    	}

      double p3 = 2.0/3.0*(alpha_h + beta_h) - 0.5*gamma_h;
      double p4 = alpha_h + beta_h;
      double p5 = alpha_h - beta_h;
      double delta[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

      for (int icv = 0; icv < ncv; icv++)
      {
        calcStrainTensor(icv,S);
        calcRotationTensor(icv,W);
        calcAnisotropyTensor(icv,b);

        // Integrity basis
        double trace_Sb = 0.0;
        for (int m = 0; m < 3; m++)
        	for (int n = 0; n < 3; n++)
        		trace_Sb += S[m][n]*b[n][m];

        double h4 = 0.0;
        for (int k = 0; k < 3; k++)
        	h4 += S[i][k]*b[k][j] + b[i][k]*S[k][j];
        h4 -= 2.0/3.0*trace_Sb*delta[i][j];

        double h5 = 0.0;
        for (int k = 0; k < 3; k++)
        	h5 += W[i][k]*b[k][j] - b[i][k]*W[k][j];

      	// Sources
        double prod = -2.0*kine[icv]*(h4 + 2.0/3.0*trace_Sb*delta[i][j] + h5 + 2.0/3.0*S[i][j]);

        double eps = beta_s*kine[icv]*omega[icv];

        double redi = - 2.0*C1*eps*b[i][j] + 2.0*kine[icv]*(p4*h4 + p5*h5
                      + p3*(S[i][j] - 1.0/3.0*diverg[icv]*delta[i][j]));

        // Total contribution
        double src = rho[icv]*(prod + redi - 2.0/3.0*eps*delta[i][j]);
        rhs[icv] += src*cv_volume[icv];

        if (flagImplicit)
        {
          int noc00 = nbocv_i[icv];
          double dproddphi = -(S[i][i] + S[j][j]);
          double depsdphi  = -1.0/3.0*beta_s*omega[icv]*delta[i][j];
          double dredidphi = -C1*beta_s*omega[icv]*(1.0 - 1.0/3.0*delta[i][j])
          		               +p4*(S[i][i] + S[j][j] - 2.0/3.0*(2.0*S[i][j] - 1.0/3.0*diverg[icv])*delta[i][j])
          		               +p3*(S[i][j] - 1.0/3.0*diverg[icv])*delta[i][j];
          double dsrcdphi = min(dproddphi + depsdphi + dredidphi, 0.0);
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
    }

    if (name == "omega")
    {
    	calcCvScalarGrad(grad_R11, R11, R11_bfa, gradreconstruction, limiterNavierS, R11, epsilonSDWLS);
    	calcCvScalarGrad(grad_R22, R22, R22_bfa, gradreconstruction, limiterNavierS, R22, epsilonSDWLS);
    	calcCvScalarGrad(grad_R33, R33, R33_bfa, gradreconstruction, limiterNavierS, R33, epsilonSDWLS);
      calcCvScalarGrad(grad_omega, omega, omega_bfa, gradreconstruction, limiterNavierS, omega, epsilonSDWLS);

      for (int icv = 0; icv < ncv; icv++)
      {
      	calcStrainTensor(icv,S);
        calcRotationTensor(icv,W);
        calcAnisotropyTensor(icv,b);

        // Production of omega
      	double trace_Sb = 0.0;
        for (int m = 0; m < 3; m++)
        	for (int n = 0; n < 3; n++)
        		trace_Sb += S[m][n]*b[n][m];

        double prod = -alpha*omega[icv]*2.0*(trace_Sb + 1.0/3.0*diverg[icv]);
        prod = max(prod, 0.0);

        // Destruction of omega
        double S_hat[3][3];
        for (int i=0; i<3; i++)
          for (int j=0; j<3; j++)
          {
            // S_hat subtracts 0.5 divergence instead of 1/3!!!
            if (i==j)   S_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i] - diverg[icv]);
            else        S_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]);
          }

        double chiOm = 0.0;
        for (int i=0; i<3; i++)
          for (int j=0; j<3; j++)
            for (int k=0; k<3; k++)
              chiOm += W[i][j]*W[j][k]*S_hat[k][i];
        chiOm = fabs(chiOm/pow(beta_s*omega[icv], 3.0));

        double fbeta = (1.0 + 85.0*chiOm)/(1.0 + 100.0*chiOm);
        double beta  = beta_0*fbeta;
        double diss  = beta*omega[icv]*omega[icv];

        // Cross-diffusion of omega
      	double grad_kine[3];
        for (int i = 0; i < 3; i++)
        	grad_kine[i] = 0.5*(grad_R11[icv][i] + grad_R22[icv][i] + grad_R33[icv][i]);

        double crossDiff = vecDotVec3d(grad_kine, grad_omega[icv]);
        double sigmad;
        if (crossDiff <= 0.0) sigmad = 0.0;
        else                  sigmad = 0.125;
        crossDiff *= sigmad/omega[icv];

        // Total contribution
        double src =  rho[icv]*(prod - diss + crossDiff);
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
    int R11_Index = getScalarTransportIndex("R11");
    int R22_Index = getScalarTransportIndex("R22");
    int R33_Index = getScalarTransportIndex("R33");
    int R12_Index = getScalarTransportIndex("R12");
    int R13_Index = getScalarTransportIndex("R13");
    int R23_Index = getScalarTransportIndex("R23");
    int omega_Index = getScalarTransportIndex("omega");

    double S[3][3], W[3][3], b[3][3];

    double p3 = 2.0/3.0*(alpha_h + beta_h) - 0.5*gamma_h;
    double p4 = alpha_h + beta_h;
    double p5 = alpha_h - beta_h;
    double delta[3][3] = {{1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};

  	calcCvScalarGrad(grad_R11, R11, R11_bfa, gradreconstruction, limiterNavierS, R11, epsilonSDWLS);
  	calcCvScalarGrad(grad_R22, R22, R22_bfa, gradreconstruction, limiterNavierS, R22, epsilonSDWLS);
  	calcCvScalarGrad(grad_R33, R33, R33_bfa, gradreconstruction, limiterNavierS, R33, epsilonSDWLS);
    calcCvScalarGrad(grad_omega, omega, omega_bfa, gradreconstruction, limiterNavierS, omega, epsilonSDWLS);

    for (int icv = 0; icv < ncv; icv++)
    {
    	calcStrainTensor(icv,S);
      calcRotationTensor(icv,W);
      calcAnisotropyTensor(icv,b);

      double trace_Sb = 0.0;
      for (int m = 0; m < 3; m++)
      	for (int n = 0; n < 3; n++)
      		trace_Sb += S[m][n]*b[n][m];

      double h4[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
      double h5[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
      double srcRij[3][3] = { {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };

      // REYNOLDS STRESSES
      for (int i = 0; i < 3; i++)
      	for (int j = i; j < 3; j++)
      	{
      		// Integrity basis
      		for (int k = 0; k < 3; k++)
      			h4[i][j] += S[i][k]*b[k][j] + b[i][k]*S[k][j];
      		h4[i][j] -= 2.0/3.0*trace_Sb*delta[i][j];

      		for (int k = 0; k < 3; k++)
      			h5[i][j] += W[i][k]*b[k][j] - b[i][k]*W[k][j];

          // Sources
          double prod = -2.0*kine[icv]*(h4[i][j] + 2.0/3.0*trace_Sb*delta[i][j] + h5[i][j] + 2.0/3.0*S[i][j]);

          double eps = beta_s*kine[icv]*omega[icv];

          double redi = - 2.0*C1*eps*b[i][j] + 2.0*kine[icv]*(p4*h4[i][j] + p5*h5[i][j]
                        + p3*(S[i][j] - 1.0/3.0*diverg[icv]*delta[i][j]));

          // Total contribution
          srcRij[i][j] = rho[icv]*(prod + redi - 2.0/3.0*eps*delta[i][j]);
      	}

      // OMEGA
      // Production
      double prod = -alpha*omega[icv]*2.0*(trace_Sb + 1.0/3.0*diverg[icv]);
      prod = max(prod, 0.0);

      // Destruction
      double S_hat[3][3];
      for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        {

          // S_hat subtracts 0.5 divergence instead 1/3, no bug!!!
          if (i==j)   S_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i] - diverg[icv]);
          else        S_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]);
        }

      double chiOm = 0.0;
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
          for (int k=0; k<3; k++)
            chiOm += W[i][j]*W[j][k]*S_hat[k][i];
      chiOm = fabs(chiOm/pow(beta_s*omega[icv], 3.0));

      double fbeta = (1.0 + 85.0*chiOm)/(1.0 + 100.0*chiOm);
      double beta  = beta_0*fbeta;
      double diss  = beta*omega[icv]*omega[icv];

      // Cross-diffusion
    	double grad_kine[3];
      for (int i = 0; i < 3; i++)
      	grad_kine[i] = 0.5*(grad_R11[icv][i] + grad_R22[icv][i] + grad_R33[icv][i]);

      double crossDiff = vecDotVec3d(grad_kine, grad_omega[icv]);
      double sigmad;
      if (crossDiff <= 0.0) sigmad = 0.0;
      else                  sigmad = 0.125;
      crossDiff *= sigmad/omega[icv];

      double srcOm =  rho[icv]*(prod - diss + crossDiff);

      // RERSIDUALS AND JACOBIANS
      rhs[icv][5+R11_Index]   += srcRij[0][0]*cv_volume[icv];
      rhs[icv][5+R22_Index]   += srcRij[1][1]*cv_volume[icv];
      rhs[icv][5+R33_Index]   += srcRij[2][2]*cv_volume[icv];
      rhs[icv][5+R12_Index]   += srcRij[0][1]*cv_volume[icv];
      rhs[icv][5+R13_Index]   += srcRij[0][2]*cv_volume[icv];
      rhs[icv][5+R23_Index]   += srcRij[1][2]*cv_volume[icv];
      rhs[icv][5+omega_Index] += srcOm*cv_volume[icv];

      if (flagImplicit)
      {
        int noc00 = nbocv_i[icv];
        A[noc00][5+R11_Index][5+R11_Index] -= -(1.0/3.0*beta_s*omega[icv] + 2.0/3.0*C1*beta_s*omega[icv])*cv_volume[icv];
        A[noc00][5+R22_Index][5+R22_Index] -= -(1.0/3.0*beta_s*omega[icv] + 2.0/3.0*C1*beta_s*omega[icv])*cv_volume[icv];
        A[noc00][5+R33_Index][5+R33_Index] -= -(1.0/3.0*beta_s*omega[icv] + 2.0/3.0*C1*beta_s*omega[icv])*cv_volume[icv];
        A[noc00][5+R12_Index][5+R12_Index] -= -C1*beta_s*omega[icv]*cv_volume[icv];
        A[noc00][5+R13_Index][5+R13_Index] -= -C1*beta_s*omega[icv]*cv_volume[icv];
        A[noc00][5+R23_Index][5+R23_Index] -= -C1*beta_s*omega[icv]*cv_volume[icv];
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

            double nVec[3], s_half[3];
            normVec3d(nVec, fa_normal[ifa]);
            vecMinVec3d(s_half, x_fa[ifa], x_cv[icv]);
            double wallDist = fabs(vecDotVec3d(s_half, nVec));

            double nuLamCV = calcMuLam(icv)/rho[icv];
            phi_fa[ifa] = 6.0*nuLamCV/(beta_0*wallDist*wallDist);
          }
        }
    }
  }

  inline void calcStrainTensor(int icv, double (*S)[3])
  {
    S[0][0] = grad_u[icv][0][0];
    S[0][1] = 0.5*(grad_u[icv][0][1] + grad_u[icv][1][0]);
    S[0][2] = 0.5*(grad_u[icv][0][2] + grad_u[icv][2][0]);

    S[1][0] = S[0][1];
    S[1][1] = grad_u[icv][1][1];
    S[1][2] = 0.5*(grad_u[icv][1][2] + grad_u[icv][2][1]);

    S[2][0] = S[0][2];
    S[2][1] = S[1][2];
    S[2][2] = grad_u[icv][2][2];
  }

  inline void calcRotationTensor(int icv, double (*W)[3])
  {
    W[0][0] = 0.0;
    W[0][1] = 0.5*(grad_u[icv][0][1] - grad_u[icv][1][0]);
    W[0][2] = 0.5*(grad_u[icv][0][2] - grad_u[icv][2][0]);

    W[1][0] = -W[0][1];
    W[1][1] = 0.0;
    W[1][2] = 0.5*(grad_u[icv][1][2] - grad_u[icv][2][1]);

    W[2][0] = -W[0][2];
    W[2][1] = -W[1][2];
    W[2][2] = 0.0;
  }


  inline void calcAnisotropyTensor(int icv, double (*b)[3])
  {
    double q2 = 2.0*kine[icv];
    b[0][0] = R11[icv]/q2 - 1.0/3.0;
    b[0][1] = R12[icv]/q2;
    b[0][2] = R13[icv]/q2;

    b[1][0] = b[0][1];
    b[1][1] = R22[icv]/q2 - 1.0/3.0;
    b[1][2] = R23[icv]/q2;

    b[2][0] = b[0][2];
    b[2][1] = b[1][2];
    b[2][2] = R33[icv]/q2 - 1.0/3.0;
  }

  virtual void finalHookScalarRansTurbModel()
  {

  }

};



#endif
