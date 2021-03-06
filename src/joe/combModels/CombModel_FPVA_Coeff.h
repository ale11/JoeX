/*
 * CombModel_FPVA_new.h
 *
 *  Created on: February 15, 2010
 *      Author: Vincent Terrapon
 *      Update: 
 */

#ifndef COMBMODEL_FPVA_COEFF_H_
#define COMBMODEL_FPVA_COEFF_H_

#define PRESSURE_SCALING_FPVA_COEFF  1      //  Correct value as per Vincent's analysis -- all calculations prior to March 2011 used (erroneously) the value 2
#define MIXING_VISCOSITY

#include "CombModel_Base.h"
#include "ChemistryTableCartesianLinear.h"
#include "ChemistryTableCartesianCubic.h"
#include "ChemistryTableAdaptiveLinear.h"

#define Tminimum 20.0


// #############################################################################################
// ------                                                                                 ------
// ----                            RansCombFPVA_new                                         ----
// ------                                                                                 ------
// #############################################################################################

/*! \brief Combustion model based on Flamelet Progress Variable Approach (FPVA). 
 *
 *  The coefficients R0, T0, E0, ... are functions of the mean scalar mixture fraction \p Zmean,
 *  the variance of the scalar mixture fraction \p Zvar and the progress variable \p Cmean.
 */
template <class Chemtable>
class RansCombFPVA_Coeff: public RansCombBase {
public:
  // member variables

  Chemtable          myChemTable;                  ///< Chemistry table based on flamelets.
  double             *ZMean;                       ///< Mixture fraction Z (mean).
  double             *ZVar;                        ///< Variance of mixture fraction.
  double             *CMean;                       ///< Progress variable C (mean).
//  double                   *rhoDkZ;                      ///< Diffusion coefficient of mixture fraction Z.
  double             *muLam;                       ///< Laminar viscosity at cell center.
  double             *LambdaOverCp;                ///< Thermal diffusion coefficient at cell center.
  double             *chi;                         ///< Scalar dissipation rate (without ZVar).
  double             *CmeanSource;                 ///< Source term for progress variable C.
  double             Schmidt_turb_Z;               ///< Turbulent Schmidt number for mixture fraction Z.
  double             Schmidt_turb_Z2;              ///< Turbulent Schmidt number for variance of mixture fraction.
  double             Schmidt_turb_C;               ///< Turbulent Schmidt number for progress variable C.
  double             Cchi;                         ///< Constant parameter for scalar dissipation rate.
  double             Cmu;                          ///< Constant parameter for scalar dissipation rate with k-omega model.
  double             Preference;                   ///< Reference pressure at which chemistry table was computed.
  int                ZMean_Index;                  ///< Index in scalar vector of ZMean.
  int                ZVar_Index;                   ///< Index in scalar vector of ZVar.
  int                CMean_Index;                  ///< Index in scalar vector of CMean.
  vector<double*>    OutputSpecies;                ///< Vector of pointers to scalars/species to be saved (e.g., for output).
  vector<string>     OutputSpeciesName;            ///< Vector with name of scalars/species to be saved (e.g., for output).  
  string             Comb_Reg;                     ///< Combustion regime (cold/hot/ignition) to set the progress variable.
  double             *dCmeanSource_dZM;            ///< Derivative of progress variable source term with respect to mean mixture fraction.
  double             *dCmeanSource_dZV;            ///< Derivative of progress variable source term with respect to mixture fraction variance.
  double             *dCmeanSource_dCM;            ///< Derivative of progress variable source term with respect to mean progress variable.
  
	// UQ
	double             *Source_Pert;						
  double             *Temp_Pert;									
	double             Source_PertFactor;      
  double             Temp_PertFactor;							
  string             UQTable;                     

private:
  double T0;
  double E0;
  double GAMMA0;
  double AGAMMA;
  double MU0;
  double AMU;
  double LOC0;
  double ALOC;

public:
  // constructors

  RansCombFPVA_Coeff()
  {    
    if (mpi_rank == 0)
      cout << "RansCombFPVA_Coeff()" << endl;

    Comb_Reg = getStringParam("COMBUSTION_REGIME", "HOT");
		
		// UQ Parameters
    UQTable = getStringParam("UQ_TABLE", "NO");
		Source_PertFactor = getDoubleParam("SOURCE_PERT_FACTOR", "0.0");
    Temp_PertFactor = getDoubleParam("TEMP_PERT_FACTOR", "0.0");

		
    // Overall constants
    Cmu = 0.09;
    Cchi            = getDoubleParam("Cchi", "2.0");
    Schmidt_turb_Z  = getDoubleParam("SCTURB_ZMean", "0.5");
    Schmidt_turb_Z2 = getDoubleParam("SCTURB_ZVar", "0.5");
    Schmidt_turb_C  = getDoubleParam("SCTURB_CMean", "0.5");    
    if (mpi_rank == 0) 
    {
      cout << "Turbulent Schmidt number for ZMean Sc_t = " << Schmidt_turb_Z  << endl;
      cout << "Turbulent Schmidt number for ZVar Sc_t2 = " << Schmidt_turb_Z2 << endl;
      cout << "Turbulent Schmidt number for CMean Sc_t = " << Schmidt_turb_C  << endl;
    }
    
    // Read name of additional scalars/species to save (e.g., for output) and register corresponding scalar
    Param * p;
    if (getParam(p, "OUTPUT_SPECIES"))
    {
      OutputSpeciesName.resize(p->getSize()-1);
      OutputSpecies.resize(p->getSize()-1);
      for (int i=0; i<p->getSize()-1; i++)
      {
        OutputSpecies[i] = NULL;
        OutputSpeciesName[i] = p->getString(i + 1);
        registerScalar(OutputSpecies[i], OutputSpeciesName[i], CV_DATA);
      }
    }
    
    // Scalar transport equations for ZMean, ZVar and CMean
    ScalarTranspEq *eq;
    eq = registerScalarTransport("ZMean", CV_DATA);
    eq->relax = getDoubleParam("RELAX_ZMean", "1.0");
    eq->diffTerm = getIntParam("SCAL_DIFF", "1");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 100;
    eq->lowerBound = -10.0;
    eq->upperBound = 10.0;
    eq->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");
    eq->coupling = "COUPLED";

    eq = registerScalarTransport("ZVar", CV_DATA);
    eq->relax = getDoubleParam("RELAX_ZVar", "1.0");
    eq->diffTerm = getIntParam("SCAL_DIFF", "1");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 100;
    eq->lowerBound = -10.0;
    eq->upperBound = 10.0;
    eq->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");
    eq->coupling = "COUPLED";

    eq = registerScalarTransport("CMean", CV_DATA);
    eq->relax = getDoubleParam("RELAX_CMean", "1.0");
    eq->diffTerm = getIntParam("SCAL_DIFF", "1");
    eq->phiZero = 1.0e-8;
    eq->phiZeroRel = 1.0e-2;
    eq->phiMaxiter = 100;
    eq->lowerBound = -10.0;
    eq->upperBound = 10.0;
    eq->reconstruction = getStringParam("SCALAR_RECONSTRUCTION", "CONSERVATIVE");
    eq->coupling = "COUPLED";
    
    // Other cell-centered quantities
//    rhoDkZ       = NULL;        registerScalar(rhoDkZ, "rhoDkZ", CV_DATA);
    muLam        = NULL;        registerScalar(muLam, "muLam" , CV_DATA);
    LambdaOverCp = NULL;        registerScalar(LambdaOverCp, "LambdaOverCp", CV_DATA);
    chi          = NULL;        registerScalar(chi, "chi", CV_DATA);
    CmeanSource  = NULL;        registerScalar(CmeanSource, "CmeanSource", CV_DATA);
		
		if (UQTable == "YES") {
			Source_Pert  = NULL;        registerScalar(Source_Pert, "Source_Pert", CV_DATA);
			Temp_Pert  = NULL;							 registerScalar(Temp_Pert, "Temp_Pert", CV_DATA);
		}

  }

  /*! \brief Destructor: unload chemistry table and clean memory. */
  virtual ~RansCombFPVA_Coeff()
  {
    myChemTable.Unload();
    
    if (scalarTranspEqVector[ZMean_Index].dpress_dphi != NULL)
    {
      delete [] scalarTranspEqVector[ZMean_Index].dpress_dphi; 
      scalarTranspEqVector[ZMean_Index].dpress_dphi = NULL;
    }
    if (scalarTranspEqVector[ZVar_Index].dpress_dphi != NULL)
    {
      delete [] scalarTranspEqVector[ZVar_Index].dpress_dphi; 
      scalarTranspEqVector[ZVar_Index].dpress_dphi = NULL;
    }
    if (scalarTranspEqVector[CMean_Index].dpress_dphi != NULL)
    {
      delete [] scalarTranspEqVector[CMean_Index].dpress_dphi; 
      scalarTranspEqVector[CMean_Index].dpress_dphi = NULL;
    }
    
    if (dCmeanSource_dZM != NULL)
    {
      delete [] dCmeanSource_dZM;
      dCmeanSource_dZM = NULL;
    }
    if (dCmeanSource_dZV != NULL)
    {
      delete [] dCmeanSource_dZV;
      dCmeanSource_dZV = NULL;
    }
    if (dCmeanSource_dCM != NULL)
    {
      delete [] dCmeanSource_dCM;
      dCmeanSource_dCM = NULL;
    }
    
    if (mpi_rank == 0)
      cout << endl << "***** FPVA_Coeff Combustion finalized *****" << endl << endl;
  }

public:
  // member functions 

  /*! \brief Read thermo input file and load chemistry table. */
  virtual void initialHookScalarRansCombModel()
  {
    myChemTable.Load(getStringParam("CHEMTABLE_FILE"));

    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");          ZMean = eq->phi;           ZMean_Index = getScalarTransportIndex("ZMean");
    eq = getScalarTransportData("ZVar");           ZVar  = eq->phi;           ZVar_Index  = getScalarTransportIndex("ZVar");
    eq = getScalarTransportData("CMean");          CMean = eq->phi;           CMean_Index = getScalarTransportIndex("CMean");
    
    Preference = getDoubleParam("REFERENCE_PRESSURE", "0.0");
    if (Preference == 0.0)
      Preference = myChemTable.GetReferencePressure();
    
    // For coupled solution
    if ((getStringParam("TIME_INTEGRATION") == "BACKWARD_EULER_COUPLED") || (getStringParam("TIME_INTEGRATION") == "BACKWARD_EULER_SEMICOUPLED")
				 || (getStringParam("TIME_INTEGRATION") == "BDF2_SEMICOUPLED"))
    {
      scalarTranspEqVector[ZMean_Index].dpress_dphi = new double[ncv_g];
      scalarTranspEqVector[ZVar_Index ].dpress_dphi = new double[ncv_g];
      scalarTranspEqVector[CMean_Index].dpress_dphi = new double[ncv_g];
      dCmeanSource_dZM = new double[ncv_g];
      dCmeanSource_dZV = new double[ncv_g];
      dCmeanSource_dCM = new double[ncv_g];
    }

    if (mpi_rank == 0)
      cout << endl << "***** FPVA_Coeff Combustion initialized *****" << endl << endl;
  }

  /*! \brief Compute for each cell the mixture properties as function of \p Zmean, \p Zvar, \p chi.
   *
   *  Compute composition Species::Yk, temperature, viscosity Mixture::mul, heat coefficient Mixture::lambda,
   *  enthalpy Mixture::H, heat capacity ratio Mixture::gama, diffusion coefficient, pressure and gas constant for the mixture.
   */
  virtual void calcStateVariables()
  {
    double E, cp;
    double kinecv = 0.0;
		double CmeanSource_n, CmeanSource_np1;
	  
		// Relaxation of the source term
		double relaxation = 1.0 ;
		if (Comb_Reg == "HOT_FILTER")
			relaxation = getDoubleParam("SOURCE_RELAXATION", "1.0");
		
    ComputeScalarDissipation(chi);
	  
    for (int icv = 0; icv < ncv; icv++)
    {
		
			// If source term filter, then save the old solution
			if (Comb_Reg == "HOT_FILTER") {
				if ((relaxation > 0.0) && (relaxation < 1.0)) 
					CmeanSource_n = CmeanSource[icv];
			}
			
      // Read and interpolate coefficients from chemistry table and save source term for progress variable
      myChemTable.LookupCoeff(RoM[icv], T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, CmeanSource[icv], ZMean[icv], ZVar[icv], CMean[icv]);

			if (UQTable == "YES") {
				Temp_Pert[icv] = Temp_PertFactor * sqrt(max( myChemTable.Lookup(ZMean[icv], ZVar[icv], CMean[icv], "TEMPCFD_VAR") ,1e-10)) ;
				Source_Pert[icv] = Source_PertFactor * sqrt(max( myChemTable.Lookup(ZMean[icv], ZVar[icv], CMean[icv], "SRC_PROG_VAR") ,1e-10)) ;
				
				// Temperature and source term clipping
				T0 = max( T0 + Temp_Pert[icv] , 0.0 ) ;
				CmeanSource[icv] = max( CmeanSource[icv] + Source_Pert[icv] , 0.0 ) ;
			}
      
			if ((Comb_Reg == "COLD") || (Comb_Reg == "COLD_SOURCE"))
        CmeanSource[icv] = 0.0;
		
			// Relaxation of the source term
			if (Comb_Reg == "HOT_FILTER") {
				if ((relaxation > 0.0) && (relaxation < 1.0)) {
					CmeanSource_np1 = CmeanSource[icv];
					CmeanSource[icv] = CmeanSource_n + (CmeanSource_np1 - CmeanSource_n)*relaxation;
				}
			}
      
      for (int i=0; i<3; i++)
        vel[icv][i] = rhou[icv][i]/rho[icv];

      if (kine != NULL)
        kinecv = kine[icv];
      
      E = rhoE[icv] / rho[icv] - 0.5 * vecDotVec3d(vel[icv], vel[icv]) - kinecv;

      // Compute mixture temperature from energy and species mass fractions
      if ( fabs(AGAMMA) < 1.0e-8 ) 
        temp[icv] = T0 + (GAMMA0 - 1.0) / RoM[icv] * (E - E0);
      else
        temp[icv] = T0 + (GAMMA0 - 1.0) / AGAMMA * (exp(AGAMMA * (E - E0) / RoM[icv]) - 1.0);
      
      // Temperature clipping
      temp[icv] = max(temp[icv], Tminimum);		

      // Compute mixture enthalpy and pressure
      enthalpy[icv] = E + RoM[icv] * temp[icv];
      press[icv] = rho[icv] * RoM[icv] * temp[icv];
          
      // Compute mixture viscosity, heat conductivity, gamma, speed of sound
      muLam[icv] = MU0 * pow(temp[icv] / T0, 0.7);
      LambdaOverCp[icv] = LOC0 * pow(temp[icv] / T0, 0.62);
      gamma[icv] = GAMMA0 + AGAMMA * (temp[icv] - T0);
      sos[icv] = sqrt(gamma[icv] * press[icv] / rho[icv]);

      if (isnan(press[icv]) || (press[icv]<= 0.0))
      {
        cout << "WARNING! :" << endl;
        cout << "x = " << x_cv[icv][0] << " / " << x_cv[icv][1]  << " / " << x_cv[icv][2] << endl;
        cout << "press = " << press[icv] << "  rho = " << rho[icv] << "  temp = " << temp[icv] << "  RoM = " << RoM[icv] << "  E = " << E <<
        "  rhoE = " << rhoE[icv] << "  1/2u^2 = " << 0.5 * vecDotVec3d(vel[icv], vel[icv]) << "  Z = " << ZMean[icv] << endl;
        throw(-1);
      }
      
      // Save additional species/scalars for output
      for (int i=0; i<OutputSpeciesName.size(); i++)
      {
        string SpeciesName = OutputSpeciesName[i];
        double *Species = getScalar(SpeciesName.c_str());  
        Species[icv] = myChemTable.Lookup(ZMean[icv], ZVar[icv], CMean[icv], SpeciesName);
      }
      
#ifdef PRESSURE_SCALING_FPVA_COEFF
      // rescale source term
      CmeanSource[icv] *= pow( press[icv]/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif      

    }
	  

	  // Implicit smoothing of the progresive variable
	  if (Comb_Reg == "HOT_FILTER") {
		  int nSmooth = getIntParam("NSMOOTH_MEANSOURCE", "1");
		  double smooth_coeff = 0.75;
		  double *CmeanSource_old = new double[ncv];
		  double *CmeanSource_sum = new double[ncv];
		  
		  for (int icv = 0; icv < ncv; icv++)
			  CmeanSource_old[icv] = CmeanSource[icv];
		  
		  for (int iSmooth = 0; iSmooth < nSmooth; iSmooth++) {
			  
			  for (int icv = 0; icv < ncv; icv++)
				  CmeanSource_sum[icv] = 0.0;
			  
			  for (int ifa = nfa_b; ifa < nfa; ifa++) {
				  int icv0 = cvofa[ifa][0]; int icv1 = cvofa[ifa][1];
				  if (icv1 < ncv) {
					  CmeanSource_sum[icv0] += CmeanSource[icv1]; 
					  CmeanSource_sum[icv1] += CmeanSource[icv0];
				  }
			  }
			  
			  for (int icv = 0; icv < ncv; icv++) {
				  int non_f = nbocv_i[icv];
				  int non_l = nbocv_i[icv+1]-1;
				  double nneigh = double(non_l-non_f);
				  CmeanSource[icv] = (CmeanSource_old[icv] + smooth_coeff*CmeanSource_sum[icv]) / (1.0 + smooth_coeff*nneigh);
			  }
			  
			  for (int ifa = 0; ifa < nfa_b; ifa++) {
				  int icv0 = cvofa[ifa][0];
				  CmeanSource[icv0] = CmeanSource_old[icv0];
			  }
				
		  }
		  delete [] CmeanSource_old;
		  delete [] CmeanSource_sum;
	  }
	  
	  
    updateCvData(vel, REPLACE_ROTATE_DATA);
    updateCvData(temp, REPLACE_DATA);
    updateCvData(RoM, REPLACE_DATA);
    updateCvData(enthalpy, REPLACE_DATA);
    updateCvData(press, REPLACE_DATA);
    updateCvData(gamma, REPLACE_DATA);
    updateCvData(sos, REPLACE_DATA);
    updateCvData(muLam, REPLACE_DATA);
    updateCvData(LambdaOverCp, REPLACE_DATA);
    updateCvData(CmeanSource, REPLACE_DATA);
//    updateCvData(rhoDkZ, REPLACE_DATA);
    
    for (int i=0; i<OutputSpeciesName.size(); i++)
    {
      string SpeciesName = OutputSpeciesName[i];
      double *Species = getScalar(SpeciesName.c_str());  
      updateCvData(Species, REPLACE_DATA);
    }    
	  
  }
  
#ifdef MIXING_VISCOSITY  
  /*! \brief Compute for each face the material properties.
   *
   *  Compute laminar viscosity mul_fa and heat conductivity lamOcp_fa at the faces.
   *  
   */
  virtual void calcMaterialProperties()
  {
    if (mu_ref > 0.0)
    {
      // internal faces
      for (int ifa = nfa_b; ifa < nfa; ifa++)
      {
        int icv0 = cvofa[ifa][0];
        int icv1 = cvofa[ifa][1];
  
        double dx0[3] = {0.0, 0.0, 0.0}, dx1[3] = {0.0, 0.0, 0.0};
        vecMinVec3d(dx0, x_fa[ifa], x_cv[icv0]);
        vecMinVec3d(dx1, x_fa[ifa], x_cv[icv1]);
        double w0 = sqrt(vecDotVec3d(dx0, dx0));
        double w1 = sqrt(vecDotVec3d(dx1, dx1));
 
        mul_fa[ifa] = (w1*muLam[icv0] + w0*muLam[icv1])/(w0+w1);
        lamOcp_fa[ifa] = (w1*LambdaOverCp[icv0] + w0*LambdaOverCp[icv1])/(w0+w1);
      }    
      
      // boundary faces computed in setBC
    }
  }
#endif
  
  // \brief Compute derivative of pressure with respect to scalar for coupled solution (Jacobi matrix)
  virtual void pressureDerivativeHookScalarRansComb() 
  {
    //  set to 0 to agree with old version where variable was passed by value, not reference...
    //  setting dpress_dphi to 0 also improves stability drastically -- should figure out why...
    for (int icv = 0; icv < ncv; icv++) {
      scalarTranspEqVector[ZMean_Index].dpress_dphi[icv] = 0.0 ;
      scalarTranspEqVector[ZVar_Index].dpress_dphi[icv] = 0.0 ;
      scalarTranspEqVector[CMean_Index].dpress_dphi[icv] = 0.0 ;
    }

    if ((Comb_Reg == "COLD") || (Comb_Reg == "COLD_SOURCE")) {
      for (int icv = 0; icv < ncv; icv++) {
        dCmeanSource_dZM[icv] = 0.0 ;
        dCmeanSource_dZV[icv] = 0.0 ;
        dCmeanSource_dCM[icv] = 0.0 ;
      }
    }
    else {   //  reacting flow
      double delta_ZM = 1.0e-5;
      double delta_ZV = 1.0e-8;
      double delta_CM = 1.0e-5;
      double pressp, pressm, CmeanSourcep, CmeanSourcem;
      for (int icv = 0; icv < ncv; icv++) {
        double E = enthalpy[icv] - RoM[icv] * temp[icv];
      
        // Derivative with respect to ZMean
        double ZMp = ZMean[icv] + delta_ZM;
        double ZMm = ZMean[icv] - delta_ZM;
        pressure_scalSource(pressp, CmeanSourcep, rho[icv], E, ZMp, ZVar[icv], CMean[icv]);
        pressure_scalSource(pressm, CmeanSourcem, rho[icv], E, ZMm, ZVar[icv], CMean[icv]);
#ifdef PRESSURE_SCALING_FPVA_COEFF
        CmeanSourcep *= pow( pressp/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
        CmeanSourcem *= pow( pressm/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif      
        dCmeanSource_dZM[icv] = (CmeanSourcep - CmeanSourcem) / (2.0 * delta_ZM);
        //scalarTranspEqVector[ZMean_Index].dpress_dphi[icv] = (pressp - pressm) / (2.0 * delta_ZM);
          
        // Derivative with respect to ZVar
        double ZVp = ZVar[icv] + delta_ZV;
        double ZVm = ZVar[icv] - delta_ZV;
        pressure_scalSource(pressp, CmeanSourcep, rho[icv], E, ZMean[icv], ZVp, CMean[icv]);
        pressure_scalSource(pressm, CmeanSourcem, rho[icv], E, ZMean[icv], ZVm, CMean[icv]);
#ifdef PRESSURE_SCALING_FPVA_COEFF
        CmeanSourcep *= pow( pressp/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
        CmeanSourcem *= pow( pressm/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif      
        dCmeanSource_dZV[icv] = (CmeanSourcep - CmeanSourcem) / (2.0 * delta_ZV);
        //scalarTranspEqVector[ZVar_Index].dpress_dphi[icv] = (pressp - pressm) / (2.0 * delta_ZV);
      
        // Derivative with respect to CMean
        double CMp = CMean[icv] + delta_CM;
        double CMm = CMean[icv] - delta_CM;
        pressure_scalSource(pressp, CmeanSourcep, rho[icv], E, ZMean[icv], ZVar[icv], CMp);
        pressure_scalSource(pressm, CmeanSourcem, rho[icv], E, ZMean[icv], ZVar[icv], CMm);
#ifdef PRESSURE_SCALING_FPVA_COEFF
        CmeanSourcep *= pow( pressp/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
        CmeanSourcem *= pow( pressm/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif      
        dCmeanSource_dCM[icv] = (CmeanSourcep - CmeanSourcem) / (2.0 * delta_CM);
        //scalarTranspEqVector[CMean_Index].dpress_dphi[icv] = (pressp - pressm) / (2.0 * delta_CM);
      }
    }

    updateCvData(scalarTranspEqVector[ZMean_Index].dpress_dphi, REPLACE_DATA);
    updateCvData(scalarTranspEqVector[ZVar_Index].dpress_dphi, REPLACE_DATA);
    updateCvData(scalarTranspEqVector[CMean_Index].dpress_dphi, REPLACE_DATA);
    updateCvData(dCmeanSource_dZM, REPLACE_DATA);
    updateCvData(dCmeanSource_dZV, REPLACE_DATA);
    updateCvData(dCmeanSource_dCM, REPLACE_DATA);
  } 
  
  void pressure_scalSource(double &pp, double &CMSource, double &rrho, double &E, double &ZM, double &ZV, double &CM)
  {
    double R0, T;
    myChemTable.LookupSelectedCoeff(R0, T0, E0, GAMMA0, AGAMMA, ZM, ZV, CM);
		
		if (UQTable == "YES") {
			double Temp_Pert = Temp_PertFactor * sqrt(max( myChemTable.Lookup(ZM, ZV, CM, "TEMPCFD_VAR") ,1e-10)) ;
			
			// Temperature clipping
			T0 = max( T0 + Temp_Pert, 0.0 ) ;
		}
		
    CMSource = myChemTable.Lookup(ZM, ZV, CM, "SRC_PROG");
		
		if (UQTable == "YES") {
			double Source_Pert = Source_PertFactor * sqrt(max( myChemTable.Lookup(ZM, ZV, CM, "SRC_PROG_VAR") ,1e-10)) ;
			
			// Source Term clipping
			CMSource = max( CMSource + Source_Pert, 0.0 ) ;
		}
		
    if ( fabs(AGAMMA) < 1.0e-8 ) 
      T = T0 + (GAMMA0 - 1.0) / R0 * (E - E0);
    else
      T = T0 + (GAMMA0 - 1.0) / AGAMMA * (exp(AGAMMA * (E - E0) / R0) - 1.0);
    
    // Temperature clipping
    T = max(T, Tminimum);
    pp = rrho * R0 * T;
  }
  
  // \brief Compute for a given temperature, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_T(double &p, double &h, double &R, double &gam, double &rho, double &T, double *Scal, int nScal)
  {
    double zm = Scal[ZMean_Index];
    double zv = Scal[ZVar_Index];
    double cm = Scal[CMean_Index];
    
    myChemTable.LookupSelectedCoeff(R, T0, E0, GAMMA0, AGAMMA, zm, zv, cm);
		
		if (UQTable == "YES") {
			double Temp_Pert = Temp_PertFactor * sqrt(max( myChemTable.Lookup(zm, zv, cm, "TEMPCFD_VAR") ,1e-10)) ;
			
			// Temperature clipping
			T0 = max( T0 + Temp_Pert, 0.0 ) ;
		}
		
    if ( fabs(AGAMMA) < 1.0e-8 ) 
      h = E0 + R / (GAMMA0 - 1.0) * (T - T0) + R * T;
    else
      h = E0 + R / AGAMMA * log(1.0 + AGAMMA * (T - T0) / (GAMMA0 - 1.0)) + R * T;
    gam = GAMMA0 + AGAMMA * (T - T0);
    p = rho * R * T;
  }
 
  // \brief Compute for a given pressure, density and scalars: pressure, enthalpy, gas constant and ratio of specific heat
  virtual void calcThermoProp_p(double &T, double &h, double &R, double &gam, double &rho, double &p, double *Scal, int nScal)
  {
    double zm = Scal[ZMean_Index];
    double zv = Scal[ZVar_Index];
    double cm = Scal[CMean_Index];
    
    myChemTable.LookupSelectedCoeff(R, T0, E0, GAMMA0, AGAMMA, zm, zv, cm);
		
		if (UQTable == "YES") {
			double Temp_Pert = Temp_PertFactor * sqrt(max( myChemTable.Lookup(zm, zv, cm, "TEMPCFD_VAR") ,1e-10)) ;
			
			// Temperature clipping
			T0 = max( T0 + Temp_Pert, 0.0 ) ;
		}
		
    T = p / (rho * R);
    if ( fabs(AGAMMA) < 1.0e-8 ) 
      h = E0 + R / (GAMMA0 - 1.0) * (T - T0) + R * T;
    else    
      h = E0 + R / AGAMMA * log(1.0 + AGAMMA * (T - T0) / (GAMMA0 - 1.0)) + R * T;
    gam = GAMMA0 + AGAMMA * (T - T0);
  }
  
  virtual double calcMuLam(int icv)
  {
    return muLam[icv];
  }
  
  virtual double calcMuLam(double temp)
  {
    cout << "### RansCombFPVA::calcMuLam(double temp) does not work! ####" << endl;
    throw(-1);
  }
  
  // \brief Compute for a given temperature the properties of the mixture at a face.
  virtual void ComputeBCProperties_T(FaZone *zone)
  {
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");       double *zM_fa = eq->phi_bfa;
    eq = getScalarTransportData("ZVar");        double *zV_fa = eq->phi_bfa;
    eq = getScalarTransportData("CMean");       double *cM_fa = eq->phi_bfa;

    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      ComputeProperties_T(h_bfa[ifa],RoM_bfa[ifa],gam_bfa[ifa],mul_fa[ifa],lamOcp_fa[ifa],T_bfa[ifa],zM_fa[ifa],zV_fa[ifa],cM_fa[ifa]);
  }
  
  // \brief Compute for a given enthalpy the properties of the mixture at a face.
  virtual void ComputeBCProperties_H(FaZone *zone)
  {    
    ScalarTranspEq *eq;
    eq = getScalarTransportData("ZMean");       double *zM_fa = eq->phi_bfa;
    eq = getScalarTransportData("ZVar");        double *zV_fa = eq->phi_bfa;
    eq = getScalarTransportData("CMean");       double *cM_fa = eq->phi_bfa;

    for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
      ComputeProperties_H(T_bfa[ifa],RoM_bfa[ifa],gam_bfa[ifa],mul_fa[ifa],lamOcp_fa[ifa],h_bfa[ifa],zM_fa[ifa],zV_fa[ifa],cM_fa[ifa]);
  }
  
  /*! \brief Compute for a given temperature the properties of the mixture.
   *
   *  The enthalpy (chemical + sensible) is a function of temperature and mixture composition. 
   *  \param[out] H Enthalpy (chemical + sensible) in [J/kg].
   *  \param[out] RoverM Gas constant R/M of the mixture in [J/(kg K)].
   *  \param[out] gam Heat capacity ratio in [-].
   *  \param[out] mul Laminar viscosity in [kg/(m s)].
   *  \param[out] lam_o_cp Heat conductivity divided by heat capacity in [kg/(m s)].
   *  \param[in]  T Temperature in [K].
   *  \param[in]  zm Mean mixture fraction.
   *  \param[in]  zv Variance of mixture fraction.
   *  \param[in]  chi Scalar dissipation rate.
   */
  void ComputeProperties_T(double &h, double &RoM, double &gam, double &mu, double &lamOcp, double T, double zm, double zv, double cm)
  {
    double dummy;

    myChemTable.LookupCoeff(RoM, T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, dummy, zm, zv, cm);

		if (UQTable == "YES") {
			double Temp_Pert = Temp_PertFactor * sqrt(max( myChemTable.Lookup(zm, zv, cm, "TEMPCFD_VAR") ,1e-10)) ;
			
			// Temperature clipping
			T0 = max( T0 + Temp_Pert, 0.0 ) ;
		}
		
    if ( fabs(AGAMMA) < 1.0e-8 ) 
      h = E0 + RoM / (GAMMA0 - 1.0) * (T - T0) + RoM * T;
    else    
      h = E0 + RoM / AGAMMA * log(1.0 + AGAMMA * (T - T0) / (GAMMA0 - 1.0)) + RoM * T;
    gam = GAMMA0 + AGAMMA * (T - T0);
    mu = MU0 * pow(T / T0, 0.7);
    lamOcp = LOC0 * pow(T / T0, 0.62);
  }

  /*! \brief Compute for a given enthalpy the properties of the mixture.
   *
   *  The enthalpy (chemical + sensible) is a function of temperature and mixture composition. 
   *  \param[out] T Temperature in [K].
   *  \param[out] RoverM Gas constant R/M of the mixture in [J/(kg K)].
   *  \param[out] gam Heat capacity ratio in [-].
   *  \param[out] mul Laminar viscosity in [kg/(m s)].
   *  \param[out] lam_o_cp Heat conductivity divided by heat capacity in [kg/(m s)].
   *  \param[in]  H Enthalpy (chemical + sensible) in [J/kg].
   *  \param[in]  zm Mean mixture fraction.
   *  \param[in]  zv Variance of mixture fraction.
   *  \param[in]  chi Scalar dissipation rate.
   */
  void ComputeProperties_H(double &T, double &RoM, double &gam, double &mu, double &lamOcp, double h, double zm, double zv, double cm)
  {
    double dummy;

    myChemTable.LookupCoeff(RoM, T0, E0, GAMMA0, AGAMMA, MU0, AMU, LOC0, ALOC, dummy, zm, zv, cm);

		if (UQTable == "YES") {
			double Temp_Pert = Temp_PertFactor * sqrt(max( myChemTable.Lookup(zm, zv, cm, "TEMPCFD_VAR") ,1e-10)) ;
			
			// Temperature clipping
			T0 = max( T0 + Temp_Pert, 0.0 ) ;
		}		
		
    T = (h - E0) / RoM;
    T = SolveNewtonTemperature_H(RoM, T0, E0, GAMMA0, AGAMMA, h, T);
    T = max(T, Tminimum);
    gam = GAMMA0 + AGAMMA * (T - T0);
    mu = MU0 * pow(T / T0, 0.7);
    lamOcp = LOC0 * pow(T / T0, 0.62);
  }
  
  double SolveNewtonTemperature_H(double R, double t0, double e0, double gamma0, double agamma, double h, double Tguess)
  {
    int    itermax = 50, iter;
    double errormax = 1.0e-7, error;
    double Tresidual, hh, cpp;

    iter = 0;
    do
    {
      ++iter;
      if ( fabs(agamma) < 1.0e-8 ) 
        hh = e0 + R / (gamma0 - 1.0) * (Tguess - t0) + R * Tguess;
      else    
        hh = e0 + R / agamma * log(1.0 + agamma * (Tguess - t0) / (gamma0 - 1.0)) + R * Tguess;
      cpp = R * (1.0 + 1.0 / (gamma0 + agamma * (Tguess - t0) - 1.0));
    
      Tresidual = Tguess - (hh - h) / cpp;

      error = Tguess - Tresidual;
      if (fabs(error) <= errormax)
        return Tresidual;

      Tguess = Tresidual;

    }
    while (iter < itermax);
    
    /* if the iteration has not converged, exit with error */
    cerr << "### Computation of temperature in Newton iteration has not converged: Tguess=" << Tguess + error 
         << ", Tresidual=" << Tresidual << ", enthalpy=" << h << ", agamma=" << agamma << ", #iterations=" << iter << " ###" << endl;
    throw(-1);
  }
  
  void ComputeScalarDissipation(double *chi)
  {
    if ((turbModel == KOM) || (turbModel == SST))
    {
      double *om = getR1("omega");
      for (int icv = 0; icv < ncv; icv++)
        chi[icv] = Cchi * Cmu * om[icv];
    }
    
    else if (turbModel == KEPS)
    {
      double *epsi = getR1("eps");
      for (int icv = 0; icv < ncv; icv++)
        chi[icv] = Cchi * epsi[icv] / kine[icv];
    }
    
    else if (turbModel == SA)
    {
      calcVorticity();
      for (int icv = 0; icv < ncv; icv++)
        chi[icv] = Cchi * sqrt(Cmu) * vortMag[icv];
    }
    
    else
    {
      for (int icv = 0; icv < ncv; icv++)
        chi[icv] = 0.0;
    }
    
    updateCvData(chi, REPLACE_DATA);
  }

  virtual void diffusivityHookScalarRansComb(const string &name)
  {
    ScalarTranspEq *eq;

    if ((name == "ZMean") || (name == "ZVar") || (name == "CMean"))
    {
      eq = getScalarTransportData(name);
 
      // internal faces
      for (int ifa = nfa_b; ifa < nfa; ifa++)
        eq->diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;    

      // boundary faces
      for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        if (zoneIsWall(zone->getName()))
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            eq->diff[ifa] = lamOcp_fa[ifa];
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
            eq->diff[ifa] = lamOcp_fa[ifa] + mut_fa[ifa]/Schmidt_turb_Z;
      }
    }
  }
  
	virtual void sourceHookScalarRansComb_new(double *rhs, double *A, const string &name, int flagImplicit)
	{
		if (name == "ZVar")
		{
			ScalarTranspEq *eq = getScalarTransportData("ZMean");       
			double *ZMean_bfa = eq->phi_bfa;  
			double (*ZMean_grad)[3] = eq->grad_phi;
			
			for (int icv = 0; icv < ncv; icv++)
			{
				double src = 2.0*InterpolateAtCellCenterFromFaceValues(mut_fa, icv)/Schmidt_turb_Z2 * vecDotVec3d(ZMean_grad[icv], ZMean_grad[icv]) - chi[icv] * rho[icv] * ZVar[icv];
				rhs[icv] += src * cv_volume[icv];
			}
			
			if (flagImplicit)
			{
				for (int icv = 0; icv < ncv; icv++)
				{
					int noc00 = nbocv_i[icv];
					double dsrcdphi = - chi[icv];
					A[noc00] -= dsrcdphi*cv_volume[icv];
				}
			}
		}
		
		if (name == "CMean")
		{
			for (int icv = 0; icv < ncv; icv++)
			{
				rhs[icv] += rho[icv] * CmeanSource[icv] * cv_volume[icv];
			}
			
			if (flagImplicit) {
				// Derivative with respect to CMean
				for (int icv = 0; icv < ncv; icv++) {
					double E = enthalpy[icv] - RoM[icv] * temp[icv];
					double delta_CM = 1.0e-5;
					double CMp = CMean[icv] + delta_CM;
					double CMm = CMean[icv] - delta_CM;
					double pressp, pressm, CmeanSourcep, CmeanSourcem;
					pressure_scalSource(pressp, CmeanSourcep, rho[icv], E, ZMean[icv], ZVar[icv], CMp);
					pressure_scalSource(pressm, CmeanSourcem, rho[icv], E, ZMean[icv], ZVar[icv], CMm);
#ifdef PRESSURE_SCALING_FPVA_COEFF
					CmeanSourcep *= pow( pressp/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
					CmeanSourcem *= pow( pressm/Preference , PRESSURE_SCALING_FPVA_COEFF ) ;
#endif      
					int noc00 = nbocv_i[icv];
					A[noc00] -= (CmeanSourcep - CmeanSourcem) / (2.0 * delta_CM) * cv_volume[icv] ;
				}
			}
		}
	}
	
	
	
  virtual void sourceHookRansCombCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
    ScalarTranspEq *eq = getScalarTransportData("ZMean");       
    double (*ZMean_grad)[3] = eq->grad_phi;
    
    
    int ZMean_Coupled_Index = getScalarTransportCoupledIndex("ZMean");
    if (ZMean_Coupled_Index == -1)
    {
      cerr << "### Error ZMean_Coupled_Index is " << ZMean_Coupled_Index << " !###" << endl;
      throw(-1);
    }
    int ZVar_Coupled_Index = getScalarTransportCoupledIndex("ZVar");
    if (ZVar_Coupled_Index == -1)
    {
      cerr << "### Error ZVar_Coupled_Index is " << ZVar_Coupled_Index << " !###" << endl;
      throw(-1);
    }
    int CMean_Coupled_Index = getScalarTransportCoupledIndex("CMean");
    if (CMean_Coupled_Index == -1)
    {
      cerr << "### Error CMean_Coupled_Index is " << CMean_Coupled_Index << " !###" << endl;
      throw(-1);
    }
    
    for (int icv = 0; icv < ncv; icv++)
    {
      rhs[icv][5+ZVar_Coupled_Index ] += (2.0 * InterpolateAtCellCenterFromFaceValues(mut_fa, icv) / Schmidt_turb_Z2 * vecDotVec3d(ZMean_grad[icv], ZMean_grad[icv]) - chi[icv] * rho[icv] * ZVar[icv]) * cv_volume[icv];
      rhs[icv][5+CMean_Coupled_Index] += rho[icv] * CmeanSource[icv] * cv_volume[icv];

      if (flagImplicit)
      {
        int noc00 = nbocv_i[icv];
        A[noc00][5+ZVar_Coupled_Index ][5+ZVar_Coupled_Index ] -= - chi[icv] * cv_volume[icv];
        A[noc00][5+CMean_Coupled_Index][5+ZMean_Coupled_Index] -= dCmeanSource_dZM[icv] * cv_volume[icv];
        A[noc00][5+CMean_Coupled_Index][5+ZVar_Coupled_Index ] -= dCmeanSource_dZV[icv] * cv_volume[icv];
        A[noc00][5+CMean_Coupled_Index][5+CMean_Coupled_Index] -= dCmeanSource_dCM[icv] * cv_volume[icv];
        A[noc00][5+CMean_Coupled_Index][0] -= (CmeanSource[icv] - ZMean[icv] * dCmeanSource_dZM[icv] - ZVar[icv] * dCmeanSource_dZV[icv] - CMean[icv] * dCmeanSource_dCM[icv]) * cv_volume[icv];
     }
    }
  }

};


#endif /* COMBMODEL_FPVA_COEFF_H_ */

