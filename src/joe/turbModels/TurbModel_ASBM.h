#ifndef RANSTURBMODEL_ASBM_H
#define RANSTURBMODEL_ASBM_H

#include "UgpWithCvCompFlow.h"
#include "myMem.h"

extern "C"{
  void asbm_(
      double*, double*, double*, double*, double*, double*,
      double*, double*, double*, int*, int*, int*);
}

class RansTurbASBM : virtual public UgpWithCvCompFlow
{
public: // constructor, destructor

  RansTurbASBM()
  {
    // General variables
    st_diag       = NULL;  registerVector(st_diag,    "st_diag",       CV_DATA);
    st_offdiag    = NULL;  registerVector(st_offdiag, "st_offdiag",    CV_DATA);
    wt_offdiag    = NULL;  registerVector(wt_offdiag, "wt_offdiag",    CV_DATA);

    as_diag       = NULL;  registerVector(as_diag,    "as_diag",       CV_DATA);
    as_offdiag    = NULL;  registerVector(as_offdiag, "as_offdiag",    CV_DATA);
    ar_diag       = NULL;  registerVector(ar_diag,    "ar_diag",       CV_DATA);
    ar_offdiag    = NULL;  registerVector(ar_offdiag, "ar_offdiag",    CV_DATA);

    // Blocking variables
    block_diag    = NULL;  registerVector(block_diag, "block_diag",    CV_DATA);
    block_offdiag = NULL;  registerVector(block_offdiag, "block_offdiag", CV_DATA);
    bphi          = NULL;  registerScalar(bphi,       "bphi",          CV_DATA);
    bphi_bfa      = NULL;  // this is a face array
    grad_bphi     = NULL;  // this array includes ghost cells

    // -------------------------------------------------------------------------------
    // Debugging variables
    marker = NULL;     // this is an integer array

    temp_muT      = NULL;     registerScalar(temp_muT,       "temp_muT",       CV_DATA);
    debug1        = NULL;     registerVector(debug1,         "debug1",         CV_DATA);
    debug2        = NULL;     registerVector(debug2,         "debug2",         CV_DATA);
    debug3        = NULL;     registerVector(debug3,         "debug3",         CV_DATA);
    etar          = NULL;     registerScalar(etar,           "etar",           CV_DATA);
    trace_aa      = NULL;     registerScalar(trace_aa,       "trace_aa",       CV_DATA);
    r_ratio       = NULL;     registerScalar(r_ratio,        "r_ratio",        CV_DATA);
    rij_diag_nd   = NULL;     registerVector(rij_diag_nd,    "rij_diag_nd",    CV_DATA);
    rij_offdiag_nd = NULL;    registerVector(rij_offdiag_nd, "rij_offdiag_nd", CV_DATA);
    // -------------------------------------------------------------------------------
  }

  virtual ~RansTurbASBM() {}

protected: // member variables

  // General variables
  int n_smooth;                   // # of times to smooth data

  double   (*st_diag)[3];         // diagonal rate of strain tensor times tau
  double   (*st_offdiag)[3];      // off-diagonal rate of strain tensor times tau
  double   (*wt_offdiag)[3];      // off-diagonal rate of rotation tensor times tau

  double   (*as_diag)[3];         // diagonal eddy axis tensor from st
  double   (*as_offdiag)[3];      // off-diagonal eddy axis tensor from st
  double   (*ar_diag)[3];         // diagonal eddy axis tensor from st and wt
  double   (*ar_offdiag)[3];      // off-diagonal eddy axis tensor from st and wt

  // Blocking variables
  double   (*block_diag)[3];      // diagonal blockage tensor
  double   (*block_offdiag)[3];   // off-diagonal blockage tensor
  double   *bphi;                 // phi for wall blocking
  double   *bphi_bfa;             // phi for wall blocking at boundaries
  double   (*grad_bphi)[3];       // gradients for wall blocking phi

  int      block_rij;             // block r_ij instead of a_ij
  int      block_frq;             // frequency of blocking computation
  int      block_maxIter;         // linear solver parameters
  double   block_zeroAbs;
  double   block_zeroRel;

  // I/O variables
  FILE *finfo;              // output file for ASBM quantities

  // -------------------------------------------------------------------------------
  // Debugging variables
  int *marker;

  double *temp_muT;
  double (*debug1)[3];
  double (*debug2)[3];
  double (*debug3)[3];
  double *etar;
  double *trace_aa;
  double *r_ratio;
  double (*rij_diag_nd)[3];
  double (*rij_offdiag_nd)[3];
  // -------------------------------------------------------------------------------

public:   // member functions

  virtual void initialHookScalarRansTurbModel()
  {
    n_smooth   = getIntParam("N_SMOOTH", "0");
    block_frq  = getIntParam("BLOCK_FRQ","10");
    block_rij  = getIntParam("BLOCK_RIJ", "0");

    if (!checkParam("LINEAR_SOLVER_BLOCK_TRESHOLDS"))
    {
      // add default values
      ParamMap::add("LINEAR_SOLVER_BLOCK_TRESHOLDS  MAX_ITER=500  ABS_RESID=1.0e-8  REL_RESID=1.0e-8");
      if (mpi_rank == 0)
        cout << "WARNING: added keyword"
             << "\"LINEAR_SOLVER_BLOCK_TRESHOLDS  MAX_ITER=500  ABS_RESID=1.0e-8  REL_RESID=1.0e-8\""
             << " to parameter map!" << endl;
    }
    block_maxIter = getParam("LINEAR_SOLVER_BLOCK_TRESHOLDS")->getInt("MAX_ITER");
    block_zeroAbs = getParam("LINEAR_SOLVER_BLOCK_TRESHOLDS")->getDouble("ABS_RESID");
    block_zeroRel = getParam("LINEAR_SOLVER_BLOCK_TRESHOLDS")->getDouble("REL_RESID");

    if (mpi_rank == 0)
    {
      if ( (finfo=fopen("asbmInfo.dat","wt")) == NULL )
      {
        cerr << "Error: cannot open file asbmInfo.dat" << endl;
        throw(-1);
      }
      else
        cout << "Opened file asbmInfo.dat" << endl;
    }

    // Initialize cell centered data
    marker = new int[ncv];
    for (int icv = 0; icv < ncv; ++icv)
    {
      marker[icv] = 0;

      as_diag[icv][0] = 1.0/3.0;
      as_diag[icv][1] = 1.0/3.0;
      as_diag[icv][2] = 1.0/3.0;

      ar_diag[icv][0] = 1.0/3.0;
      ar_diag[icv][1] = 1.0/3.0;
      ar_diag[icv][2] = 1.0/3.0;
    }

    // Initialize face centered data
    bphi_bfa  = new double[nfa_b];
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    {
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
        {
          if (param->getString() == "WALL")
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              bphi_bfa[ifa] = 1.0;
            }
          }
          else
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              int icv0 = cvofa[ifa][0];
              bphi_bfa[ifa] = bphi[icv0];
            }
          }
        }
      }
    }

    // Initialize blocking gradients
    grad_bphi = new double[ncv_g][3];
    calcCvScalarGrad(grad_bphi, bphi, bphi_bfa, gradreconstruction, limiterNavierS, bphi, epsilonSDWLS);
  }

  void calcRsCenterASBM(){
    // ====================================================================
    //        variable declaration and initialization
    // ====================================================================

    // Input variables
    double WFT[3][3] = {{0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
    double ST[3][3], WT[3][3], BL[3][3];
    double EMPTY[3][3] = {{0.0, 0.0, 0.0},
                          {0.0, 0.0, 0.0},
                          {0.0, 0.0, 0.0}};

    // In and out variables
    double AS[3][3], AR[3][3];

    // Output variables
    double REY[3][3] = {{0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0},
                        {0.0, 0.0, 0.0}};
    double DIM[3][3], CIR[3][3];

    // Error flags
    int ierr = 0, myNerr = 0, Nerr = 0, maxitrs = 99;

    // ====================================================================
    // compute inputs
    // ====================================================================
    for (int icv = 0; icv < ncv; icv++)
    {
      double div = diverg[icv];
      double tau = turbTS[icv];

      st_diag[icv][0] = (grad_u[icv][0][0] - 1.0/3.0*div)*tau;
      st_diag[icv][1] = (grad_u[icv][1][1] - 1.0/3.0*div)*tau;
      st_diag[icv][2] = (grad_u[icv][2][2] - 1.0/3.0*div)*tau;

      st_offdiag[icv][0] = 0.5*(grad_u[icv][0][1] + grad_u[icv][1][0])*tau;
      st_offdiag[icv][1] = 0.5*(grad_u[icv][0][2] + grad_u[icv][2][0])*tau;
      st_offdiag[icv][2] = 0.5*(grad_u[icv][1][2] + grad_u[icv][2][1])*tau;

      wt_offdiag[icv][0] = 0.5*(grad_u[icv][0][1] - grad_u[icv][1][0])*tau;
      wt_offdiag[icv][1] = 0.5*(grad_u[icv][0][2] - grad_u[icv][2][0])*tau;
      wt_offdiag[icv][2] = 0.5*(grad_u[icv][1][2] - grad_u[icv][2][1])*tau;
    }

    for (int i = 0; i < n_smooth; i++)
    {
      updateCvData(st_diag, REPLACE_ROTATE_DATA);
      updateCvData(st_offdiag, REPLACE_ROTATE_DATA);
      updateCvData(wt_offdiag, REPLACE_ROTATE_DATA);

      // Smooth inputs
      smoothingVec(st_diag);
      smoothingVec(st_offdiag);
      smoothingVec(wt_offdiag);
    }

    // ====================================================================
    // Compute Reynolds stresses
    // ====================================================================
    for (int icv = 0; icv < ncv; icv++)
    {
      // Anisotropic rate of strain tensor, sij = 0.5*(dui_dxj + duj_dxi) - 1/3*divg
      ST[0][0] = st_diag[icv][0];     ST[0][1] = st_offdiag[icv][0];     ST[0][2] = st_offdiag[icv][1];
      ST[1][0] = ST[0][1];            ST[1][1] = st_diag[icv][1];        ST[1][2] = st_offdiag[icv][2];
      ST[2][0] = ST[0][2];            ST[2][1] = ST[1][2];               ST[2][2] = st_diag[icv][2];

      // Rate of mean rotation tensor, C TO FORTRAN, NEED TRANSPOSE: wij = -0.5*(dui_dxj - duj_dxi)
      WT[0][0] = 0.0;                 WT[0][1] = -wt_offdiag[icv][0];    WT[0][2] = -wt_offdiag[icv][1];
      WT[1][0] = -WT[0][1];           WT[1][1] = 0.0;                    WT[1][2] = -wt_offdiag[icv][2];
      WT[2][0] = -WT[0][2];           WT[2][1] = -WT[1][2];              WT[2][2] = 0.0;

      // Blockage tensor
      BL[0][0] = block_diag[icv][0];  BL[0][1] = block_offdiag[icv][0];  BL[0][2] = block_offdiag[icv][1];
      BL[1][0] = BL[0][1];            BL[1][1] = block_diag[icv][1];     BL[1][2] = block_offdiag[icv][2];
      BL[2][0] = BL[0][2];            BL[2][1] = BL[1][2];               BL[2][2] = block_diag[icv][2];

      // Strained and rotated A's
      AS[0][0] = as_diag[icv][0];     AS[0][1] = as_offdiag[icv][0];     AS[0][2] = as_offdiag[icv][1];
      AS[1][0] = AS[0][1];            AS[1][1] = as_diag[icv][1];        AS[1][2] = as_offdiag[icv][2];
      AS[2][0] = AS[0][2];            AS[2][1] = AS[1][2];               AS[2][2] = as_diag[icv][2];

      AR[0][0] = ar_diag[icv][0];     AR[0][1] = ar_offdiag[icv][0];     AR[0][2] = ar_offdiag[icv][1];
      AR[1][0] = AR[0][1];            AR[1][1] = ar_diag[icv][1];        AR[1][2] = ar_offdiag[icv][2];
      AR[2][0] = AR[0][2];            AR[2][1] = AR[1][2];               AR[2][2] = ar_diag[icv][2];

      // Call the ASBM
      ierr = 0;

      asbm_(&(ST[0][0]), &(WT[0][0]), &(WFT[0][0]), &(BL[0][0]), &(AS[0][0]), &(AR[0][0]),
            &(REY[0][0]), &(DIM[0][0]), &(CIR[0][0]), &maxitrs, &block_rij, &ierr);

      if (ierr != 0)
      {
        myNerr ++;
        cout << "ASBM error number " << ierr << " in cell " << icv << endl;
        cout << "Cell-x: " << x_cv[icv][0]
             << " Cell-y: " << x_cv[icv][1]
             << " Cell-z: " << x_cv[icv][2] << endl;
      }

      /*//REY[0][0] = rij_diag[icv][0]/(-2.0*rho[icv]*kine[icv]);
      REY[1][1] = rij_diag[icv][1]/(-2.0*rho[icv]*kine[icv]);
      REY[2][2] = rij_diag[icv][2]/(-2.0*rho[icv]*kine[icv]);
      //REY[0][1] = rij_offdiag[icv][0]/(-2.0*rho[icv]*kine[icv]);
      //REY[0][2] = rij_offdiag[icv][1]/(-2.0*rho[icv]*kine[icv]);
      REY[1][2] = rij_offdiag[icv][2]/(-2.0*rho[icv]*kine[icv]);

      rij_diag[icv][0] = -REY[0][0]*2.0*kine[icv]*rho[icv];
      //rij_diag[icv][1] = -REY[1][1]*2.0*kine[icv]*rho[icv];
      //rij_diag[icv][2] = -REY[2][2]*2.0*kine[icv]*rho[icv];

      rij_offdiag[icv][0] = -REY[0][1]*2.0*kine[icv]*rho[icv];
      rij_offdiag[icv][1] = -REY[0][2]*2.0*kine[icv]*rho[icv];
      //rij_offdiag[icv][2] = -REY[1][2]*2.0*kine[icv]*rho[icv];*/

      //REY[0][0] = 1.0/3.0 - temp_muT[icv]/(rho[icv]*kine[icv])*(grad_u[icv][0][0] - 1.0/3.0*diverg[icv]);
      //REY[1][1] = 1.0/3.0 - temp_muT[icv]/(rho[icv]*kine[icv])*(grad_u[icv][1][1] - 1.0/3.0*diverg[icv]);
      //REY[2][2] = 1.0/3.0 - temp_muT[icv]/(rho[icv]*kine[icv])*(grad_u[icv][2][2] - 1.0/3.0*diverg[icv]);
      //REY[0][1] = -temp_muT[icv]/(rho[icv]*kine[icv])*0.5*(grad_u[icv][0][1] + grad_u[icv][1][0]);
      //REY[0][2] = -temp_muT[icv]/(rho[icv]*kine[icv])*0.5*(grad_u[icv][0][2] + grad_u[icv][2][0]);
      //REY[1][2] = -temp_muT[icv]/(rho[icv]*kine[icv])*0.5*(grad_u[icv][1][2] + grad_u[icv][2][1]);

      rij_diag[icv][0] = -REY[0][0]*2.0*kine[icv]*rho[icv];
      rij_diag[icv][1] = -REY[1][1]*2.0*kine[icv]*rho[icv];
      rij_diag[icv][2] = -REY[2][2]*2.0*kine[icv]*rho[icv];

      rij_offdiag[icv][0] = -REY[0][1]*2.0*kine[icv]*rho[icv];
      rij_offdiag[icv][1] = -REY[0][2]*2.0*kine[icv]*rho[icv];
      rij_offdiag[icv][2] = -REY[1][2]*2.0*kine[icv]*rho[icv];

      as_diag[icv][0] = AS[0][0];       as_offdiag[icv][0] = AS[0][1];
      as_diag[icv][1] = AS[1][1];       as_offdiag[icv][1] = AS[0][2];
      as_diag[icv][2] = AS[2][2];       as_offdiag[icv][2] = AS[1][2];

      ar_diag[icv][0] = AR[0][0];       ar_offdiag[icv][0] = AR[0][1];
      ar_diag[icv][1] = AR[1][1];       ar_offdiag[icv][1] = AR[0][2];
      ar_diag[icv][2] = AR[2][2];       ar_offdiag[icv][2] = AR[1][2];

      // debugging
      debug1[icv][0] = DIM[0][0];
      debug1[icv][1] = DIM[1][0];
      debug1[icv][2] = DIM[2][0];

      debug2[icv][0] = DIM[0][1];
      debug2[icv][1] = DIM[1][1];
      debug2[icv][2] = DIM[2][1];

      debug3[icv][0] = DIM[0][2];
      debug3[icv][1] = DIM[1][2];
      debug3[icv][2] = DIM[2][2];

      trace_aa[icv] = CIR[0][0];
      r_ratio[icv] = CIR[1][0];
      etar[icv] = CIR[2][0];

      if (rij_offdiag[icv][0] != rij_offdiag[icv][0])
        marker[icv] = 1;

      rij_diag_nd[icv][0] = REY[0][0];
      rij_diag_nd[icv][1] = REY[1][1];
      rij_diag_nd[icv][2] = REY[2][2];

      rij_offdiag_nd[icv][0] = REY[0][1];
      rij_offdiag_nd[icv][1] = REY[0][2];
      rij_offdiag_nd[icv][2] = REY[1][2];
    }

    updateCvData(rij_diag, REPLACE_ROTATE_DATA);
    updateCvData(rij_offdiag, REPLACE_ROTATE_DATA);

    // Smooth Outputs
    for (int i = 0; i < n_smooth; i++)
    {
      smoothingVec(rij_diag);
      smoothingVec(rij_offdiag);

      updateCvData(rij_diag, REPLACE_ROTATE_DATA);
      updateCvData(rij_offdiag, REPLACE_ROTATE_DATA);
    }

    // Interpolate cell centered stresses to cell faces
    interpolateReStressToFace();

    // Write nan's to screen
    MPI_Barrier(mpi_comm);
    if (mpi_rank != 0)
    {
      int dummy;
      MPI_Status status;
      MPI_Recv(&dummy,1,MPI_INT,mpi_rank-1,1234,mpi_comm,&status);
    }

    for (int icv = 0; icv < ncv; icv++)
      if (marker[icv] == 1)
      {
        writeDiagnostics(icv);
        marker[icv] = 0;
      }

    if ( mpi_rank < mpi_size-1 )
    {
      int dummy = 1;
      MPI_Send(&dummy,1,MPI_INT,mpi_rank+1,1234,mpi_comm);
    }
    MPI_Barrier(mpi_comm);

    // ====================================================================
    // Output error file
    // ====================================================================
    MPI_Reduce(&myNerr,&Nerr,1,MPI_INT,MPI_SUM,0,mpi_comm);
    if (mpi_rank == 0)
    {
      if (finfo == NULL)
      {
        cerr << "Error: cannot write to file asbmInfo.dat" << endl;
        throw(-1);
      }
      else
      {
        fprintf(finfo,"%8d\t%8d\n",step,Nerr);
      }
    }
  }

  void calcBlockTensor()
  {
    // =========================================================================================
    // compute A_block and rhs_block
    // =========================================================================================
    char scalName[] = "bphi";

    double *A_block   = new double[nbocv_s];
    double *rhs_block = new double[ncv];
    double *buffer    = new double[ncv_g];

    for (int icv = 0; icv < ncv; ++icv)     rhs_block[icv] = 0.0;
    for (int noc = 0; noc < nbocv_s; noc++) A_block[noc]   = 0.0;
    for (int icv = 0; icv < ncv_g; ++icv)   buffer[icv]    = 0.0;

    // ..........................................................................................
    // compute the viscous terms
    // ..........................................................................................
    calcViscousFluxScalarAux(rhs_block, A_block, bphi_bfa, grad_bphi);

    // ..........................................................................................
    // compute the source terms
    // ..........................................................................................
    for (int icv = 0; icv < ncv; icv++)
    {
      int noc00 = nbocv_i[icv];
      A_block[noc00] += cv_volume[icv]/(turbLS[icv]*turbLS[icv]);
    }

    // ..........................................................................................
    // under-relaxation
    // ..........................................................................................
    //for (int icv = 0; icv < ncv; icv++){
    //    int noc = nbocv_i[icv];
    //    A_block[noc] /= 1.0;
    //    rhs_block[icv] += (1.0 - 1.0)*A_block[noc]*phi_asbm[icv];
    //}

    // =========================================================================================
    // solve the linear system
    // =========================================================================================
    // //norm of A matrix
    //myAnorm = 0.0; Anorm = 0.0;
    //for (int noc = 0; noc < nbocv_s; noc++)
    //  myAnorm += A_block[noc];
    //MPI_Reduce(&myAnorm,&Anorm,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    //if (mpi_rank == 0){
    //    if (asbmInfo.is_open())
    //      asbmInfo << "Anorm after: " << Anorm << endl;
    //    else
    //      cout << "Unable to write to asbmInfo.dat file." << endl;
    //}

    solveLinSysScalar(buffer, A_block, rhs_block, block_zeroAbs, block_zeroRel, block_maxIter, scalName);

    // clipping
    for (int icv = 0; icv < ncv; ++icv)
      bphi[icv] = min(max(buffer[icv],0.0),1.0);

    updateCvData(bphi, REPLACE_DATA);
    
    // compute residual
    double thisRes, myRes = 0.0, Res = 0.0;
    for (int icv = 0; icv < ncv; icv++)
    {
      thisRes = rhs_block[icv];

      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv+1] - 1;

      thisRes -= A_block[noc_f]*bphi[icv];
      for (int noc = noc_f + 1; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_v[noc];
        thisRes -= A_block[noc]*bphi[icv_nbr];
      }

      myRes += fabs(thisRes);
    }
    MPI_Reduce(&myRes,&Res,1,MPI_DOUBLE,MPI_SUM,0,mpi_comm);
    
    // output for info
    if (mpi_rank == 0)
    {
      if (finfo == NULL)
      {
        cerr << "Error: cannot write to file asbmInfo.dat" << endl;
        throw(-1);
      }
      else
      {
        fprintf(finfo,"bphi resid: %12.6e\n",Res);
      }
    }
    
    delete [] A_block;
    delete [] rhs_block;
    delete [] buffer;

    // =========================================================================================
    // set the boundary conditions
    // =========================================================================================
    for (list<FaZone>::iterator zone = faZoneList.begin(); zone != faZoneList.end(); zone++)
    {
      if (zone->getKind() == FA_ZONE_BOUNDARY)
      {
        Param *param;
        if (getParam(param, zone->getName()))
        {
          if (param->getString() == "WALL")
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              bphi_bfa[ifa] = 1.0;
            }
          }
          else
          {
            for (int ifa = zone->ifa_f; ifa <= zone->ifa_l; ifa++)
            {
              int icv0 = cvofa[ifa][0];
              bphi_bfa[ifa] = bphi[icv0];
            }
          }
        }
      }
    }

    // =========================================================================================
    // compute gradients, with boundary values
    // =========================================================================================
    calcCvScalarGrad(grad_bphi, bphi, bphi_bfa, gradreconstruction, limiterNavierS, bphi, epsilonSDWLS);

    // =========================================================================================
    // compute the blockage tensor
    // =========================================================================================
    double div_bphi;
    for (int icv = 0; icv < ncv; icv++)
    {
      div_bphi = 0.0;
      div_bphi += grad_bphi[icv][0]*grad_bphi[icv][0] +
                  grad_bphi[icv][1]*grad_bphi[icv][1] +
                  grad_bphi[icv][2]*grad_bphi[icv][2];

      if (div_bphi > 1.0e-5){
        block_diag[icv][0] = grad_bphi[icv][0]*grad_bphi[icv][0]/div_bphi*bphi[icv];
        block_diag[icv][1] = grad_bphi[icv][1]*grad_bphi[icv][1]/div_bphi*bphi[icv];
        block_diag[icv][2] = grad_bphi[icv][2]*grad_bphi[icv][2]/div_bphi*bphi[icv];

        block_offdiag[icv][0] = grad_bphi[icv][0]*grad_bphi[icv][1]/div_bphi*bphi[icv];
        block_offdiag[icv][1] = grad_bphi[icv][0]*grad_bphi[icv][2]/div_bphi*bphi[icv];
        block_offdiag[icv][2] = grad_bphi[icv][1]*grad_bphi[icv][2]/div_bphi*bphi[icv];
      }
      else {
        block_diag[icv][0] = 0.0;
        block_diag[icv][1] = 0.0;
        block_diag[icv][2] = 0.0;

        block_offdiag[icv][0] = 0.0;
        block_offdiag[icv][1] = 0.0;
        block_offdiag[icv][2] = 0.0;
      }
    }

    updateCvData(block_diag, REPLACE_DATA);
    updateCvData(block_offdiag, REPLACE_DATA);

  }

  void smoothingVec(double (*input)[3])
  {
    double vol_tot, vol_nbr;
    for (int icv = 0; icv < ncv; icv++)
    {
      vol_tot = cv_volume[icv];

      input[icv][0] *= vol_tot;
      input[icv][1] *= vol_tot;
      input[icv][2] *= vol_tot;

      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;

      for (int noc = noc_f + 1; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_v[noc];

        vol_nbr = cv_volume[icv_nbr];
        vol_tot += vol_nbr;

        input[icv][0] += input[icv_nbr][0]*vol_nbr;
        input[icv][1] += input[icv_nbr][1]*vol_nbr;
        input[icv][2] += input[icv_nbr][2]*vol_nbr;
      }

      input[icv][0] /= vol_tot;
      input[icv][1] /= vol_tot;
      input[icv][2] /= vol_tot;
    }
  }

  void smoothingScal(double *input)
  {
    double vol_tot, vol_nbr;
    for (int icv = 0; icv < ncv; icv++)
    {
      vol_tot = cv_volume[icv];

      input[icv] *= vol_tot;

      int noc_f = nbocv_i[icv];
      int noc_l = nbocv_i[icv + 1] - 1;

      for (int noc = noc_f + 1; noc <= noc_l; noc++)
      {
        int icv_nbr = nbocv_v[noc];

        vol_nbr = cv_volume[icv_nbr];
        vol_tot += vol_nbr;

        input[icv] += input[icv_nbr]*vol_nbr;
      }
      input[icv] /= vol_tot;
    }
  }

  void writeDiagnostics(int icv)
  {
    cout << "Writing Diagnostics"<<endl;
    cout << "step: " << step << endl;
    cout << "icv: " << icv << endl;

    cout << "st_diag: " << st_diag[icv][0] << endl;
    cout << "st_diag: " << st_diag[icv][1] << endl;
    cout << "st_diag: " << st_diag[icv][2] << endl;
    cout << "st_offdiag: " << st_offdiag[icv][0] << endl;
    cout << "st_offdiag: " << st_offdiag[icv][1] << endl;
    cout << "st_offdiag: " << st_offdiag[icv][2] << endl;

    cout << "wt_offdiag: " << wt_offdiag[icv][0] << endl;
    cout << "wt_offdiag: " << wt_offdiag[icv][1] << endl;
    cout << "wt_offdiag: " << wt_offdiag[icv][2] << endl;

    cout << "as_diag: " << as_diag[icv][0] << endl;
    cout << "as_diag: " << as_diag[icv][1] << endl;
    cout << "as_diag: " << as_diag[icv][2] << endl;
    cout << "as_offdiag: " << as_offdiag[icv][0] << endl;
    cout << "as_offdiag: " << as_offdiag[icv][1] << endl;
    cout << "as_offdiag: " << as_offdiag[icv][2] << endl;

    cout << "ar_diag: " << ar_diag[icv][0] << endl;
    cout << "ar_diag: " << ar_diag[icv][1] << endl;
    cout << "ar_diag: " << ar_diag[icv][2] << endl;
    cout << "ar_offdiag: " << ar_offdiag[icv][0] << endl;
    cout << "ar_offdiag: " << ar_offdiag[icv][1] << endl;
    cout << "ar_offdiag: " << ar_offdiag[icv][2] << endl;

    cout << "block_diag: " << block_diag[icv][0] << endl;
    cout << "block_diag: " << block_diag[icv][1] << endl;
    cout << "block_diag: " << block_diag[icv][2] << endl;
    cout << "block_offdiag: " << block_offdiag[icv][0] << endl;
    cout << "block_offdiag: " << block_offdiag[icv][1] << endl;
    cout << "block_offdiag: " << block_offdiag[icv][2] << endl;

    cout << "bphi: " << bphi[icv] << endl;

    cout << "rij_diag: " << rij_diag[icv][0] << endl;
    cout << "rij_diag: " << rij_diag[icv][1] << endl;
    cout << "rij_diag: " << rij_diag[icv][2] << endl;
    cout << "rij_offdiag: " << rij_offdiag[icv][0] << endl;
    cout << "rij_offdiag: " << rij_offdiag[icv][1] << endl;
    cout << "rij_offdiag: " << rij_offdiag[icv][2] << endl;

    /*cout << "turbTS: " << turbTS[icv] << endl;
    cout << "tturb: " << tturb[icv] << endl;
    cout << "tkol: " << tkol[icv] << endl;
    cout << "trel: " << trel[icv] << endl;
    cout << "kine: " << kine[icv] << endl;
    cout << "eps:  " << eps[icv] << endl;
    cout << "v2  : " << v2[icv] << endl;*/
    cout << "str : " << strMag[icv] << endl << endl;
  }

  virtual void finalHookScalarRansTurbModel()
  {
  	if (marker    != NULL) {delete [] marker;       marker    = NULL;}
    if (bphi_bfa  != NULL) {delete [] bphi_bfa;     bphi_bfa  = NULL;}
    if (grad_bphi != NULL) {delete [] grad_bphi;    grad_bphi = NULL;}
    if (mpi_rank == 0) fclose(finfo);
  }

};

//######################################################//
//                                                      //
// ASBM with keps Model                                  //
//                                                      //
//######################################################//

class RansTurbASBMkeps : virtual public RansTurbASBM
{
public:
  RansTurbASBMkeps()
  {
    if (mpi_rank == 0)
      cout << "RansTurbASBMkeps()" << endl;

    turbModel = ASBMkeps;

    C_MU  = getDoubleParam("C_MU",  "0.22");
    SIG_K = getDoubleParam("SIG_K", "1.0");
    SIG_D = getDoubleParam("SIG_D", "1.3");
    CEPS1 = getDoubleParam("CEPS1", "1.4");
    CEPS2 = getDoubleParam("CEPS2", "1.9");
    CETA  = getDoubleParam("CETA",  "70.0");
    CL    = getDoubleParam("CL",    "0.23");

    LIMIT_PK     = getIntParam("LIMIT_PK",     "0");
    LIMIT_TL     = getIntParam("LIMIT_TL",     "1");
    VEL_SCALE    = getIntParam("VEL_SCALE",    "0");
    RIJ_BASED_PK = getIntParam("RIJ_BASED_PK", "0");

    ScalarTranspEq *eq;
    eq = registerScalarTransport("kine",  CV_DATA);
    eq->relax = getDoubleParam("RELAX_kine", "0.7");
    eq->phiZero = 1.0e-10;
    eq->phiMaxiter = 500;
    eq->lowerBound = 1.0e-10;
    eq->upperBound = 1.0e10;
    eq->turbSchmidtNumber = SIG_K;

    eq = registerScalarTransport("eps", CV_DATA);
    eq->relax = getDoubleParam("RELAX_eps", "0.7");
    eq->phiZero = 1.0e-8;
    eq->phiMaxiter = 500;
    eq->lowerBound = 1.0e-4;
    eq->upperBound = 1.0e10;
    eq->turbSchmidtNumber = SIG_D;

    omega   = NULL;   registerScalar(omega,   "omega"  , CV_DATA);
    v2      = NULL;   registerScalar(v2,      "v2"     , CV_DATA);
    strMag  = NULL;   registerScalar(strMag,  "strMag" , CV_DATA);
    diverg  = NULL;   registerScalar(diverg,  "diverg" , CV_DATA);
    turbTS  = NULL;   registerScalar(turbTS,  "turbTS" , CV_DATA);
    turbLS  = NULL;   registerScalar(turbLS,  "turbLS" , CV_DATA);
    muT     = NULL;   registerScalar(muT,     "muT"    , CV_DATA);

    tturb  = NULL;       registerScalar(tturb, "tturb", CV_DATA);
    tkol   = NULL;       registerScalar(tkol, "tkol", CV_DATA);
    trel   = NULL;       registerScalar(trel, "trel", CV_DATA);
    lturb  = NULL;       registerScalar(lturb, "lturb", CV_DATA);
    lkol   = NULL;       registerScalar(lkol, "lkol", CV_DATA);
    lrel   = NULL;       registerScalar(lrel, "lrel", CV_DATA);
  }

  virtual ~RansTurbASBMkeps() {}

public:

  double *eps, *omega, *v2;   // turb scalars
  double *kine_bfa, *eps_bfa; // turb scalars at boundary
  double *muT;                // turb visc at cell center

  double C_MU, SIG_K, SIG_D, CEPS1, CEPS2, CETA, CL; // model constants

  double *tturb, *tkol, *trel, *lkol, *lrel, *lturb;

  int LIMIT_PK, LIMIT_TL, VEL_SCALE, RIJ_BASED_PK;

public:
  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0)
      cout << "initialHook() for ASBMkeps" << endl;

    ScalarTranspEq *eq;
    eq = getScalarTransportData("kine");
    kine = eq->phi;
    kine_bfa = eq->phi_bfa;

    eq = getScalarTransportData("eps");
    eps = eq->phi;
    eps_bfa = eq->phi_bfa;

    for (int ifa = 0; ifa < nfa; ifa++)
    	nonLinear[ifa] = 1.0;
    if (mpi_rank == 0)
    	cout << "NON-LINEAR DOMAIN ACTIVATED." << endl;

    RansTurbASBM::initialHookScalarRansTurbModel();

    // Output to screen
    if (mpi_rank == 0)
    {
      cout << "ASBM properties" << endl;
      switch(LIMIT_PK)
      {
      case 0: cout << "    PK limiter: no" << endl; break;
      case 1: cout << "    PK limiter: yes" << endl; break;
      }
      switch(LIMIT_TL)
      {
      case 0: cout << "    TL limiter: no" << endl; break;
      case 1: cout << "    TL limiter: yes" << endl; break;
      }
      switch(RIJ_BASED_PK)
      {
      case 0: cout << "    PK based on ASBM stresses: no" << endl; break;
      case 1: cout << "    Pk based on ASBM stresses: yes" << endl; break;
      }
      switch(VEL_SCALE)
      {
      case 0:  cout<<"    Velocity scale = VY          "<<endl; break;
      case 1:  cout<<"    Velocity scale = VNXY        "<<endl; break;
      case 11: cout<<"    Velocity scale = VNXYXY      "<<endl; break;
      case 2:  cout<<"    Velocity scale = VNXZ        "<<endl; break;
      case 3:  cout<<"    Velocity scale = 0.66 TKE    "<<endl; break;
      default: cout<<"    Velocity scale = v2 from v2f "<<endl; break;
      }
      switch(block_rij)
      {
      case 0: cout <<"    Blocking rij: no" << endl; break;
      case 1: cout <<"    Blocking rij: yes" << endl; break;
      }
    }
  }

  virtual void calcRansTurbViscMuet()
  {
  /*if (step == 50000){
      for (int icv = 0; icv < ncv; icv++)
        eps[icv] = kine[icv]*omega[icv];
      updateCvData(eps, REPLACE_DATA);
    }*/

    // update velocity gradients, timescale, lengthscale
    calcGradVel();
    calcStrainRateAndDivergence();
    calcTurbTimeScale();
    calcTurbLengthScale();

    // compute wall blocking tensor
    if ( step%block_frq == 0 )
      calcBlockTensor();

    // compute Reynolds stresses
    calcRsCenterASBM();

    // compute wall-normal stress
    if(VEL_SCALE==0)
      for (int icv = 0; icv < ncv; icv++)
      {
        double fluc_norm = fabs(rij_diag[icv][1]);
        v2[icv] = min(fluc_norm/rho[icv],2./3.*kine[icv]);
      }

    if(VEL_SCALE==1)
      for (int icv = 0; icv < ncv; icv++)
      {
        double nx = vel[icv][0]/sqrt(vel[icv][0]*vel[icv][0]+vel[icv][1]*vel[icv][1]+1.e-12);
        double ny = vel[icv][1]/sqrt(vel[icv][0]*vel[icv][0]+vel[icv][1]*vel[icv][1]+1.e-12);
        double fluc_norm = fabs(rij_diag[icv][0]*ny-rij_diag[icv][1]*nx);
        v2[icv] = min(fluc_norm/rho[icv],2./3.*kine[icv]);
      }

    if(VEL_SCALE==11)
      for (int icv = 0; icv < ncv; icv++)
      {
        double nx = vel[icv][0]/sqrt(vel[icv][0]*vel[icv][0]+vel[icv][1]*vel[icv][1]+1.e-12);
        double ny = vel[icv][1]/sqrt(vel[icv][0]*vel[icv][0]+vel[icv][1]*vel[icv][1]+1.e-12);
        double fluc_norm = fabs(rij_diag[icv][0]*ny*ny+rij_diag[icv][1]*nx*nx-2.0*nx*ny*rij_offdiag[icv][0]);
        v2[icv] = min(fluc_norm/rho[icv],2./3.*kine[icv]);
      }

    if(VEL_SCALE==2)
      for (int icv = 0; icv < ncv; icv++)
      {
        double nx = vel[icv][0]/sqrt(vel[icv][0]*vel[icv][0]+vel[icv][2]*vel[icv][2]+1.e-12);
        double nz = vel[icv][2]/sqrt(vel[icv][0]*vel[icv][0]+vel[icv][2]*vel[icv][2]+1.e-12);
        double fluc_norm = fabs(rij_diag[icv][0]*nz-rij_diag[icv][2]*nx);
        //double fluc_norm = fabs(rij_diag[icv][0]*nz*nz+rij_diag[icv][2]*nx*nx-2*nx*nz*rij_offdiag[icv][1]);
        v2[icv] = min(fluc_norm/rho[icv],2./3.*kine[icv]);
      }

    if(VEL_SCALE==3)
      for (int icv = 0; icv < ncv; icv++)
      {
        v2[icv] = 2./3.*kine[icv];
      }

    updateCvData(v2, REPLACE_DATA);

    // calculate turbulent viscosity
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
      double turbTS_fa = (w1*turbTS[icv0] + w0*turbTS[icv1])*invwtot;
      double v2_fa     = (w1*v2[icv0] + w0*v2[icv1])*invwtot;

      mut_fa[ifa] = min(max(C_MU*rho_fa*v2_fa*turbTS_fa, 0.0), 1.0);
    }

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
            mut_fa[ifa] = min(max(C_MU*rho[icv0]*v2[icv0]*turbTS[icv0], 0.0), 1.0);
          }
      }
    }

    for (int icv=0; icv<ncv; icv++)
      muT[icv] = InterpolateAtCellCenterFromFaceValues(mut_fa, icv);
  }

  /*
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
    if (name == "eps")
      for (int icv = 0; icv < ncv; icv++)
      {
        double ce1 = CEPS1*(1.0 + 0.0425*pow(kine[icv]/v2[icv], 0.5));

        double src = (ce1*getTurbProd(icv) - CEPS2*rho[icv]*eps[icv])/turbTS[icv];
        rhs[icv]  += src*cv_volume[icv];

        // d(ce2*rho*eps/TS)/d(rho*eps)
        if (flagImplicit)
        {
          double dsrcdphi = -CEPS2/turbTS[icv];
          int noc00 = nbocv_i[icv];
          A[noc00] -= dsrcdphi*cv_volume[icv];
        }
      }
  }

  virtual void sourceHookRansTurbCoupled(double **rhs, double ***A, int nScal, int flagImplicit)
  {
    int kine_Index = 5+getScalarTransportIndex("kine");
    int eps_Index =  5+getScalarTransportIndex("eps");

    double dsrcdphi;

    // mu_t*str*str - rho*k/k*eps
    for (int icv = 0; icv < ncv; icv++)
    {
      double kine_src  = getTurbProd(icv) - rho[icv]*eps[icv];
      double ce1 = CEPS1*(1.0 + 0.045*pow(kine[icv]/v2[icv], 0.5));
      double eps_src = (ce1*getTurbProd(icv) - CEPS2*rho[icv]*eps[icv])/turbTS[icv];

      rhs[icv][kine_Index] += kine_src*cv_volume[icv];
      rhs[icv][eps_Index]  += eps_src*cv_volume[icv];

      if (flagImplicit)
      {
        int noc00 = nbocv_i[icv];

        dsrcdphi = -eps[icv]/kine[icv];
        A[noc00][kine_Index][kine_Index] -= dsrcdphi*cv_volume[icv];
        dsrcdphi = -1.0;
        A[noc00][kine_Index][eps_Index]  -= dsrcdphi*cv_volume[icv];

        dsrcdphi = -CEPS2/turbTS[icv];
        A[noc00][eps_Index][eps_Index] -= dsrcdphi*cv_volume[icv];
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
    for (int icv = 0; icv < ncv; icv++)
    {
      double nu = calcMuLam(icv)/rho[icv];

      double tau    = kine[icv]/eps[icv];
      double tauKol = 6.0*sqrt(nu/eps[icv]);
      double tauRel = 0.6*kine[icv]/(max(sqrt(3.0)*v2[icv]*C_MU*strMag[icv],1.0e-14));

      switch (LIMIT_TL)
      {
      case 0: turbTS[icv] = tau;                         break;
      case 1: turbTS[icv] = min(max(tau,tauKol),tauRel); break;
      }
    }
    updateCvData(turbTS, REPLACE_DATA);
  }

  void calcTurbLengthScale()
  {
    for (int icv = 0; icv < ncv; icv++)
    {
      double nu = calcMuLam(icv)/rho[icv];

      double len    = pow(kine[icv],1.5)/eps[icv];
      double lenKol = CETA*pow(nu,0.75)/pow(eps[icv],0.25);
      double lenRel = pow(kine[icv],1.5)/(max(sqrt(3.0)*v2[icv]*C_MU*strMag[icv],1.0e-14));

      switch (LIMIT_TL)
      {
      case 0: turbLS[icv] = CL*len;                         break;
      case 1: turbLS[icv] = CL*max(min(len,lenRel),lenKol); break;
      }
    }
    updateCvData(turbLS, REPLACE_DATA);
  }

  double getTurbProd(int icv)
  {
    double Pk, mu_t;

    switch(RIJ_BASED_PK)
    {
    case 0:
      mu_t = min(max(C_MU*rho[icv]*v2[icv]*turbTS[icv], 0.0), 1.0);
      Pk = mu_t*strMag[icv]*strMag[icv] - 2.0/3.0*rho[icv]*kine[icv]*diverg[icv];
      break;
    case 1:
      Pk = rij_diag[icv][0]*grad_u[icv][0][0]    + rij_offdiag[icv][0]*grad_u[icv][0][1] + rij_offdiag[icv][1]*grad_u[icv][0][2] +
           rij_offdiag[icv][0]*grad_u[icv][1][0] + rij_diag[icv][1]*grad_u[icv][1][1]    + rij_offdiag[icv][2]*grad_u[icv][1][2] +
           rij_offdiag[icv][1]*grad_u[icv][2][0] + rij_offdiag[icv][2]*grad_u[icv][2][1] + rij_diag[icv][2]*grad_u[icv][2][2];
      break;
    }

    Pk = max(Pk,0.0);
    if (LIMIT_PK == 1)
      Pk = min(Pk, 20.0*rho[icv]*eps[icv]);
    return Pk;
  }

  virtual void finalHookScalarRansTurbModel()
  {
    RansTurbASBM::finalHookScalarRansTurbModel();
  }

};

//######################################################//
//                                                      //
// ASBM with Wilcox k-omega Model                       //
//                                                      //
//######################################################//

class RansTurbASBMkom : virtual public RansTurbASBM
{
public:   // constructors

  RansTurbASBMkom()
  {
    if (mpi_rank == 0)
      cout << "RansTurbASBMkom()" << endl;

    turbModel = ASBMkom;

    betaStar = getDoubleParam("betaStar", "0.09");
    sigma_k  = getDoubleParam("sigma_k",  "0.6");
    sigma_om = getDoubleParam("sigma_om", "0.5");
    alfa     = getDoubleParam("alfa",     "0.52");
    sigmad0  = getDoubleParam("sigmad0",  "0.125");
    beta0    = getDoubleParam("beta0",    "0.0708");
    cLim     = getDoubleParam("cLim",     "0.875");
    CETA     = getDoubleParam("CETA",     "70.0");
    CL       = getDoubleParam("CL",       "0.03");

    RIJ_BASED_PK = getIntParam("RIJ_BASED_PK", "0");
    LIMIT_TL     = getIntParam("LIMIT_TL",     "1");

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
    diverg   = NULL;       registerScalar(diverg, "diverg", CV_DATA);
    turbTS   = NULL;       registerScalar(turbTS, "turbTS", CV_DATA);
    turbLS   = NULL;       registerScalar(turbLS, "turbLS", CV_DATA);
    muT      = NULL;       registerScalar(muT, "muT", CV_DATA);
    wallDist = NULL;       registerScalar(wallDist, "wallDist",  CV_DATA);

    tturb  = NULL;       registerScalar(tturb, "tturb", CV_DATA);
    tkol   = NULL;       registerScalar(tkol, "tkol", CV_DATA);
    trel   = NULL;       registerScalar(trel, "trel", CV_DATA);
    lturb  = NULL;       registerScalar(lturb, "lturb", CV_DATA);
    lkol   = NULL;       registerScalar(lkol, "lkol", CV_DATA);
    lrel   = NULL;       registerScalar(lrel, "lrel", CV_DATA);
  }

  virtual ~RansTurbASBMkom() {}

public:

  // model variables
  double *omega;           ///< specific dissipation
  double (*grad_kine)[3];  ///< gradient of tke
  double (*grad_omega)[3]; ///< gradient of omega
  double *kine_bfa;        ///< tke at the boundaries
  double *omega_bfa;       ///< omega at the boundaries
  double *muT;             ///< turbulent viscosity at cell center
  double *wallDist;        ///< wall distance

  double *tturb, *tkol, *trel, *lkol, *lrel, *lturb;
  int LIMIT_TL, RIJ_BASED_PK;

  // model constants
  double betaStar, sigma_om, sigma_k, alfa, sigmad0, beta0, cLim;
  double CL, CETA;

public:

  virtual void initialHookScalarRansTurbModel()
  {
    if (mpi_rank == 0)
      cout << "initialHook() for ASBMkom" << endl;

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

    for (int ifa = 0; ifa < nfa; ifa++)
      nonLinear[ifa] = 1.0;
    if (mpi_rank == 0)
      cout << "NON-LINEAR DOMAIN ACTIVATED." << endl;

    RansTurbASBM::initialHookScalarRansTurbModel();

    // Output to screen
    if (mpi_rank == 0)
    {
      cout << "ASBM properties" << endl;
      switch(LIMIT_TL)
      {
      case 0: cout << "    TL limiter: no" << endl; break;
      case 1: cout << "    TL limiter: yes" << endl; break;
      }
      switch(RIJ_BASED_PK)
      {
      case 0: cout << "    PK based on ASBM stresses: no" << endl; break;
      case 1: cout << "    Pk based on ASBM stresses: yes" << endl; break;
      }
      switch(block_rij)
      {
      case 0: cout <<"    Blocking rij: no" << endl; break;
      case 1: cout <<"    Blocking rij: yes" << endl; break;
      }
    }

  }

  virtual void calcRansTurbViscMuet()
  {
    // update velocity gradients, timescale, lengthscale
    calcGradVel();
    calcStrainRateAndDivergence();
    calcTurbTimeScale();
    calcTurbLengthScale();

    // eddy-viscosity at cell center
    /*for (int icv = 0; icv < ncv; icv++)
    {
      double omega_tilde = max(omega[icv], cLim*strMag[icv]/sqrt(betaStar));
      temp_muT[icv] = min(rho[icv]*kine[icv]/omega_tilde, 100.0);
    }*/
    
    // compute wall blocking tensor
    if ( step%block_frq == 0 )
      calcBlockTensor();

    // compute Reynolds stresses
    calcRsCenterASBM();

    // eddy-viscosity at cell center
    for (int icv = 0; icv < ncv; icv++)
    {
      double omega_tilde = max(omega[icv], cLim*strMag[icv]/sqrt(betaStar));
      muT[icv] = min(rho[icv]*kine[icv]/omega_tilde, 100.0);
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
            mut_fa[ifa] = 0.0;        // set mut zero at walls
        else
          for (int ifa=zone->ifa_f; ifa<=zone->ifa_l; ifa++)
          {
            int icv0 = cvofa[ifa][0];
            mut_fa[ifa] = muT[icv0];  // zero order extrapolation for others
          }
      }
  }

  /*
   * calculate diffusivity scalars
   */
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
        double Pk = getTurbProd(icv);

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
      double OM[3][3], STR_hat[3][3];

      calcCvScalarGrad(grad_kine,  kine,  kine_bfa,  gradreconstruction, limiterNavierS, kine,  epsilonSDWLS);
      calcCvScalarGrad(grad_omega, omega, omega_bfa, gradreconstruction, limiterNavierS, omega, epsilonSDWLS);

      for (int icv=0; icv<ncv; icv++)
      {
        // Production of tke
        double Pk = getTurbProd(icv);

        // Destruction of omega
        for (int i=0; i<3; i++)
          for (int j=0; j<3; j++)
          {
            OM[i][j] = 0.5*(grad_u[icv][i][j] - grad_u[icv][j][i]);
            // STR_hat subtracts 0.5 divergence instead of 1/3!!!
            if (i==j)   STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i] - diverg[icv]);
            else        STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]);
          }

        double chiOm = 0.0;
        for (int i=0; i<3; i++)
          for (int j=0; j<3; j++)
            for (int k=0; k<3; k++)
              chiOm += OM[i][j]*OM[j][k]*STR_hat[k][i];
        chiOm = fabs(chiOm/pow(betaStar*omega[icv], 3.0));

        double fbeta = (1.0 + 85.0*chiOm)/(1.0 + 100.0*chiOm);
        double beta = beta0*fbeta;

        // Cross-diffusion of omega
        double sigmad;
        double crossDiff = vecDotVec3d(grad_kine[icv], grad_omega[icv]);
        if (crossDiff <= 0.0) sigmad = 0.0;
        else                  sigmad = sigmad0;

        double src =  alfa*omega[icv]/kine[icv]*Pk
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
      // Production of tke
      double Pk = getTurbProd(icv);

      // Destruction of omega
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
        {
          OM[i][j] = 0.5*(grad_u[icv][i][j] - grad_u[icv][j][i]);

          // STR_hat subtracts 0.5 divergence instead 1/3, no bug!!!
          if (i==j)   STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i] - diverg[icv]);
          else        STR_hat[i][j] = 0.5*(grad_u[icv][i][j] + grad_u[icv][j][i]);
        }

      double chiOm = 0.0;
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
          for (int k=0; k<3; k++)
            chiOm += OM[i][j]*OM[j][k]*STR_hat[k][i];
      chiOm = fabs(chiOm/pow(betaStar*omega[icv], 3.0));

      double fbeta = (1.0 + 85.0*chiOm)/(1.0 + 100.0*chiOm);
      double beta = beta0*fbeta;

      // Cross-diffusion of omega
      double sigmad;
      double crossDiff = vecDotVec3d(grad_kine[icv], grad_omega[icv]);
      if (crossDiff <= 0.0) sigmad = 0.0;
      else                  sigmad = sigmad0;

      double src =  alfa*omega[icv]/kine[icv]*Pk
                  - beta*rho[icv]*omega[icv]*omega[icv]
                  + sigmad*rho[icv]/omega[icv]*crossDiff;

      rhs[icv][5+kine_Index]  += (Pk - betaStar*rho[icv]*omega[icv]*kine[icv])*cv_volume[icv];
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
  void calcTurbTimeScale()
  {
    for (int icv = 0; icv < ncv; icv++)
    {
      double nu = calcMuLam(icv)/rho[icv];

      double tau    = 1.0/max(betaStar*omega[icv],betaStar*cLim*strMag[icv]/sqrt(betaStar));
      double tauKol = 6.0*sqrt(nu/(betaStar*kine[icv]*omega[icv]));
      double tauRel = 0.6*sqrt(3.0)/(max(2.0*betaStar*strMag[icv],1.0e-14));

      switch (LIMIT_TL)
      {
      case 0: turbTS[icv] = tau;                         break;
      case 1: turbTS[icv] = min(max(tau,tauKol),tauRel); break;
      }
    }
    updateCvData(turbTS,REPLACE_DATA);
  }

  void calcTurbLengthScale()
  {
    for (int icv = 0; icv < ncv; icv++)
    {
      double nu = calcMuLam(icv)/rho[icv];

      double len    = pow(kine[icv],0.5)/(betaStar*omega[icv]);
      double lenKol = CETA*pow(nu,0.75)/pow(betaStar*kine[icv]*omega[icv],0.25);
      double lenRel = pow(kine[icv],0.5)*sqrt(3.0)/(max(2.0*betaStar*strMag[icv],1.0e-14));

      switch (LIMIT_TL)
      {
      case 0: turbLS[icv] = CL*len;                         break;
      case 1: turbLS[icv] = CL*max(min(len,lenRel),lenKol); break;
      }
    }
    updateCvData(turbLS,REPLACE_DATA);
  }

  double getTurbProd(int icv)
  {
    double Pk;

    switch(RIJ_BASED_PK)
    {
    case 0:
      Pk = muT[icv]*strMag[icv]*strMag[icv] - 2./3.*rho[icv]*kine[icv]*diverg[icv];
      break;
    case 1:
      Pk = rij_diag[icv][0]*grad_u[icv][0][0]    + rij_offdiag[icv][0]*grad_u[icv][0][1] + rij_offdiag[icv][1]*grad_u[icv][0][2] +
           rij_offdiag[icv][0]*grad_u[icv][1][0] + rij_diag[icv][1]*grad_u[icv][1][1]    + rij_offdiag[icv][2]*grad_u[icv][1][2] +
           rij_offdiag[icv][1]*grad_u[icv][2][0] + rij_offdiag[icv][2]*grad_u[icv][2][1] + rij_diag[icv][2]*grad_u[icv][2][2];
      Pk = min(Pk, 20.0*betaStar*rho[icv]*omega[icv]*kine[icv]);
      break;
    }

    Pk = max(Pk, 0.0);
    return Pk;
  }

  virtual void finalHookScalarRansTurbModel()
  {
    RansTurbASBM::finalHookScalarRansTurbModel();
  }

};

#endif
