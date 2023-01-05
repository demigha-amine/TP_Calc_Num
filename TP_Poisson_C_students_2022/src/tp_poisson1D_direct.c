/******************************************/
/* tp2_poisson1D_direct.c                 */
/* This file contains the main function   */
/* to solve the Poisson 1D problem        */
/******************************************/
#include "lib_poisson1D.h"
#include "atlas_headers.h"
#include <time.h>
int main(int argc,char *argv[])
/* ** argc: Nombre d'arguments */
/* ** argv: Valeur des arguments */
{
  int ierr;
  int jj;
  int nbpoints, la;
  int ku, kl, kv, lab;
  int *ipiv;
  int info;
  int NRHS;
  double T0, T1;
  double *RHS, *EX_SOL, *X;
  double **AAB;
  double *AB;

  double temp, relres;

  NRHS=1;
  nbpoints=10;
  la=nbpoints-2;
  T0=-5.0;
  T1=5.0;

  printf("--------- Poisson 1D ---------\n\n");
  RHS=(double *) malloc(sizeof(double)*la);
  EX_SOL=(double *) malloc(sizeof(double)*la);
  X=(double *) malloc(sizeof(double)*la);

  // TODO : you have to implement those functions
  set_grid_points_1D(X, &la);
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
  set_analytical_solution_DBC_1D(EX_SOL, X, &la, &T0, &T1);
  
  write_vec(RHS, &la, "RHS.dat");
  write_vec(EX_SOL, &la, "EX_SOL.dat");
  write_vec(X, &la, "X_grid.dat");

  kv=1;
  ku=1;
  kl=1;
  lab=kv+kl+ku+1;

  AB = (double *) malloc(sizeof(double)*lab*la);

  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB.dat");

  // set_GB_operator_colMajor_poisson1D_Id(AB, &lab, &la, &kv);
  // write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "AB_Id.dat");

  double *RHS_dgbmv = (double *) malloc(sizeof(double)*la);

  cblas_dgbmv(CblasColMajor,CblasNoTrans,la,la,kl,ku,1.0,AB+1,lab,EX_SOL,1,0.0,RHS_dgbmv,1);
  write_vec(RHS_dgbmv, &la, "RHS_dgbmv.dat");
  
  printf("methode de validation \n");
  cblas_daxpy(la, -1, RHS_dgbmv, 1, RHS, 1);

  //calculer la norme
  double norm = cblas_dnrm2(la, RHS, 1);
  printf("Norm = %f\n",norm);
  printf("La methode est validé\n");


  printf("Solution with LAPACK\n");
  /* LU Factorization */
  info=0;
  ipiv = (int *) calloc(la, sizeof(int));

   /* LU for tridiagonal matrix  (can replace dgbtrf_) */
  printf("Factorisation LU sans DGBTRF \n");

  clock_t debut_lu = clock();
  ierr = dgbtrftridiag(&la, &la, &kl, &ku, AB, &lab, ipiv, &info, &kv);
  clock_t fin_lu = clock();

  printf("Temps d'execution LU sans DGBTRF  = %f seconds\n", (double)(fin_lu - debut_lu) / CLOCKS_PER_SEC);


  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "notre_LU.dat");
  set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);

  
  clock_t debut_trf = clock();

  //dgbtrf_(&la, &la, &kl, &ku, AB, &lab, ipiv, &info);
  info = LAPACKE_dgbtrf(LAPACK_COL_MAJOR, la, la, kl, ku, AB, lab, ipiv);

  clock_t fin_trf = clock();
  printf("Temps d'execution DGBTRF  = %f seconds\n", (double)(fin_trf - debut_trf) / CLOCKS_PER_SEC);

  
  write_GB_operator_colMajor_poisson1D(AB, &lab, &la, "LU.dat");
 
  
  /* Solution (Triangular) */
  set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);

  if (info==0){

    clock_t debut_trs = clock();

    //dgbtrs_("N", &la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);
    info = LAPACKE_dgbtrs(LAPACK_COL_MAJOR, 'N', la, kl, ku, NRHS, AB, lab, ipiv, RHS, la);

    clock_t fin_trs = clock();

    write_xy(RHS, X, &la, "TRS.dat");

    printf("Temps d'execution dgbtrf + dgbtrs = %f seconds\n", (double)((fin_trf - debut_trf) + (fin_trs - debut_trs)) / CLOCKS_PER_SEC);

    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
  }else{
    printf("\n INFO = %d\n",info);
  }

  /* It can also be solved with dgbsv */
  // TODO : use dgbsv
  /* It can also be solved with dgbsv (dgbtrf+dgbtrs) */
  /*  ierr = dgbsv_(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);*/
  
    set_dense_RHS_DBC_1D(RHS,&la,&T0,&T1);
    set_GB_operator_colMajor_poisson1D(AB, &lab, &la, &kv);


    if (info==0){
    clock_t debut_sv = clock();

    LAPACK_dgbsv(&la, &kl, &ku, &NRHS, AB, &lab, ipiv, RHS, &la, &info);

    clock_t fin_sv = clock();
      write_xy(RHS, X, &la, "SV.dat");

    printf("Temps d'execution dgbsv = %f seconds\n",(double) (fin_sv - debut_sv) / CLOCKS_PER_SEC);


    if (info!=0){printf("\n INFO DGBTRS = %d\n",info);}
    }else{
     printf("\n INFO = %d\n",info);
    }




  write_xy(RHS, X, &la, "SOL.dat");

  /* Relative forward error */
  temp = cblas_ddot(la, RHS, 1, RHS,1);
  temp = sqrt(temp);
  cblas_daxpy(la, -1.0, RHS, 1, EX_SOL, 1);
  relres = cblas_ddot(la, EX_SOL, 1, EX_SOL,1);
  relres = sqrt(relres);
  relres = relres / temp;
  
  printf("\nThe relative forward error is relres = %e\n",relres);

  free(RHS);
  free(EX_SOL);
  free(X);
  free(AB);
  printf("\n\n--------- End -----------\n");

}
