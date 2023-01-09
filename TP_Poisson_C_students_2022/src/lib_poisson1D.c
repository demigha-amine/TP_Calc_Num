/**********************************************/
/* lib_poisson1D.c                            */
/* Numerical library developed to solve 1D    */ 
/* Poisson problem (Heat equation)            */
/**********************************************/
#include "lib_poisson1D.h"

void set_GB_operator_colMajor_poisson1D(double* AB, int *lab, int *la, int *kv)
{

  for (int j=0; j<(*la); j++)
  {
    int indice  = j * (*lab); //l'indice j*n

    //pour la premiere colonne = 0 (kv = 1)
    if (*kv>=0)
    {
      for (int i=0; i< *kv; i++)  AB[indice+i]=0.0;
    }
    
    //pour les autres elements
    AB[indice+ *kv]=-1.0;
    AB[indice+ *kv+1]=2.0;
    AB[indice+ *kv+2]=-1.0;
  }
  //pour le premier element
  AB[0]=0.0;
  if(*kv == 1) AB[1]=0.0;
  
  //pour le dernier element = 0
  AB[(*lab)*(*la)-1]=0.0; //indice = n*m - 1 
}

void set_GB_operator_colMajor_poisson1D_Id(double* AB, int *lab, int *la, int *kv)
{

  for (int j=0; j<(*la); j++)
  {
    int indice  = j * (*lab); //l'indice j*n

    //pour la premiere colonne = 0 (kv = 1)
    if (*kv>=0)
    {
      for (int i=0; i< *kv; i++)  AB[indice+i]=0.0;
    }
    
    //pour les autres elements
    AB[indice+ *kv]=0.0;
    AB[indice+ *kv+1]=1.0;
    AB[indice+ *kv+2]=0.0;
  }
  //pour le premier element apres la colonne kv
  AB[1]=0.0;

  //pour le dernier element = 0
  AB[(*lab)*(*la)-1]=0.0; //indice = n*m - 1 
}

void set_dense_RHS_DBC_1D(double* RHS, int* la, double* BC0, double* BC1){
  //premier element = T0
  RHS[0] = *BC0;
  //dernier element = T1
  RHS[*la-1] = *BC1;

  //les autres elements = 0
  for (int i=1; i<(*la-1); i++)
    RHS[i] = 0.0;

}  

void set_analytical_solution_DBC_1D(double* EX_SOL, double* X, int* la, double* BC0, double* BC1){
  //diff entre T1 - T0
  double diff = *BC1 - *BC0;

  //T (x) = T 0 + x(T 1 - T 0 )
  for(int i=0; i<*la; i++)
    EX_SOL[i] = *BC0 + X[i] * diff ; 

}  

void set_grid_points_1D(double* x, int* la){
  //le pas h
  double h = 1.0 / (*la + 1);

  //from ]0,1[
  for (int i=0; i<(*la); i++)
    x[i] = (i+1) * h;
}

void write_GB_operator_rowMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*lab);ii++){
      for (jj=0;jj<(*la);jj++){
	fprintf(file,"%lf\t",AB[ii*(*la)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB_operator_colMajor_poisson1D(double* AB, int* lab, int* la, char* filename){
  FILE * file;
  int ii,jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (ii=0;ii<(*la);ii++){
      for (jj=0;jj<(*lab);jj++){
	fprintf(file,"%lf\t",AB[ii*(*lab)+jj]);
      }
      fprintf(file,"\n");
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_GB2AIJ_operator_poisson1D(double* AB, int* la, char* filename){
  FILE * file;
  int jj;
  file = fopen(filename, "w");
  //Numbering from 1 to la
  if (file != NULL){
    for (jj=1;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj,jj+1,AB[(*la)+jj]);
    }
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+1,jj+1,AB[2*(*la)+jj]);
    }
    for (jj=0;jj<(*la)-1;jj++){
      fprintf(file,"%d\t%d\t%lf\n",jj+2,jj+1,AB[3*(*la)+jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  }
}

void write_vec(double* vec, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\n",vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  

void write_xy(double* vec, double* x, int* la, char* filename){
  int jj;
  FILE * file;
  file = fopen(filename, "w");
  // Numbering from 1 to la
  if (file != NULL){
    for (jj=0;jj<(*la);jj++){
      fprintf(file,"%lf\t%lf\n",x[jj],vec[jj]);
    }
    fclose(file);
  }
  else{
    perror(filename);
  } 
}  


void eig_poisson1D(double* eigval, int *la){
  int i;
  double scal;
  for (i=0; i< *la; i++){
    scal = (1.0*i+1.0)* M_PI_2 *(1.0/(*la+1));
    eigval[i] = sin(scal);
    eigval[i] = 4 * eigval[i] * eigval[i];
  } 
}

double eigmax_poisson1D(int *la){
  double eigmax;
  eigmax = sin(*la * M_PI_2 *(1.0/(*la+1)));
  eigmax = 4 * eigmax * eigmax;
  return eigmax;
}

double eigmin_poisson1D(int *la){
  double eigmin;
  eigmin = sin(M_PI_2 * (1.0/(*la+1)));
  eigmin = 4 * eigmin * eigmin;
  return eigmin;
}

double richardson_alpha_opt(int *la){
    return 2.0 / (eigmax_poisson1D(la) + eigmin_poisson1D(la));
}

void richardson_alpha(double *AB, double *RHS, double *X, double *alpha_rich, int *lab, 
int *la,int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){

  double res_norm;

    cblas_dcopy(*la, RHS, 1, resvec, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X,1, 1.0, resvec, 1);
    cblas_daxpy(*la, *alpha_rich, resvec, 1, X, 1);
    res_norm = cblas_dnrm2(*la, resvec, 1);
    res_norm /= cblas_dnrm2(*la, RHS, 1);

  // x(k+1) = x(k) + alpha(b - Ax(k))
  while (res_norm > *tol && *nbite < *maxit) {

    cblas_dcopy(*la, RHS, 1, resvec, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X,1, 1.0, resvec, 1);
    cblas_daxpy(*la, *alpha_rich, resvec, 1, X, 1);
    res_norm = cblas_dnrm2(*la, resvec, 1);
    res_norm /= cblas_dnrm2(*la, RHS, 1);
    //printf("%d %1.6f\n", *nbit, res_norm);
    

    (*nbite)++;
  }

  //printf("res_norm Richardson_ALPHA =  %1.6f\n",res_norm);
  printf("Nombre iteration avec Richardson_ALPHA = %d \n", *nbite);

}


void extract_MB_jacobi_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){
  int i,j;

  for (i = 0; i < (*la); i++)
  {
    //recuperer l'index 
    j = indexABCol(0,i,lab);
    //recuperer l'inverse du diagonale de AB dans MB
    MB[j + *kv] = 1 / AB[j + *kv];
  }
}

void extract_MB_gauss_seidel_tridiag(double *AB, double *MB, int *lab, int *la,int *ku, int*kl, int *kv){

  int i;
  //la taille du matrice MB_Gauss_seidel
  int size = (*la) * (*la);
  //matrice temporaire
  double *temp = (double *) calloc(size, sizeof(double)); 

  for (i = 0; i < *la; i++) 
  { 
    //recuperer les indices
    int k = indexABCol(i,i,la); 
    int kk = indexABCol(1,i,lab);
    temp[k] = AB[kk]; 
  }

  for (i = 0; i < *la - 1; i++) 
  {
    int k = indexABCol(i,i+1,la);
    int kk = indexABCol(2,i,lab);
    temp[k] = AB[kk];
  }

  int *ipiv = malloc(sizeof(int) * (*la));
  //factorisation LU
  LAPACKE_dgetrf(CblasRowMajor, *la, *la, temp, *la, ipiv);
  //inverse de temp
  LAPACKE_dgetri(CblasRowMajor, *la, temp, *la, ipiv);
  free(ipiv);

  //recopier temp dans MB
  cblas_dcopy(size, temp, 1, MB, 1);
}

//avec extraction jacobi
void richardson_MB_Jacobi(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,
int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
  double res_norm;

    cblas_dcopy(*la, RHS, 1, resvec, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X,1, 1.0, resvec, 1);
    res_norm = cblas_dnrm2(*la, resvec, 1);
    res_norm /= cblas_dnrm2(*la, RHS, 1);

  // x(k+1) = x(k) + alpha(b - Ax(k))
  while (res_norm > *tol && *nbite < *maxit) {

    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, 1.0, MB, *lab, resvec,1, 1.0, X, 1);

    cblas_dcopy(*la, RHS, 1, resvec, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X,1, 1.0, resvec, 1);
    res_norm = cblas_dnrm2(*la, resvec, 1);
    res_norm /= cblas_dnrm2(*la, RHS, 1);
    //printf("%d %1.6f\n", *nbit, res_norm);
    

    (*nbite)++;
  }

  //printf("res_norm Richardson_MB_JACOBI =  %1.6f\n",res_norm);
  printf("Nombre iteration avec Richardson_MB_JACOBI = %d \n", *nbite);

}

//avec extraction gauss
void richardson_MB_Gauss(double *AB, double *RHS, double *X, double *MB, int *lab, int *la,
int *ku, int*kl, double *tol, int *maxit, double *resvec, int *nbite){
    double res_norm;

    cblas_dcopy(*la, RHS, 1, resvec, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X,1, 1.0, resvec, 1);
    res_norm = cblas_dnrm2(*la, resvec, 1);
    res_norm /= cblas_dnrm2(*la, RHS, 1);

  // x(k+1) = x(k) + alpha(b - Ax(k))
  while (res_norm > *tol && *nbite < *maxit) {

    cblas_dgemv(CblasColMajor, CblasNoTrans, *la, *la, 1.0, MB, *la, resvec,1, 1.0, X, 1);

    cblas_dcopy(*la, RHS, 1, resvec, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, *la, *la, *kl, *ku, -1.0, AB, *lab, X,1, 1.0, resvec, 1);
    res_norm = cblas_dnrm2(*la, resvec, 1);
    res_norm /= cblas_dnrm2(*la, RHS, 1);
    //printf("%d %1.6f\n", *nbit, res_norm);
    

    (*nbite)++;
  }

  //printf("res_norm =  %1.6f\n",res_norm);
  printf("Nombre iteration avec Richardson_MB_GAUSS = %d \n", *nbite);

}

int indexABCol(int i, int j, int *lab){
  return j*(*lab)+i;
}

int dgbtrftridiag(int *la, int*n, int *kl, int *ku, double *AB, int *lab, int *ipiv, int *info, int* kv){
  int i, j, k, kk = 3; //kk : nbr d colonnes

      if (*kv>=0){
        kk = 4; 
        for (i=0;i< *kv;i++){
            AB[i]=0.0;
        }
      }

      AB[*kv+2]/=AB[*kv+1]; //modification

    for (j=1;j<(*la);j++)
    {
      if (*kv>=0){
        for (i=0;i< *kv;i++)
        {
            k = indexABCol(i,j,lab); //recupere l'index
            AB[k]=0.0;
        }
      }

      k = indexABCol(0,j,lab);
      AB[k+ *kv+1] -= AB[k+ *kv] * AB[(k-kk)+ *kv+2];
      AB[k+ *kv+2] /= AB[k+ *kv+1];
    }

  return *info;
}

