#ifndef TC_VEC3D_H
#define TC_VEC3D_H

#include <math.h> // cmath on g++ ?

void set_vec3d(int * v1,const int * v2);
void set_vec3d(double * v1,const double * v2);
void add_vec3d(double * v1,const double * v2);
void subtract_vec3d(double * v1,const double * v2);
void normalize_vec3d(double * v1);
double dotproduct_vec3d(const double * v1,const double * v2);
void divide_vec3d(double * v1,const double scalar);

inline double vecDotVec3d(const double * v1, const double * v2) {

  double dotproduct = v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];

  return(dotproduct);
}

inline double normVec3d(double * v1) {

  double mag = vecDotVec3d(v1, v1);

  mag = sqrt(mag);

  v1[0] /= mag;
  v1[1] /= mag;
  v1[2] /= mag;

  return mag;
}

inline double normVec3d(double *vNorm, const double *v1) {

  double mag = vecDotVec3d(v1, v1);

  mag = sqrt(mag);

  vNorm[0] = v1[0]/mag;
  vNorm[1] = v1[1]/mag;
  vNorm[2] = v1[2]/mag;

  return mag;
}

inline void vecMinVec3d(double *vRes, const double *v1, const double *v2) {

  vRes[0] = v1[0] - v2[0];
  vRes[1] = v1[1] - v2[1];
  vRes[2] = v1[2] - v2[2];
}

inline void matMultMatR3(double (*res)[3], double (*m1)[3], double (*m2)[3])  // m*m matrices
{
  for (int i=0; i<3; i++)
  for (int j=0; j<3; j++)
  {
    res[i][j] = 0.0;
    for (int k=0; k<3; k++)
      res[i][j] += m1[i][k]*m2[k][j];
  }
}

inline void transMatR3(double (*m)[3])
{
  double temp;

  temp    = m[0][1];
  m[0][1] = m[1][0];
  m[1][0] = temp;

  temp    = m[0][2];
  m[0][2] = m[2][0];
  m[2][0] = temp;

  temp    = m[1][2];
  m[1][2] = m[2][1];
  m[2][1] = temp;
}

inline void similMatR3(double (*res)[3], double (*m1)[3], double (*m2)[3])
{
  for (int k=0; k<3; k++)
    for (int l=0; l<3; l++)
    {
      res[k][l] = 0.0;
      for (int i=0; i<3; i++)
        for (int j=0; j<3; j++)
          res[k][l] += m1[i][k]*m1[j][l]*m2[i][j];
    }
}

inline void matMultMatR5(double (*res)[5], double (*m1)[5], double (*m2)[5])  // m*m matrices
{
  for (int i=0; i<5; i++)
  for (int j=0; j<5; j++)
  {
    res[i][j] = 0.0;
    for (int k=0; k<5; k++)
      res[i][j] += m1[i][k]*m2[k][j];
  }
}

inline void matMultDiagMatR5(double (*res)[5], double (*mat)[5], double *diag)  // m*m matrices
{
  for (int i=0; i<5; i++)
  for (int j=0; j<5; j++)
    res[i][j] = mat[i][j]*diag[j];
}

inline void matMultVecR5(double *res, double (*mat)[5], double *vec)  // m*m matrices !!!!!!
{
  for (int i=0; i<5; i++)
  {
    res[i] = 0.0;

    for (int j=0; j<5; j++)
      res[i] += mat[i][j]*vec[j];
  }
}







#endif




