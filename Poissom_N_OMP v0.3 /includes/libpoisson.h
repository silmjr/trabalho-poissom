#ifndef LIBPOISSON_H_INCLUDED
#define LIBPOISSON_H_INCLUDED

typedef struct node
{
    double up;
    double rh;
    double dn;
    double lf;
} nodeSides;

typedef struct material
{
    double f;       /* fonte */
    double perm;    /* permeabilidade das rochas */
    double shi;     /* variável auxiliar valor de shi */
} nodeMaterial;

double walltime( double* );
void errorMsg(char error_text[]);
nodeSides* criaVetorNode(int lin,int col);
nodeSides** criaMatrizNode(int lin,int col, nodeSides* v_aux);
nodeSides** montaMatrizNode(int lin,int col);
double* criaVetor(int lin,int col);
double** criaMatriz(int lin,int col, double* v_aux);
double** montaMatriz(int lin,int col);
nodeMaterial* criaVetorMaterial(int lin,int col);
nodeMaterial** criaMatrizMaterial(int lin,int col, nodeMaterial* v_aux);
nodeMaterial** montaMatrizMaterial(int lin,int col);

void canto_d_l(const int i, const int j,    nodeMaterial **pMat,
               nodeSides **beta,
               nodeSides **q,
               nodeSides **q_old,
               nodeSides **l_old,
               double **p);

void canto_d_lArray(const int i, const int j, const int N,
               nodeMaterial *pMat,
               nodeSides *beta,
               nodeSides *q,
               nodeSides *q_old,
               nodeSides *l_old,
               double *p, double *Media);

void canto_u_l(const int i, const int j,    nodeMaterial **pMat,
               nodeSides **beta,
               nodeSides **q,
               nodeSides **q_old,
               nodeSides **l_old,
               double **p);

void canto_u_lArray(const int i, const int j, const int N,
                    nodeMaterial *pMat,
                    nodeSides *beta,
                    nodeSides *q,
                    nodeSides *q_old,
                    nodeSides *l_old,
                    double *p, double *Media);

void canto_d_r(const int i, const int j,    nodeMaterial **pMat,
               nodeSides **beta,
               nodeSides **q,
               nodeSides **q_old,
               nodeSides **l_old,
               double **p);

void canto_d_rArray(const int i, const int j, const int N,
                    nodeMaterial *pMat,
                    nodeSides *beta,
                    nodeSides *q,
                    nodeSides *q_old,
                    nodeSides *l_old,
                    double *p, double *Media);

void canto_u_r(const int i, const int j,    nodeMaterial **pMat,
               nodeSides **beta,
               nodeSides **q,
               nodeSides **q_old,
               nodeSides **l_old,
               double **p);

void canto_u_rArray(const int i, const int j, const int N,
                    nodeMaterial *pMat,
                    nodeSides *beta,
                    nodeSides *q,
                    nodeSides *q_old,
                    nodeSides *l_old,
                    double *p, double *Media);

void fronteira_u(const int n, const int j,
                 nodeMaterial **pMat,
                 nodeSides **beta,
                 nodeSides **q,
                 nodeSides **q_old,
                 nodeSides **l_old,
                 double **p);

void fronteira_uArray(const int N, const int j,
                      nodeMaterial *pMat,
                      nodeSides *beta,
                      nodeSides *q,
                      nodeSides *q_old,
                      nodeSides *l_old,
                      double *p, double *Media);

void fronteira_d(const int n, const int j,
                 nodeMaterial **pMat,
                 nodeSides **beta,
                 nodeSides **q,
                 nodeSides **q_old,
                 nodeSides **l_old,
                 double **p);

void fronteira_dArray(const int N, const int j,
                      nodeMaterial *pMat,
                      nodeSides *beta,
                      nodeSides *q,
                      nodeSides *q_old,
                      nodeSides *l_old,
                      double *p, double *Media);

void fronteira_r(const int i, const int n,
                 nodeMaterial **pMat,
                 nodeSides **beta,
                 nodeSides **q,
                 nodeSides **q_old,
                 nodeSides **l_old,
                 double **p);

void fronteira_rArray(const int i, const int N,
                      nodeMaterial *pMat,
                      nodeSides *beta,
                      nodeSides *q,
                      nodeSides *q_old,
                      nodeSides *l_old,
                      double *p, double *Media);

void fronteira_l(const int i, const int n,
                 nodeMaterial **pMat,
                 nodeSides **beta,
                 nodeSides **q,
                 nodeSides **q_old,
                 nodeSides **l_old,
                 double **p);

void fronteira_lArray(const int i, const int N,
                      nodeMaterial *pMat,
                      nodeSides *beta,
                      nodeSides *q,
                      nodeSides *q_old,
                      nodeSides *l_old,
                      double *p, double *Media);

void internos(const int n,
              nodeMaterial **pMat,
              nodeSides **beta,
              nodeSides **q,
              nodeSides **q_old,
              nodeSides **l_old,
              double **p);

double internosArray(const int N,
              nodeMaterial *pMat,
              nodeSides *beta,
              nodeSides *q,
              nodeSides *q_old,
              nodeSides *l_old,
              double *p, int threadID, int localSize, int amountThreads);

double lagrangeUpdate(const int n,
                      nodeSides **beta,
                      nodeSides **q,
                      nodeSides **q_old,
                      nodeSides **l,
                      nodeSides **l_old,
                      double **p);

double lagrangeUpdateArray(const int N,
                      nodeSides *beta,
                      nodeSides *q,
                      nodeSides *q_old,
                      nodeSides *l,
                      nodeSides *l_old,
                      double *p, double *p_old, double *Media);

double mediaZero(const int n,
                 double Media,
                 nodeSides **l,
                 double **p,
                 double **p_old);

double mediaZeroArray(const int N,
                 double Media,
                 nodeSides *l,
                 double *p,
                 double *p_old);


#endif // LIBPOISSON_H_INCLUDED

