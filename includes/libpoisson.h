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
double* criaVetor(int lin,int col);
double** criaMatriz(int lin,int col, double* v_aux);
nodeMaterial* criaVetorMaterial(int lin,int col);
nodeMaterial** criaMatrizMaterial(int lin,int col, nodeMaterial* v_aux);

void canto_d_l(const int i, const int j,    nodeMaterial **pMat,
               nodeSides **beta,
               nodeSides **q,
               nodeSides **q_old,
               nodeSides **l_old,
               double **p);

void canto_u_l(const int i, const int j,    nodeMaterial **pMat,
               nodeSides **beta,
               nodeSides **q,
               nodeSides **q_old,
               nodeSides **l_old,
               double **p);

void canto_d_r(const int i, const int j,    nodeMaterial **pMat,
               nodeSides **beta,
               nodeSides **q,
               nodeSides **q_old,
               nodeSides **l_old,
               double **p);

void canto_u_r(const int i, const int j,    nodeMaterial **pMat,
               nodeSides **beta,
               nodeSides **q,
               nodeSides **q_old,
               nodeSides **l_old,
               double **p);

void fronteira_u(const int i, const int j,    nodeMaterial **pMat,
                 nodeSides **beta,
                 nodeSides **q,
                 nodeSides **q_old,
                 nodeSides **l_old,
                 double **p);

void fronteira_d(const int i, const int j,    nodeMaterial **pMat,
                 nodeSides **beta,
                 nodeSides **q,
                 nodeSides **q_old,
                 nodeSides **l_old,
                 double **p);

void fronteira_r(const int i, const int j,    nodeMaterial **pMat,
                 nodeSides **beta,
                 nodeSides **q,
                 nodeSides **q_old,
                 nodeSides **l_old,
                 double **p);

void fronteira_l(const int i, const int j,    nodeMaterial **pMat,
                 nodeSides **beta,
                 nodeSides **q,
                 nodeSides **q_old,
                 nodeSides **l_old,
                 double **p);

void internos(const int i, const int j,    nodeMaterial **pMat,
            nodeSides **beta,
            nodeSides **q,
            nodeSides **q_old,
            nodeSides **l_old,
            double **p);


#endif // LIBPOISSON_H_INCLUDED

