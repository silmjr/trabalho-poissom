#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "../includes/libpoisson.h"

/* Função que retorna o tempo em função do relógio */
double walltime( double *t0 )
{

    double mic, time;
    double mega = 0.000001;
    struct timeval tp;
    struct timezone tzp;
    static long base_sec = 0;
    static long base_usec = 0;

    (void) gettimeofday(&tp,&tzp);

    if (base_sec == 0)
    {
        base_sec = tp.tv_sec;
        base_usec = tp.tv_usec;
    }

    time = (double) (tp.tv_sec - base_sec);
    mic = (double) (tp.tv_usec - base_usec);
    time = (time + mic * mega) - *t0;
    return(time);
}

void errorMsg(char error_text[])
/* Standard error handler */
{
    fprintf(stderr,"Run-time error...\n");
    fprintf(stderr,"%s\n",error_text);
    fprintf(stderr,"...now exiting to system...\n");
    exit(1);
}

/*Função que aloca as matrizes de nós, como arrays, de tamanho linXcol*/
nodeSides* criaVetorNode(int lin,int col)
{
    nodeSides* v_aux;

    v_aux = (nodeSides*) malloc(lin*col*sizeof(nodeSides));

    if(v_aux==NULL)
        errorMsg("allocation failure in dvector()");

    return v_aux;
}

/*Função que aloca as matrizes de nós de tamanho linXcol*/
nodeSides** criaMatrizNode(int lin,int col, nodeSides* v_aux)
{
    int i;
    nodeSides** p_aux;

    p_aux = (nodeSides**) malloc(lin*sizeof(nodeSides*));

    if(p_aux==NULL)
        errorMsg("allocation failure in dvector()");

    for(i=0; i<lin; i++)
        p_aux[i] = (nodeSides*) &v_aux[i*col];

    return p_aux;
}

/*Função que aloca as matrizes, como arrays, de tamanho linXcol*/
double* criaVetor(int lin,int col)
{
    double* v_aux;

    v_aux = (double*) malloc(lin*col*sizeof(double));

    if(v_aux==NULL)
        errorMsg("allocation failure in dvector()");

    return v_aux;

}

/*Função que aloca as matrizes de tamanho linXcol*/
double** criaMatriz(int lin,int col, double* v_aux)
{
    int i;
    double** p_aux;

    p_aux = (double**) malloc(lin*sizeof(double*));

    if(p_aux==NULL)
        errorMsg("allocation failure in dvector()");

    for(i=0; i<lin; i++)
        p_aux[i] = (double*) &v_aux[i*col];

    return p_aux;
}

/*Função que aloca a matriz de parametros materiais, como
 *arrays, de tamanho linXcol
 */
nodeMaterial* criaVetorMaterial(int lin,int col)
{
    nodeMaterial* v_aux;

    v_aux = (nodeMaterial*) malloc(lin*col*sizeof(nodeMaterial));

    if(v_aux==NULL)
        errorMsg("allocation failure in dvector()");

    return v_aux;
}

/*Função que aloca a matriz de parametros materiais
 *de tamanho linXcol
 */
nodeMaterial** criaMatrizMaterial(int lin,int col, nodeMaterial* v_aux)
{
    int i;
    nodeMaterial** p_aux;

    p_aux = (nodeMaterial**) malloc(lin*sizeof(nodeMaterial*));

    if(p_aux==NULL)
        errorMsg("allocation failure in dvector()");

    for(i=0; i<lin; i++)
        p_aux[i] = (nodeMaterial*) &v_aux[i*col];

    return p_aux;
}

/* Função para o canto inferior esquerdo*/
void canto_d_l(const int i, const int j,
               nodeMaterial **pMat,
               nodeSides **beta,
               nodeSides **q,
               nodeSides **q_old,
               nodeSides **l_old,
               double **p)

{
    double AuxU,AuxR,DU,DR;	/* Auxiliares para cada lado das células */

    AuxU = pMat[i][j].shi/(1+beta[i][j].up*pMat[i][j].shi);
    AuxR = pMat[i][j].shi/(1+beta[i][j].rh*pMat[i][j].shi);
    DU = AuxU*(beta[i][j].up*q_old[i][j+1].dn+l_old[i][j+1].dn);
    DR = AuxR*(beta[i][j].rh*q_old[i+1][j].lf+l_old[i+1][j].lf);
    p[i][j] = (pMat[i][j].f + DU + DR)/(AuxU + AuxR);
    //printf("Cando_d_l p = %f\n",p[i][j]);
    q[i][j].up = AuxU*p[i][j] - DU;
    q[i][j].rh = AuxR*p[i][j] - DR;
}

/* Função para o canto superior esquerdo*/
void canto_u_l(const int i, const int j,    nodeMaterial **pMat,
               nodeSides **beta,
               nodeSides **q,
               nodeSides **q_old,
               nodeSides **l_old,
               double **p)
{
    double AuxD, AuxR,DD, DR;	/* Auxiliares para cada lado das células */

    AuxD = pMat[i][j].shi/(1+beta[i][j].dn*pMat[i][j].shi);
    AuxR = pMat[i][j].shi/(1+beta[i][j].rh*pMat[i][j].shi);
    DD = AuxD*(beta[i][j].dn*q_old[i][j-1].up+l_old[i][j-1].up);
    DR = AuxR*(beta[i][j].rh*q_old[i+1][j].lf+l_old[i+1][j].lf);
    p[i][j] = (pMat[i][j].f + DD + DR)/(AuxD + AuxR);
    q[i][j].dn = AuxD*p[i][j] - DD;
    q[i][j].rh = AuxR*p[i][j] - DR;
}

/* Função para o canto inferior direito*/
void canto_d_r(const int i, const int j,    nodeMaterial **pMat,
               nodeSides **beta,
               nodeSides **q,
               nodeSides **q_old,
               nodeSides **l_old,
               double **p)
{

    double AuxU, AuxL,DU,DL;	/* Auxiliares para cada lado das células */

    AuxU = pMat[i][j].shi/(1+beta[i][j].up*pMat[i][j].shi);
    AuxL = pMat[i][j].shi/(1+beta[i][j].lf*pMat[i][j].shi);
    DU = AuxU*(beta[i][j].up*q_old[i][j+1].dn+l_old[i][j+1].dn);
    DL = AuxL*(beta[i][j].lf*q_old[i-1][j].rh+l_old[i-1][j].rh);
    p[i][j] = (pMat[i][j].f + DU + DL)/(AuxU + AuxL);
    q[i][j].up = AuxU*p[i][j] - DU;
    q[i][j].lf = AuxL*p[i][j] - DL;
}

/* Funcao para o canto superior dereito*/
void canto_u_r(const int i, const int j,    nodeMaterial **pMat,
               nodeSides **beta,
               nodeSides **q,
               nodeSides **q_old,
               nodeSides **l_old,
               double **p)
{

    double AuxD,AuxL,DD,DL;	/* Auxiliares para cada lado das células */

    AuxD = pMat[i][j].shi/(1+beta[i][j].dn*pMat[i][j].shi);
    AuxL = pMat[i][j].shi/(1+beta[i][j].lf*pMat[i][j].shi);
    DD = AuxD*(beta[i][j].dn*q_old[i][j-1].up+l_old[i][j-1].up);
    DL = AuxL*(beta[i][j].lf*q_old[i-1][j].rh+l_old[i-1][j].rh);
    p[i][j] = (pMat[i][j].f + DD + DL)/(AuxD + AuxL);
    q[i][j].dn = AuxD*p[i][j] - DD;
    q[i][j].lf = AuxL*p[i][j] - DL;
}  //OKKKKKKKKKK

/* Função para a fronteira superior U */
void fronteira_u(const int i, const int j,    nodeMaterial **pMat,
                 nodeSides **beta,
                 nodeSides **q,
                 nodeSides **q_old,
                 nodeSides **l_old,
                 double **p)
{

    double AuxD, AuxR, AuxL, DD, DR, DL;	/* Auxiliares para cada lado das celulas */

    AuxL = pMat[i][j].shi/(1+beta[i][j].lf*pMat[i][j].shi);
    AuxR = pMat[i][j].shi/(1+beta[i][j].rh*pMat[i][j].shi);
    AuxD = pMat[i][j].shi/(1+beta[i][j].dn*pMat[i][j].shi);
    DL = AuxL*(beta[i][j].lf*q_old[i-1][j].rh+l_old[i-1][j].rh);
    DR = AuxR*(beta[i][j].rh*q_old[i+1][j].lf+l_old[i+1][j].lf);
    DD = AuxD*(beta[i][j].dn*q_old[i][j-1].up+l_old[i][j-1].up);
    p[i][j] = (pMat[i][j].f + DD + DL + DR)/(AuxD + AuxL+AuxR);
    q[i][j].lf = AuxL*p[i][j] - DL;
    q[i][j].rh = AuxR*p[i][j] - DR;
    q[i][j].dn = AuxD*p[i][j] - DD;
}

/* Função para a fronteira inferior D */
void fronteira_d(const int i, const int j,    nodeMaterial **pMat,
                 nodeSides **beta,
                 nodeSides **q,
                 nodeSides **q_old,
                 nodeSides **l_old,
                 double **p)
{
    double AuxU,AuxR, AuxL,DU,DR, DL;	/* Auxiliares para cada lado das celulas */

    AuxL = pMat[i][j].shi/(1+beta[i][j].lf*pMat[i][j].shi);
    AuxR = pMat[i][j].shi/(1+beta[i][j].rh*pMat[i][j].shi);
    AuxU = pMat[i][j].shi/(1+beta[i][j].up*pMat[i][j].shi);
    DL = AuxL*(beta[i][j].lf*q_old[i-1][j].rh+l_old[i-1][j].rh);
    DR = AuxR*(beta[i][j].rh*q_old[i+1][j].lf+l_old[i+1][j].lf);
    DU = AuxU*(beta[i][j].up*q_old[i][j+1].dn+l_old[i][j+1].dn);
    p[i][j] = (pMat[i][j].f + DU + DL + DR)/(AuxU + AuxL+AuxR);
    q[i][j].lf = AuxL*p[i][j] - DL;
    q[i][j].rh = AuxR*p[i][j] - DR;
    q[i][j].up = AuxU*p[i][j] - DU;

}

/* Função para a fronteira dereita R */
void fronteira_r(const int i, const int j,    nodeMaterial **pMat,
                 nodeSides **beta,
                 nodeSides **q,
                 nodeSides **q_old,
                 nodeSides **l_old,
                 double **p)
{

    double AuxU,AuxD, AuxL,DU, DD, DL;	/* Auxiliares para cada lado das celulas */

    AuxU = pMat[i][j].shi/(1+beta[i][j].up*pMat[i][j].shi);
    AuxD = pMat[i][j].shi/(1+beta[i][j].dn*pMat[i][j].shi);
    AuxL = pMat[i][j].shi/(1+beta[i][j].lf*pMat[i][j].shi);
    DU = AuxU*(beta[i][j].up*q_old[i][j+1].dn+l_old[i][j+1].dn);
    DD = AuxD*(beta[i][j].dn*q_old[i][j-1].up+l_old[i][j-1].up);
    DL = AuxL*(beta[i][j].lf*q_old[i-1][j].rh+l_old[i-1][j].rh);
    p[i][j] = (pMat[i][j].f + DU + DL + DD)/(AuxU + AuxL+AuxD);
    q[i][j].up = AuxU*p[i][j] - DU;
    q[i][j].dn = AuxD*p[i][j] - DD;
    q[i][j].lf = AuxL*p[i][j] - DL;

}

/* Função para a fronteira esquerda L */
void fronteira_l(const int i, const int j,    nodeMaterial **pMat,
                 nodeSides **beta,
                 nodeSides **q,
                 nodeSides **q_old,
                 nodeSides **l_old,
                 double **p)
{

    double AuxU, AuxD, AuxR,DU, DD, DR; /* Auxiliares para cada lado das células */

    AuxU = pMat[i][j].shi/(1+beta[i][j].up*pMat[i][j].shi);
    AuxD = pMat[i][j].shi/(1+beta[i][j].dn*pMat[i][j].shi);
    AuxR = pMat[i][j].shi/(1+beta[i][j].rh*pMat[i][j].shi);
    DU = AuxU*(beta[i][j].up*q_old[i][j+1].dn+l_old[i][j+1].dn);
    DD = AuxD*(beta[i][j].dn*q_old[i][j-1].up+l_old[i][j-1].up);
    DR = AuxR*(beta[i][j].rh*q_old[i+1][j].lf+l_old[i+1][j].lf);
    p[i][j] = (pMat[i][j].f + DU + DR + DD)/(AuxU + AuxR+AuxD);
    q[i][j].up = AuxU*p[i][j] - DU;
    q[i][j].dn = AuxD*p[i][j] - DD;
    q[i][j].rh = AuxR*p[i][j] - DR;
}

/* Função para a fronteira esquerda L */
void internos(const int i, const int j,    nodeMaterial **pMat,
            nodeSides **beta,
            nodeSides **q,
            nodeSides **q_old,
            nodeSides **l_old,
            double **p)
{
    double AuxU, AuxD, AuxR, AuxL, DU, DD, DR, DL; /* Auxiliares para cada lado das células */

    AuxU = pMat[i][j].shi/(1+beta[i][j].up*pMat[i][j].shi);
    AuxR = pMat[i][j].shi/(1+beta[i][j].rh*pMat[i][j].shi);
    AuxD = pMat[i][j].shi/(1+beta[i][j].dn*pMat[i][j].shi);
    AuxL = pMat[i][j].shi/(1+beta[i][j].lf*pMat[i][j].shi);

    DL = AuxL*(beta[i][j].lf*q_old[i-1][j].rh+l_old[i-1][j].rh);
    DD = AuxD*(beta[i][j].dn*q_old[i][j-1].up+l_old[i][j-1].up);
    DR = AuxR*(beta[i][j].rh*q_old[i+1][j].lf+l_old[i+1][j].lf);
    DU = AuxU*(beta[i][j].up*q_old[i][j+1].dn+l_old[i][j+1].dn);

    p[i][j] = (pMat[i][j].f + DU + DD + DL + DR)/(AuxU + AuxD + AuxL+AuxR);

    q[i][j].up = AuxU*p[i][j] - DU;
    q[i][j].rh = AuxR*p[i][j] - DR;
    q[i][j].dn = AuxD*p[i][j] - DD;
    q[i][j].lf = AuxL*p[i][j] - DL;
}

