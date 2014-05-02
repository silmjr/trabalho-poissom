#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <math.h>

#include "includes/libpoisson.h"

#define N 258

//int main(int argc, char *argv[])
int main(void)
{
    /*
    if(argc != 2)
    {
        errorMsg("Erro na execucao!\nEx. ./[nome_executavel] [tamanho_da_malha+2]\n");
        return (1);
    }
    else
    {
        N = atoi(argv[1]);
    }
    */

    /*
     * Fluxos atuais e antigos em cada um dos lados da célula espacial
     * up (Uper), dn (Down), lf (Left), rh (Right)
     */
    nodeSides *q, **mq, *q_old, **mq_old;

    q = criaVetorNode(N, N);
    mq = criaMatrizNode(N, N, q);
    q_old = criaVetorNode(N, N);
    mq_old = criaMatrizNode(N, N, q_old);

    /*
     *  Multiplicadores de Lagrange em cada um dos lados da célula espacial
     */
    nodeSides *l, **ml, *l_old, **ml_old;

    l = criaVetorNode(N, N);
    ml = criaMatrizNode(N, N, l);
    l_old = criaVetorNode(N, N);
    ml_old = criaMatrizNode(N, N, l_old);

    /*
     *	Betas da condição de Robin em cada um dos lados da célula espacial
     */

    nodeSides *beta, **mbeta;

    beta = criaVetorNode(N, N);
    mbeta = criaMatrizNode(N, N, beta);

    /*
     *	Pressões atuais e antigas em cada uma das células
     */

    double *p, *p_old, **mp, **mp_old;

    p = criaVetor(N, N);
    mp = criaMatriz(N, N, p);
    p_old = criaVetor(N, N);
    mp_old = criaMatriz(N, N, p_old);

    /*
     *  Parametros materiais da grade
     */

    nodeMaterial *pMat, **mpMat;
    pMat = criaVetorMaterial(N, N);
    mpMat = criaMatrizMaterial(N, N, pMat);

    double lsize = 25600.00; /* dimensão da regiao */

    /* Variáveis para auxiliar na contagem do tempo*/
    double start, stop;
    double startTime, elapsedTime;
    double clockZero = 0.0;
    clock_t ticks1, ticks2;
    /* Variáveis auxiliares na implementação */
    int i, j, k, n;
    double  h, aux,
            c = 1.,                 /* Valor para o calculo de beta */
            Media, erro,            /* Media das pressoes e erro na norma */
            sum1, sum2,             /* Auxiliares somadores */
            Keff;                   /* Valor medio da permeabilidade */
    double **mp_aux;                /* Ponteiro para troca */
    nodeSides **mq_aux;

    start = omp_get_wtime();
    startTime = walltime( &clockZero );
    ticks1 = clock();

    /*
    * Inicialização das variaveis
    */

    for(i=1; i < N-1; i++)
        for(j=1; j < N-1; j++)
        {
            mpMat[i][j].f = 0.0;
            ml[i][j].dn = ml[i][j].lf = ml[i][j].rh = ml[i][j].up = 0.0;
            mq[i][j].dn = mq[i][j].lf = mq[i][j].rh = mq[i][j].up = 0.0;
            mp[i][j]=0.0;
        }

    /*
     * Cálculo de algumas variáveis auxiliares para o valor de permeabilidade
     * en blocos verticais
     */
    n = N-2;
    h = lsize/n;
    k = n/2;

    for(i=1; i<=k; i++)
        for(j=1; j<N-1; j++)
        {
            mpMat[i][j].perm = 1.0e-10;
            mpMat[i+k][j].perm = 1.0e-11;
        }

    /*
     * Calcula os Beta da Condicao de Robin
     */

    Keff = (2*mpMat[1][1].perm*mpMat[1][2].perm)/(mpMat[1][1].perm + mpMat[1][2].perm);
    mbeta[1][1].up = c*h/Keff;
    Keff = (2*mpMat[1][1].perm*mpMat[2][1].perm)/(mpMat[1][1].perm + mpMat[2][1].perm);
    mbeta[1][1].rh = c*h/Keff;

    Keff = (2*mpMat[1][n].perm*mpMat[1][n-1].perm)/(mpMat[1][n].perm + mpMat[1][n-1].perm);
    mbeta[1][n].dn = c*h/Keff;
    Keff = (2*mpMat[1][n].perm*mpMat[2][n].perm)/(mpMat[1][n].perm + mpMat[2][n].perm);
    mbeta[1][n].rh = c*h/Keff;

    Keff = (2*mpMat[n][1].perm*mpMat[n][2].perm)/(mpMat[n][1].perm + mpMat[n][2].perm);
    mbeta[n][1].up = c*h/Keff;
    Keff = (2*mpMat[n][1].perm*mpMat[n-1][1].perm)/(mpMat[n][1].perm + mpMat[n-1][1].perm);
    mbeta[n][1].lf = c*h/Keff;

    Keff = (2*mpMat[n][n].perm*mpMat[n][n-1].perm)/(mpMat[n][n].perm + mpMat[n][n-1].perm);
    mbeta[n][n].dn = c*h/Keff;
    Keff = (2*mpMat[n][n].perm*mpMat[n-1][n].perm)/(mpMat[n][n].perm + mpMat[n-1][n].perm);
    mbeta[n][n].lf = c*h/Keff;

    for (i=2; i<n; i++)
    {
        Keff = (2*mpMat[i][1].perm*mpMat[i][2].perm)/(mpMat[i][1].perm + mpMat[i][2].perm);
        mbeta[i][1].up= c*h/Keff;
        Keff = (2*mpMat[i][1].perm*mpMat[i-1][1].perm)/(mpMat[i][1].perm + mpMat[i-1][1].perm);
        mbeta[i][1].lf= c*h/Keff;
        Keff = (2*mpMat[i][1].perm*mpMat[i+1][1].perm)/(mpMat[i][1].perm + mpMat[i+1][1].perm);
        mbeta[i][1].rh= c*h/Keff;

        Keff = (2*mpMat[i][n].perm*mpMat[i][n-1].perm)/(mpMat[i][n].perm + mpMat[i][n-1].perm);
        mbeta[i][n].dn= c*h/Keff;
        Keff = (2*mpMat[i][n].perm*mpMat[i-1][n].perm)/(mpMat[i][n].perm + mpMat[i-1][n].perm);
        mbeta[i][n].lf= c*h/Keff;
        Keff = (2*mpMat[i][n].perm*mpMat[i+1][n].perm)/(mpMat[i][n].perm + mpMat[i+1][n].perm);
        mbeta[i][n].rh= c*h/Keff;


        Keff = (2*mpMat[1][i].perm*mpMat[1][i+1].perm)/(mpMat[1][i].perm + mpMat[1][i+1].perm);
        mbeta[1][i].up= c*h/Keff;
        Keff = (2*mpMat[1][i].perm*mpMat[1][i-1].perm)/(mpMat[1][i].perm + mpMat[1][i-1].perm);
        mbeta[1][i].dn= c*h/Keff;
        Keff = (2*mpMat[1][i].perm*mpMat[2][i].perm)/(mpMat[1][i].perm + mpMat[2][i].perm);
        mbeta[1][i].rh= c*h/Keff;


        Keff = (2*mpMat[n][i].perm*mpMat[n][i+1].perm)/(mpMat[n][i].perm + mpMat[n][i+1].perm);
        mbeta[n][i].up= c*h/Keff;
        Keff = (2*mpMat[n][i].perm*mpMat[n][i-1].perm)/(mpMat[n][i].perm + mpMat[n][i-1].perm);
        mbeta[n][i].dn= c*h/Keff;
        Keff = (2*mpMat[n][i].perm*mpMat[n-1][i].perm)/(mpMat[n][i].perm + mpMat[n-1][i].perm);
        mbeta[n][i].lf= c*h/Keff;

    }


    for(i=2; i<n; i++)
        for(j=2; j<n; j++)
        {
            Keff = (2*mpMat[i][j].perm*mpMat[i][j+1].perm)/(mpMat[i][j].perm + mpMat[i][j+1].perm);
            mbeta[i][j].up= c*h/Keff;
            Keff = (2*mpMat[i][j].perm*mpMat[i][j-1].perm)/(mpMat[i][j].perm + mpMat[i][j-1].perm);
            mbeta[i][j].dn= c*h/Keff;
            Keff = (2*mpMat[i][j].perm*mpMat[i+1][j].perm)/(mpMat[i][j].perm + mpMat[i+1][j].perm);
            mbeta[i][j].rh= c*h/Keff;
            Keff = (2*mpMat[i][j].perm*mpMat[i-1][j].perm)/(mpMat[i][j].perm + mpMat[i-1][j].perm);
            mbeta[i][j].lf= c*h/Keff;
        }

    /*
     * Inicializando valores da fonte
     */
    mpMat[1][1].f=1.0e-7;
    mpMat[n][n].f=-1.0e-7;

    /*
     * calculo de parâmetros que nao dependem das iterações
     */

    aux = 1/h;

    for (i=1; i<=n; i++)
        for (j=1; j<=n; j++)
        {
            mpMat[i][j].shi=2*mpMat[i][j].perm*aux;
            mpMat[i][j].f*=h;
        }

    /*
     * Ciclo até convergência do problema
     */
    k = 0;   // quantidade de iterações

    do
    {
        mp_aux = mp;
        mp = mp_old;
        mp_old = mp_aux;

        mq_aux = mq;
        mq = mq_old;
        mq_old = mq_aux;

        mq_aux = ml;
        ml = ml_old;
        ml_old = mq_aux;

        k++;
        //printf("Iteração %d \n",k);

        /*Cálculo da pressão e dos fluxos em cada elemento */

        /*Canto inferior esquerdo [1][1]*/
        canto_d_l(1, 1, mpMat, mbeta, mq, mq_old, ml_old, mp);

        /*Canto superior esquerdo [1][N-2]*/
        canto_u_l(1, n, mpMat, mbeta, mq, mq_old, ml_old, mp);

        /*Canto inferior direito [N-2][1]*/
        canto_d_r(n, 1, mpMat, mbeta, mq, mq_old, ml_old, mp);

        /*Canto superior direito [N-2][N-2]*/
        canto_u_r(n, n, mpMat, mbeta, mq, mq_old, ml_old, mp);

        /*Fronteira U [2...N-3][N-2]*/
        for (i=2; i<n; i++)
            fronteira_u(i, n, mpMat, mbeta, mq, mq_old, ml_old, mp);

        /*Fronteira D [2...N-3][1]*/
        for (i=2; i<n; i++)
            fronteira_d(i, 1, mpMat, mbeta, mq, mq_old, ml_old, mp);

        /*Fronteira R [N-2][2...N-3]*/
        for (j=2; j<n; j++)
            fronteira_r(n, j, mpMat, mbeta, mq, mq_old, ml_old, mp);

        /*Fronteira L [1][2...N-3]*/
        for (j=2; j<n; j++)
            fronteira_l(1, j, mpMat, mbeta, mq, mq_old, ml_old, mp);

        /*Elementos internos [2..N-3][2..N-3]*/
        for (i=2; i<n; i++)
            for (j=2; j<n; j++)
                internos(i, j, mpMat, mbeta, mq, mq_old, ml_old, mp);

        /*
        	 * Atualização dos multiplicadores de lagrange
        	 */

        Media = 0.0;

        for (i=1; i<=n; i++)
            for (j=1; j<=n; j++)
            {
                ml[i][j].up = mbeta[i][j].up*(mq[i][j].up + mq_old[i][j+1].dn) + ml_old[i][j+1].dn;
                ml[i][j].dn = mbeta[i][j].dn*(mq[i][j].dn + mq_old[i][j-1].up) + ml_old[i][j-1].up;
                ml[i][j].rh = mbeta[i][j].rh*(mq[i][j].rh + mq_old[i+1][j].lf) + ml_old[i+1][j].lf;
                ml[i][j].lf = mbeta[i][j].lf*(mq[i][j].lf + mq_old[i-1][j].rh) + ml_old[i-1][j].rh;
                //printf("P[%d][%d]=%f\n",i,j,mp[i][j]);
                Media += mp[i][j];
            }

        //printf("M = %f\n", Media);

        Media /= (n*n);

        //printf("M = %f\n", Media);

        sum1 = 0.;
        sum2 = 0.;

        /* Impondo a média zero na distriubição de pressões
         * Início de cálculo de verificação de convergência
         */
        for (i=1; i<=n; i++)
            for (j=1; j<=n; j++)
            {
                mp[i][j] -= Media;
                ml[i][j].up -= Media;
                ml[i][j].dn -= Media;
                ml[i][j].rh -= Media;
                ml[i][j].lf -= Media;

                aux = mp[i][j] - mp_old[i][j];
                sum1 += aux*aux;
                sum2 += mp[i][j] * mp[i][j];
            }

        /*Erro relativo entre a pressão atual e anterior*/
        erro = sqrt(sum1/sum2);

    }while(erro > 1e-5);

    free(mpMat);
    free(pMat);
    free(mp_old);
    free(p_old);
    free(mp);
    free(p);
    free(mbeta);
    free(beta);
    free(ml_old);
    free(l_old);
    free(ml);
    free(l);
    free(mq_old);
    free(q_old);
    free(mq);
    free(q);

    ticks2 = clock();
    elapsedTime = walltime( &startTime );
    stop = omp_get_wtime();

    /*Dados finais mostrados na tela*/

    printf("\n------ MALHA %d x %d",n,n);
    printf("\n------ PRESSAO ------- ITERACAO %d ------ ERRO %6.4E", k, erro);
    printf("\n------ CPU TIME %ld, %ld, %g,  %g, %.4g", ticks2, ticks1, (ticks1+1.-1.)/CLOCKS_PER_SEC, (ticks2+1.-1.)/CLOCKS_PER_SEC, (ticks2+1.-1.)/CLOCKS_PER_SEC - (ticks1+1.-1.)/CLOCKS_PER_SEC);
    printf("\n------ TOTAL TIME  %.6f\n\n",elapsedTime);
    printf("\n------ TOTAL TIME  %.6f\n\n",stop - start);

    return 0; /* fim da funcao main */


    printf("Oi\n");

    return 0;


}
