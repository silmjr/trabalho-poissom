#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

 typedef struct  {
    /*
     * Fluxos atuais e antigos em cada um dos lados da célula espacial
     * U (Uper), D (Down), L (Left), R (Right)
    */
      double *qu;
      double *qu_old;
      double *qd;
      double *qd_old;
      double *qr;
      double *qr_old;
      double *ql;
      double *ql_old;
}Fluxos;

 typedef struct {
    /*
      *  Multiplicadores de Lagrange em cada um dos lados da célula espacial
    */
      double *lu;
      double *lu_old;
      double *ld;
      double *ld_old;
      double *lr;
      double *lr_old;
      double *ll;
      double *ll_old;
}Multiplicadores;

 typedef struct  {
    /*
     *	Pressões atuais e antigas em cada uma das células
     */
     double *p;
     double *p_old;

    /*
     *	Betas da condição de Robin em cada um dos lados da célula espacial
     */
     double *betau;
     double *betad;
     double *betar;
     double *betal;


     double *f;     /* fonte */
     double *perm;  /* permeabilidade das rochas */
     double *shi;	 /* variável auxiliar valor de shi */


}Fisicos;


/*
* Medida do tempo de processamento levando em consideração o tempo de relógio
*/
double walltime( double* );

/*
 *	Funcoes para cálculo do método de decomposicão de domínio
 */
/*void canto_d_l(const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void canto_u_l(const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void canto_d_r(const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void canto_u_r(const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);*/
// Todos os cantos agora serão calculados em apenas uma chamada de função
void cantos(const int, const int, const int, Fluxos, Multiplicadores, Fisicos);
void fronteiras(const int, const int, const int, Fluxos, Multiplicadores, Fisicos);
/*void fronteira_l(const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void fronteira_r(const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void fronteira_u(const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void fronteira_d(const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void internos(const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);*/





double size=25600.00; /* dimensão da regiao /*/
int main(int argc, char *argv[])
{

    // Declaração das estruturas que substituiram as variaveis globais.

     Fluxos f; // Declaração dos Fluxos de cada celula
     Multiplicadores m;//Declaração dos Multiplicadores de lagrange
     Fisicos fis; // Declaração das variaveis fisicas, como permeabilidade, os betas e etc

      int N;

   int j,k,i,  /* Variáveis auxiliares na implementação */
	    n;      /* Quantidade real de células espaciais n = N-2 */

	double startTime, elapsedTime; /* Variáveis para auxiliar na contagem do tempo*/
	double clockZero = 0.0;
	clock_t ticks1, ticks2;

   double	aux, h,
          	AuxU, AuxD, AuxR, AuxL, /* Auxiliares para cada lado das celulas */
          	DU, DD, DR, DL,         /* Auxiliares para cada lado das celulas */
            M, erro,                /* Media das pressoes e erro na norma */
            sum1, sum2,             /* Auxiliares somadores */
            c = 1.,                 /* Valor para o calculo de beta */
            Keff;                   /* Valor medio da permeabilidade */

	double *p_aux;

	if(argc != 2)
	{
        printf("Erro na execucao!\nEx. ./[nome_executavel] [tamanho_da_malha+2]\n");
		  return (1);
	}
	else
	{
		N = atoi(argv[1]);
	}

	/* No lugar de chamar a função repetidas vezes achei melhor fazer a alocação explicitamente*/
	/*Cria ponteiro com linhas = N colunas = N*/
	f.qu = (double*) malloc((N*N)*sizeof(double));
	f.qu_old = (double*) malloc((N*N)*sizeof(double));
	f.qd = (double*) malloc((N*N)*sizeof(double));
	f.qd_old = (double*) malloc((N*N)*sizeof(double));
	f.qr = (double*) malloc((N*N)*sizeof(double));
	f.qr_old = (double*) malloc((N*N)*sizeof(double));
	f.ql = (double*) malloc((N*N)*sizeof(double));
	f.ql_old = (double*) malloc((N*N)*sizeof(double));

	m.lu = (double*) malloc((N*N)*sizeof(double));
	m.lu_old = (double*) malloc((N*N)*sizeof(double));
	m.ld = (double*) malloc((N*N)*sizeof(double));
	m.ld_old = (double*) malloc((N*N)*sizeof(double));
	m.lr = (double*) malloc((N*N)*sizeof(double));
	m.lr_old = (double*) malloc((N*N)*sizeof(double));
	m.ll = (double*) malloc((N*N)*sizeof(double));
	m.ll_old = (double*) malloc((N*N)*sizeof(double));

	fis.p = (double*) malloc((N*N)*sizeof(double));
	fis.p_old = (double*) malloc((N*N)*sizeof(double));
	fis.betau = (double*) malloc((N*N)*sizeof(double));
	fis.betad = (double*) malloc((N*N)*sizeof(double));
	fis.betar = (double*) malloc((N*N)*sizeof(double));
	fis.betal = (double*) malloc((N*N)*sizeof(double));
	fis.f = (double*) malloc((N*N)*sizeof(double));
	fis.perm = (double*) malloc((N*N)*sizeof(double));
	fis.shi = (double*) malloc((N*N)*sizeof(double));

	startTime = walltime( &clockZero );
  	ticks1 = clock();
	int lin,soma;

	/*
	* Inicialização das variaveis
	*/
	for(j=1; j<N-1; j++){
		lin = j*(N-2);
		for(k=1; k<N-1; k++)
		{	soma = lin + k;
		 	fis.f[soma]=0.0;
			m.lu_old[soma] = m.ld_old[soma] = m.lr_old[soma] = m.ll_old[soma] = 0.0;
			f.qu_old[soma] = f.qd_old[soma] = f.qr_old[soma] = 0.0;
		 	f.ql_old[soma]=0.0;
			fis.p_old[soma]=0.0;
		 }
	}
	/*
	 * Cálculo de algumas variáveis auxiliares para o valor de permeabilidade
	 * en blocos verticais
	 */
	n = N-2;
	h = size/n;
	i = n/2;

	for(j=1; j<=i; j++)
		for(k=1; k<N-1; k++)
		{
			fis.perm[(j*n) + k] = 1.0e-10;
			fis.perm[(j+i)*n + k] = 1.0e-11;
		}

	/*
	 * Calcula os Beta da Condicao de Robin
  	 */

	Keff = (2*fis.perm[n+1]*fis.perm[n+2])/(fis.perm[n+2]+fis.perm[n+2]);
	fis.betau[n+1] = c*h/Keff;
	Keff = (2*fis.perm[n+1]*fis.perm[(2*n)+1])/(fis.perm[n+1]+fis.perm[(2*n)+1]);
	fis.betar[n+1] = c*h/Keff;

	Keff = (2*fis.perm[n+n]*fis.perm[n+(n-1)])/(fis.perm[n+n]+fis.perm[n+(n-1)]);
	fis.betad[n+n] = c*h/Keff;
	Keff = (2*fis.perm[n+n]*fis.perm[(2*n)+n])/(fis.perm[n+n]+fis.perm[(2*n)+n]);
	fis.betar[n+n] = c*h/Keff;

	Keff = (2*fis.perm[(n*n)+1]*fis.perm[(n*n)+2])/(fis.perm[(n*n)+1]+fis.perm[(n*n)+2]);
	fis.betau[(n*n)+1] = c*h/Keff;
	Keff = (2*fis.perm[(n*n)+1]*fis.perm[((n-1)*n)+1])/(fis.perm[(n*n)+1]+fis.perm[((n-1)*n)+1]);
	fis.betal[(n*n)+1] = c*h/Keff;

	Keff = (2*fis.perm[(n*n)+n]*fis.perm[(n*n)+(n-1)])/(fis.perm[(n*n)+n]+fis.perm[(n*n)+(n-1)]);
	fis.betad[(n*n)+n] = c*h/Keff;
	Keff = (2*fis.perm[(n*n)+n]*fis.perm[((n-1)*n)+n])/(fis.perm[(n*n)+n]+fis.perm[((n-1)*n)+n]);
	fis.betal[(n*n)+n] = c*h/Keff;
    //Auxialiares para calculo de linha
    int lin1, linmenos,linmais;

	for (i=2; i<n; i++)
	{
	    lin1 = i*n;
	    linmenos = (i-1) * n;
	    linmais = (i+1) * n;
		Keff = (2*fis.perm[(lin1)+1]*fis.perm[(lin1)+2])/(fis.perm[(lin1)+1]+fis.perm[(lin1)+2]);
		fis.betau[(lin1)+1] = c*h/Keff;
		Keff = (2*fis.perm[(lin1)+1]*fis.perm[(linmenos)+1])/(fis.perm[(lin1)+1]+fis.perm[(linmenos)+1]);
		fis.betal[(lin1)+1] = c*h/Keff;
		Keff = (2*fis.perm[(lin1)+1]*fis.perm[(linmais)+1])/(fis.perm[lin1+1]+fis.perm[linmais+1]);
		fis.betar[lin1+1] = c*h/Keff;

		Keff = (2*fis.perm[lin1+n]*fis.perm[lin1+(n-1)])/(fis.perm[lin1+n]+fis.perm[lin1+(n-1)]);
		fis.betad[lin1+n] = c*h/Keff;
		Keff = (2*fis.perm[lin1+n]*fis.perm[linmenos+n])/(fis.perm[lin1+n]+fis.perm[linmenos+n]);
		fis.betal[lin1+n] = c*h/Keff;
		Keff = (2*fis.perm[lin1+n]*fis.perm[linmais+n])/(fis.perm[lin1+n]+fis.perm[linmais+n]);
		fis.betar[lin1+n] = c*h/Keff;

		Keff = (2*fis.perm[(n)+i]*fis.perm[(n)+(i+1)])/(fis.perm[(n)+i]+fis.perm[(n)+(i+1)]);
		fis.betau[(n)+(i)] = c*h/Keff;
		Keff = (2*fis.perm[(n)+i]*fis.perm[(n)+(n-1)])/(fis.perm[(n)+i]+fis.perm[(n)+(i-1)]);
		fis.betad[(n)+(i)] = c*h/Keff;
		Keff = (2*fis.perm[(n)+i]*fis.perm[(2*n)+i])/(fis.perm[(n)+i]+fis.perm[(2*n)+i]);
		fis.betar[(n)+i] = c*h/Keff;

		Keff = (2*fis.perm[(n*n)+i]*fis.perm[(n*n)+(i+1)])/(fis.perm[(n*n)+i]+fis.perm[(n*n)+(i+1)]);
		fis.betau[(n*n)+i] = c*h/Keff;
		Keff = (2*fis.perm[(n*n)+i]*fis.perm[(n*n)+(i-1)])/(fis.perm[(n*n)+i]+fis.perm[(n*n)+(i-1)]);
		fis.betad[(n*n)+i] = c*h/Keff;
		Keff = (2*fis.perm[(n*n)+i]*fis.perm[((n-1)*n)+i])/(fis.perm[(n*n)+i]+fis.perm[((n-1)*n)+i]);
		fis.betal[(n*n)+i] = c*h/Keff;
	}


	for(j=2; j<n; j++){
	  lin = n*j;
	  linmais = (j+1) * n ;
	  linmenos = (j-1) * n;
		for(k=2; k<n; k++)
		{
		   soma = lin + k;
		   Keff = (2*fis.perm[soma]*fis.perm[soma+1])/(fis.perm[soma]+fis.perm[soma+1]);
		   fis.betau[soma] = c*h/Keff;
		   Keff = (2*fis.perm[soma]*fis.perm[soma-1])/(fis.perm[soma]+fis.perm[soma-1]);
		   fis.betad[soma] = c*h/Keff;
		   Keff = (2*fis.perm[soma]*fis.perm[(linmais)+k])/(fis.perm[soma]+fis.perm[(linmais)+k]);
		   fis.betar[soma] = c*h/Keff;
		   Keff = (2*fis.perm[soma]*fis.perm[linmenos+ k])/(fis.perm[soma]+fis.perm[(linmenos)+k]);
		   fis.betal[soma] = c*h/Keff;
	    }
	}
	/*
	 * Inicializando valores da fonte
	 */
	fis.f[n+1]=1.0e-7;
	fis.f[(n*n)+n]=-1.0e-7;

	/*
	 * calculo de parâmetros que nao dependem das iterações
	 */

	aux = 1/h;

	for (j=1; j<=n; j++){
	  lin = n*j;
		for (k=1; k<=n; k++)
		{
		  soma = lin + k;
	  		fis.shi[soma]=2*fis.perm[soma]*aux;
			fis.f[soma]*=h;
		}
	}
	/*
    * Ciclo até convergência do problema
    */
  	i = 1;

	for (; ;)
	{


		//Todos os cantos e fornteiras são calculados simultaneos agora
		cantos(n,n,n,f,m,fis);

		for (j=2; j<n; j++)
            fronteiras(j,j,n,f,m,fis);


		/*Elementos internos [2..N-3][2..N-3]*/
		for (j=2; j<n; j++){
		 lin = j*n;
		 linmais = (j+1) * n;
		 linmenos = (j-1) * n;
			for (k=2; k<n; k++)
			{
                soma = lin + k;
		  		AuxL = fis.shi[soma]/(1+fis.betal[soma]*fis.shi[soma]);
	  			AuxR = fis.shi[soma]/(1+fis.betar[soma]*fis.shi[soma]);
	  			AuxU = fis.shi[soma]/(1+fis.betau[soma]*fis.shi[soma]);
				AuxD = fis.shi[soma]/(1+fis.betad[soma]*fis.shi[soma]);
				DL = AuxL*(fis.betal[soma]*f.qr_old[linmenos + k]+m.lr_old[linmenos + k]);
				DR = AuxR*(fis.betar[soma]*f.ql_old[linmais + k]+m.ll_old[((j+1)*n) + k]);
				DU = AuxU*(fis.betau[soma]*f.qd_old[soma+1]+m.ld_old[soma+1]);
				DD = AuxD*(fis.betad[soma]*f.qu_old[soma-1]+m.lu_old[soma-1]);

				fis.p[soma] = (fis.f[soma] + DU + DD + DR + DL)/(AuxU + AuxR + AuxD + AuxL);
				f.ql[soma] = AuxL*fis.p[soma] - DL;
				f.qr[soma] = AuxR*fis.p[soma] - DR;
				f.qu[soma] = AuxU*fis.p[soma] - DU;
				f.qd[soma] = AuxD*fis.p[soma] - DD;
			}

		}
			/*
			 * Atualização dos multiplicadores de lagrange
			 */

			M = 0.0;

			for (j=1; j<=n; j++){
			  lin = j*n;
			  linmais = (j+1)*n;
			  linmenos = (j-1)*n;
				for (k=1; k<=n; k++)
			  	{
				  soma = lin + k;
					m.lu[soma] = fis.betau[soma]*(f.qu[soma] + f.qd_old[soma+1]) + m.ld_old[soma+1];
					m.ld[soma] = fis.betad[soma]*(f.qd[soma] + f.qu_old[soma-1]) + m.lu_old[soma-1];
					m.lr[soma] = fis.betar[soma]*(f.qr[soma] + f.ql_old[linmais + k]) + m.ll_old[linmais + k];
					m.ll[soma] = fis.betal[soma]*(f.ql[soma] + f.qr_old[linmenos + k]) + m.lr_old[linmenos + k];
					M += fis.p[soma];
				}
			}

			M = M / (n*n);

			sum1 = 0.;
			sum2 = 0.;

			/* Impondo a média zero na distriubição de pressões
			 * Início de cálculo de verificação de convergência
			 */
			for (j=1; j<=n; j++){
				lin = n * j;
				for (k=1; k<=n; k++)
		  		{
					soma = lin + k;
					fis.p[soma] -= M;
			  		m.lu[soma] -= M;
			  		m.ld[soma] -= M;
			  		m.lr[soma] -= M;
			  		m.ll[soma] -= M;

					aux = fis.p[soma] - fis.p_old[soma];
			  		sum1 += aux*aux;
			  		sum2 += fis.p[soma]*fis.p[soma];
				}
			}

		/*Erro relativo entre a pressão atual e anterior*/
		erro = sqrt(sum1/sum2);


		/*Critério de convergência satisfeito - fim do cálculo*/
		if (erro < 1e-5) break;

		/*
		 * Critério de convergência não satisfeito atualização das pressões,
		 * Fluxose multiplicadores antigos para uma nova iteração
		 */


		p_aux = fis.p;
		fis.p = fis.p_old;
		fis.p_old = p_aux;

		p_aux = f.qu;
		f.qu = f.qu_old;
		f.qu_old = p_aux;

		p_aux = f.qd;
		f.qd = f.qd_old;
		f.qd_old = p_aux;

		p_aux = f.ql;
		f.ql = f.ql_old;
		f.ql_old = p_aux;

		p_aux = f.qr;
		f.qr = f.qr_old;
		f.qr_old = p_aux;

		p_aux = m.lu;
		m.lu = m.lu_old;
		m.lu_old = p_aux;

		p_aux = m.ld;
		m.ld = m.ld_old;
		m.ld_old = p_aux;

		p_aux = m.ll;
		m.ll = m.ll_old;
		m.ll_old = p_aux;

		p_aux = m.lr;
		m.lr = m.lr_old;
		m.lr_old = p_aux;

		i++;

	} /* fin do ciclo "infinito" for(;;) */

	free(fis.p);
	free(fis.p_old);
	free(f.qu);
	free(f.qu_old);
	free(f.qd);
	free(f.qd_old);
	free(f.qr);
	free(f.qr_old);
	free(f.ql);
	free(f.ql_old);

	free(m.lu);
	free(m.lu_old);
	free(m.ld);
	free(m.ld_old);
	free(m.lr);
	free(m.lr_old);
	free(m.ll);
	free(m.ll_old);

   ticks2 = clock();
  	elapsedTime = walltime( &startTime );

	/*Dados finais mostrados na tela*/

   printf("\n------ MALHA %d x %d",n,n);
	printf("\n------ PRESSAO ------- ITERACAO %d ------ ERRO %6.4E", i, erro);
  	printf("\n------ CPU TIME %ld, %ld, %g,  %g, %.4g", ticks2, ticks1, (ticks1+1.-1.)/CLOCKS_PER_SEC, (ticks2+1.-1.)/CLOCKS_PER_SEC, (ticks2+1.-1.)/CLOCKS_PER_SEC - (ticks1+1.-1.)/CLOCKS_PER_SEC);
  	printf("\n------ TOTAL TIME  %.6f\n\n",elapsedTime);

   return 0; /* fim da funcao main */
}

void cantos(const int _j, const int _k, const int n, Fluxos f, Multiplicadores m, Fisicos fis)
{
        int lin = n*n;
        double AuxU, AuxD,AuxL,AuxR,DL,DD,DU,DR;	/* Auxiliares para cada lado das células */

        /*Canto inferior esquerdo [1][1]*/
		//canto_d_l(1, 1,n,f,m,fis);

        AuxU = fis.shi[n + 1]/(1+fis.betau[n + 1]*fis.shi[n + 1]);
        AuxR = fis.shi[n + 1]/(1+fis.betar[n + 1]*fis.shi[n + 1]);
        DU = AuxU*(fis.betau[n + 1]*f.qd_old[n + 2]+m.ld_old[n + 2]);
        DR = AuxR*(fis.betar[n + 1]*f.ql_old[(2*n) +1]+m.ll_old[(2*n) + 1]);
        fis.p[n + 1] = (fis.f[n + 1] + DU + DR)/(AuxU + AuxR);
        f.qu[n + 1] = AuxU*fis.p[n + 1] - DU;
        f.qr[n + 1] = AuxR*fis.p[n + 1] - DR;

		/*Canto superior esquerdo [1][N-2]*/
		//canto_u_l(1, n,n,f,m,fis);

		AuxD = fis.shi[n+_k]/(1+fis.betad[n+_k]*fis.shi[n+_k]);
        AuxR = fis.shi[n+_k]/(1+fis.betar[n+_k]*fis.shi[n+_k]);
        DD = AuxD*(fis.betad[n+_k]*f.qu_old[n+(_k-1)]+m.lu_old[n+(_k-1)]);
        DR = AuxR*(fis.betar[n+_k]*f.ql_old[(2*n)+_k]+m.ll_old[(2*n)+_k]);
        fis.p[n+_k] = (fis.f[n+_k] + DD + DR)/(AuxD + AuxR);
        f.qd[n+_k] = AuxD * fis.p[n+_k] - DD;
        f.qr[n+_k] = AuxR * fis.p[n+_k] - DR;

		/*Canto inferior direito [N-2][1]*/
		//canto_d_r(n, 1,n,f,m,fis);

		AuxU = fis.shi[lin + 1]/(1+fis.betau[lin + 1]*fis.shi[lin + 1]);
        AuxL = fis.shi[lin + 1]/(1+fis.betal[lin + 1]*fis.shi[lin + 1]);
        DU = AuxU*(fis.betau[lin + 1]*f.qd_old[lin + 2]+m.ld_old[lin + 2]);
        DL = AuxL*(fis.betal[lin + 1]*f.qr_old[((n-1)*n)+1]+m.lr_old[((n-1)*n)+1]);
        fis.p[lin + 1] = (fis.f[lin + 1] + DU + DL)/(AuxU + AuxL);
        f.qu[lin + 1] = AuxU*fis.p[lin + 1] - DU;
        f.ql[lin + 1] = AuxL*fis.p[lin + 1] - DL;

		/*Canto superior direito [N-2][N-2]*/
		//canto_u_r(n, n,n,f,m,fis);

		AuxD = fis.shi[lin + _k]/(1+fis.betad[lin + _k]*fis.shi[lin + _k]);
        AuxL = fis.shi[lin + _k]/(1+fis.betal[lin + _k]*fis.shi[lin + _k]);
        DD = AuxD*(fis.betad[lin + _k]*f.qu_old[lin + (_k-1)]+m.lu_old[lin + (_k-1)]);
        DL = AuxL*(fis.betal[lin + _k]*f.qr_old[((n-1)*n) + _k]+m.lr_old[((n-1)*n)+_k]);
        fis.p[lin + _k] = (fis.f[lin + _k] + DD + DL)/(AuxD + AuxL);
        f.qd[lin + _k] = AuxD*fis.p[lin + _k] - DD;
        f.ql[lin + _k] = AuxL*fis.p[lin + _k] - DL;
}

void fronteiras(const _j, const _k, const int n, Fluxos f, Multiplicadores m, Fisicos fis)
{
        int lin = _j*n;
        int nn = n*n;
        double AuxL, AuxU, AuxD, AuxR, DL, DU, DD, DR; /* Auxiliares para cada lado das células */

        /*Fronteira U [2...N-3][N-2]*/
        /*fronteira_u(j, n,n,f,m,fis);*/
        AuxL = fis.shi[lin+n]/(1+fis.betal[lin+n]*fis.shi[lin+n]);
        AuxR = fis.shi[lin+n]/(1+fis.betar[lin+n]*fis.shi[lin+n]);
        AuxD = fis.shi[lin+n]/(1+fis.betad[lin+n]*fis.shi[lin+n]);
        DL = AuxL*(fis.betal[lin+n]*f.qr_old[((_j-1)*n)+n]+m.lr_old[((_j-1)*n)+n]);
        DR = AuxR*(fis.betar[lin+n]*f.ql_old[((_j+1)*n)+n]+m.ll_old[((_j+1)*n)+n]);
        DD = AuxD*(fis.betad[lin+n]*f.qu_old[lin+n-1]+m.lu_old[lin+n-1]);
        fis.p[lin+n] = (fis.f[lin+n] + DD + DL + DR)/(AuxD + AuxL + AuxR);
        f.ql[lin+n] = AuxL*fis.p[lin+n] - DL;
        f.qr[lin+n] = AuxR*fis.p[lin+n] - DR;
        f.qd[lin+n] = AuxD*fis.p[lin+n] - DD;

		/*Fronteira D [2...N-3][1]*/
        //fronteira_d(j, 1,n,f,m,fis);
        AuxL = fis.shi[lin+1]/(1+fis.betal[lin+1]*fis.shi[lin+1]);
        AuxR = fis.shi[lin+1]/(1+fis.betar[lin+1]*fis.shi[lin+1]);
        AuxU = fis.shi[lin+1]/(1+fis.betau[lin+1]*fis.shi[lin+1]);
        DL = AuxL*(fis.betal[lin+1]*f.qr_old[((_j-1)*n)+1]+m.lr_old[((_j-1)*n)+1]);
        DR = AuxR*(fis.betar[lin+1]*f.ql_old[((_j+1)*n)+1]+m.ll_old[((_j+1)*n)+1]);
        DU = AuxU*(fis.betau[lin+1]*f.qd_old[lin+2]+m.ld_old[lin+2]);
        fis.p[lin+1] = (fis.f[lin+1] + DU + DL + DR)/(AuxU + AuxL + AuxR);
        f.ql[lin+1] = AuxL*fis.p[lin+1] - DL;
        f.qr[lin+1] = AuxR*fis.p[lin+1] - DR;
        f.qu[lin+1] = AuxU*fis.p[lin+1] - DU;

		/*Fronteira R [N-2][2...N-3]*/
		//fronteira_r(n, j,n,f,m,fis);
		AuxU = fis.shi[nn+_k]/(1+fis.betau[nn+_k]*fis.shi[nn+_k]);
        AuxD = fis.shi[nn+_k]/(1+fis.betad[nn+_k]*fis.shi[nn+_k]);
        AuxL = fis.shi[nn+_k]/(1+fis.betal[nn+_k]*fis.shi[nn+_k]);
        DU = AuxU*(fis.betau[nn+_k]*f.qd_old[nn+_k+1]+m.ld_old[nn+_k+1]);
        DD = AuxD*(fis.betad[nn+_k]*f.qu_old[nn+_k-1]+m.lu_old[nn+_k-1]);
        DL = AuxL*(fis.betal[nn+_k]*f.qr_old[((n-1)*n) +_k]+m.lr_old[((n-1)*n)+_k]);
        fis.p[nn+_k] = (fis.f[nn+_k] + DU + DL + DD)/(AuxU + AuxL + AuxD);
        f.qu[nn+_k] = AuxU*fis.p[nn+_k] - DU;
        f.qd[nn+_k] = AuxD*fis.p[nn+_k] - DD;
        f.ql[nn+_k] = AuxL*fis.p[nn+_k] - DL;
		/*Fronteira L [1][2...N-3]*/
        //fronteira_l(1, j,n,f,m,fis);

        AuxU = fis.shi[n + _k]/(1+fis.betau[n + _k]*fis.shi[n + _k]);
        AuxD = fis.shi[n + _k]/(1+fis.betad[n + _k]*fis.shi[n + _k]);
        AuxR = fis.shi[n + _k]/(1+fis.betar[n + _k]*fis.shi[n + _k]);
        DU = AuxU*(fis.betau[n + _k]*f.qd_old[n + _k+1]+m.ld_old[n + _k+1]);
        DD = AuxD*(fis.betad[n + _k]*f.qu_old[n + _k-1]+m.lu_old[n + _k-1]);
        DR = AuxR*(fis.betar[n + _k]*f.ql_old[(2*n) + _k]+m.ll_old[((2)*n) + _k]);
        fis.p[n + _k] = (fis.f[n + _k] + DU + DR + DD)/(AuxU + AuxR + AuxD);
        f.qu[n + _k] = AuxU*fis.p[n + _k] - DU;
        f.qd[n + _k] = AuxD*fis.p[n + _k] - DD;
        f.qr[n + _k] = AuxR*fis.p[n + _k] - DR;
}


/* Função para os elementos interiores */
void internos(const int _j, const int _k, const int n, Fluxos f,  Multiplicadores m, Fisicos fis)
{
	int lin = _j * n;
	double AuxU, AuxD, AuxR, AuxL,DU, DD, DR, DL;	/* Auxiliares para cada lado das celulas */

  	AuxL = fis.shi[lin+_k]/(1+fis.betal[lin+_k]*fis.shi[lin+_k]);
  	AuxR = fis.shi[lin+_k]/(1+fis.betar[lin+_k]*fis.shi[lin+_k]);
  	AuxU = fis.shi[lin+_k]/(1+fis.betau[lin+_k]*fis.shi[lin+_k]);
  	AuxD = fis.shi[lin+_k]/(1+fis.betad[lin+_k]*fis.shi[lin+_k]);
  	DL = AuxL*(fis.betal[lin+_k]*f.qr_old[((_j-1)*n)+_k]+m.lr_old[((_j-1)*n)+_k]);
  	DR = AuxR*(fis.betar[lin+_k]*f.ql_old[((_j+1)*n)+_k]+m.ll_old[((_j+1)*n)+_k]);
  	DU = AuxU*(fis.betau[lin+_k]*f.qd_old[lin+_k+1]+m.ld_old[lin+_k+1]);
  	DD = AuxD*(fis.betad[lin+_k]*f.qu_old[lin+_k-1]+m.lu_old[lin+_k-1]);

	fis.p[lin+_k] = (fis.f[lin+_k] + DU + DD + DR + DL)/(AuxU + AuxR + AuxD + AuxL);
  	f.ql[lin+_k] = AuxL*fis.p[lin+_k] - DL;
 	f.qr[lin+_k] = AuxR*fis.p[lin+_k] - DR;
  	f.qu[lin+_k] = AuxU*fis.p[lin+_k] - DU;
 	f.qd[lin+_k] = AuxD*fis.p[lin+_k] - DD;

}

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


