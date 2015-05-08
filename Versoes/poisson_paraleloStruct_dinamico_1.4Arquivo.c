#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

#define PamM 2e-11
#define S 0.5


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
// Todos os cantos agora serão calculados em apenas uma chamada de função
void canto_d_l( const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void canto_u_l( const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void canto_d_r( const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void canto_u_r( const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void fronteira_l( const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void fronteira_r( const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void fronteira_u( const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void fronteira_d( const int, const int, const int, Fluxos ,  Multiplicadores ,  Fisicos);
void internos( const int _j, const int _k, const int n, Fluxos f,  Multiplicadores m, Fisicos fis);


int main(int argc, char *argv[])
{
	double size=25600.00; /* dimensão da regiao /*/
    // Declaração das estruturas que substituiram as variaveis globais.

     Fluxos f; // Declaração dos Fluxos de cada celula
     Multiplicadores m;//Declaração dos Multiplicadores de lagrange
     Fisicos fis; // Declaração das variaveis fisicas, como permeabilidade, os betas e etc
    /* register */double *aux_reg;
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

	double *p_aux, **mat;
	float t;
	FILE *arq;
	FILE *ref;

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
	//--------------------------------------------------------------------------------
    mat = (double**) malloc(N*sizeof(double*));

    int Imat;

	if(mat==NULL)
	{
		printf("Erro na alocacao de memoria\n");
	}

	for(Imat=0;Imat<N;Imat++)
	{
		mat[Imat] = (double*) malloc(N*sizeof(double));

		if(mat[Imat]==NULL)
		{
			printf("Erro na alocacao de memoria\n");
		}

	}
	//--------------------------------------------------------------------------------


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


	ref = fopen("comp128x.txt", "w");
	arq = fopen("perm.dat", "r");

    	for (k=1; k<=n; k++)
			{
				fscanf(arq, "%d", &i);
				for (j=1; j<=n; j++)
					{
						fscanf(arq,"%g", &t);
						fis.perm[j*n + k] = PamM*exp(S*t);
						//printf("%d %d %e\n",k,j,t);
		//				printf("%d %d :%e\n",k,j,fis.perm[j*n + k]);
						//printf("%e \n", fis.perm[n+1]);
				  }
				fscanf(arq, "%d \n", &i);
			}

		fclose(arq);

	/*
	 * Calcula os Beta da Condicao de Robin
  	 */

	Keff = (2*fis.perm[n+1]*fis.perm[n+2])/(fis.perm[n+2]+fis.perm[n+1]);
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
		Keff = (2*fis.perm[(n)+i]*fis.perm[(n)+(i-1)])/(fis.perm[(n)+i]+fis.perm[(n)+(i-1)]);
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
	  		//printf("1 %d %d %e \n", j, k, fis.shi[soma]);
			fis.f[soma]*=h;
			//printf("2 %d %d %e \n",j,k,fis.f[soma]);
		}
	}
	/*
    * Ciclo até convergência do problema
    */
  	i = 1;
  	int soma1;


	for (; ;)
	{

	  	/*Cálculo da pressão e dos Fluxo sem cada elemento */

		/*Canto inferior esquerdo [1][1]*/
		canto_d_l(1, 1,n,f,m,fis);

		/*Canto superior esquerdo [1][N-2]*/
		canto_u_l( 1, n,n,f,m,fis);

		/*Canto inferior direito [N-2][1]*/
		canto_d_r( n, 1,n,f,m,fis);

		/*Canto superior direito [N-2][N-2]*/
		canto_u_r( n, n,n,f,m,fis);

		M = 0.0;
		sum1 = 0.;
		sum2 = 0.0;

		for (j=2; j<n; j++){
         /*Fronteira U [2...N-3][N-2]*/
			fronteira_u( j, n,n,f,m,fis);
		/*Fronteira D [2...N-3][1]*/
			fronteira_d( j, 1,n,f,m,fis);
		/*Fronteira R [N-2][2...N-3]*/
			fronteira_r( n, j,n,f,m,fis);
		/*Fronteira L [1][2...N-3]*/
			fronteira_l( 1, j,n,f,m,fis);

        }

        for (j=2; j<n; j++){
            /*Elementos internos [2..N-3][2..N-3]*/
			for (k=2; k<n; k++){
					internos( j, k, n, f,  m, fis);

			}
		}


			for(j=1;j<=n;j++){
				aux_reg = &fis.p[j*n];
				for (k=1;k<=n; k++)
				M += *(aux_reg+k);
			}
			/*
			 * Atualização dos multiplicadores de lagrange
			 */


			M = M / (n*n);

			for (j=1; j<=n; j++){
			  lin = j*n;
			    register double
				*betauReg = &fis.betau[lin],*betadReg = &fis.betad[lin],*betarReg = &fis.betar[lin],*betalReg = &fis.betal[lin],
				*luReg = &m.lu[lin],*ldReg = &m.ld[lin],*lrReg = &m.lr[lin],*llReg = &m.ll[lin],
				*pReg = &fis.p[lin];

				register double
				*quOldReg = &f.qu_old[lin],*qdOldReg = &f.qd_old[lin],*qrOldReg = &f.qr_old[(j-1)*n],*qlOldReg = &f.ql_old[(j+1)*n],
				*luOldReg = &m.lu_old[lin],*ldOldReg = &m.ld_old[lin],*lrOldReg = &m.lr_old[(j-1)*n],*llOldReg = &m.ll_old[(j+1)*n],
				*quReg = &f.qu[lin],*qdReg = &f.qd[lin],*qrReg = &f.qr[lin],*qlReg = &f.ql[lin];

				for (k=1; k<=n; k++)
			  	{
				  	soma = lin + k;
					*(luReg+k) = *(betauReg+k)*(*(quReg+k) + *(qdOldReg+k+1)) + *(ldOldReg+k+1);
					*(ldReg+k) = *(betadReg+k)*(*(qdReg+k) + *(quOldReg+k-1)) + *(luOldReg+k-1);
					*(lrReg+k) = *(betarReg+k)*(*(qrReg+k) + *(qlOldReg+k)) + *(llOldReg+k);
					*(llReg+k) = *(betalReg+k)*(*(qlReg+k) + *(qrOldReg+k)) + *(lrOldReg+k);

					/* Impondo a média zero na distriubição de pressões
			 		* Início de cálculo de verificação de convergência
			 		*/
					*(pReg+k) -= M;
			  		*(luReg+k) -= M;
			  		*(ldReg+k) -= M;
			  		*(lrReg+k) -= M;
			  		*(llReg+k) -= M;

					aux = *(pReg+k) - fis.p_old[soma];
			  		sum1 += aux*aux;
			  		sum2 += *(pReg+k)**(pReg+k);

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

	} /* fim do ciclo "infinito" for(;;) */

	for(k=1 ;k<=n; k++)
    {
        fprintf(ref,"%d ",k);
        for(j=1; j<=n; j++)
            fprintf(ref,"%e ",fis.p[j*n+k]);
        fprintf(ref,"%d \n",192837465);
    }

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
void canto_d_l( const int _j, const int _k,const int n,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{
	int lin = _j*n;
	register double AuxU,AuxR,DU,DR;	/* Auxiliares para cada lado das células */
	register double *shiReg = &fis.shi[lin],
	*betauReg = &fis.betau[lin],*betadReg = &fis.betad[lin],*betarReg = &fis.betar[lin],*betalReg = &fis.betal[lin],
	*pReg = &fis.p[lin];

	double *fReg = &fis.f[lin],
	*qdOldReg = &f.qd_old[lin],*qlOldReg = &f.ql_old[(_j+1)*n],
	*ldOldReg = &m.ld_old[lin],*llOldReg = &m.ll_old[(_j+1)*n],
	*quReg = &f.qu[lin],*qdReg = &f.qd[lin],*qrReg = &f.qr[lin],*qlReg = &f.ql[lin],
	*luReg = &m.lu[lin],*ldReg = &m.ld[lin],*lrReg = &m.lr[lin],*llReg = &m.ll[lin];

	AuxU = *(shiReg + _k)/(1+*(betauReg + _k)* *(shiReg + _k));
   	AuxR = *(shiReg + _k)/(1+*(betarReg + _k)	* *(shiReg + _k));
   	DU = AuxU*(*(betauReg + _k)* *(qdOldReg + _k+1)+ *(ldOldReg + _k+1));
   	DR = AuxR*(*(betarReg + _k)* *(qlOldReg + _k)+ *(llOldReg+_k));
   	*(pReg + _k ) = (*(fReg + _k) + DU + DR)/(AuxU + AuxR);
   	*(quReg+_k) = AuxU**(pReg + _k ) - DU;
   	*(qrReg+_k) = AuxR**(pReg + _k ) - DR;



}

/* Função para o canto superior esquerdo*/
void canto_u_l( const int _j, const int _k, const int n,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{
	int lin = _j * n;
	register double AuxD, AuxR,DD, DR;	/* Auxiliares para cada lado das células */

	register double *shiReg = &fis.shi[lin],
	*betarReg = &fis.betar[lin],*betadReg = &fis.betad[lin],
	*pReg = &fis.p[lin];

	double *fReg = &fis.f[lin],
	*qlOldReg = &f.ql_old[(_j+1)*n],*quOldReg = &f.qu_old[(_j-1)*n],
	*luOldReg = &m.lu_old[lin],*llOldReg = &m.ll_old[(_j+1)*n],
	*qdReg = &f.qd[lin],*qrReg = &f.qr[lin];

	AuxD = *(shiReg + _k )/(1+ *(betadReg+_k)* *(shiReg + _k ));
   	AuxR = *(shiReg + _k )/(1+ *(betarReg+_k)* *(shiReg + _k ));
   	DD = AuxD*(*(betadReg +_k)**(quOldReg+_k)+*(luOldReg + _k-1));
   	DR = AuxR*(*(betarReg)**(qlOldReg+_k)+ *(llOldReg + _k));
   	*(pReg + _k ) = ((*(fReg + _k )) + DD + DR)/(AuxD + AuxR);
   	*(qdReg+_k) = AuxD * *(pReg + _k ) - DD;
   	*(qrReg+_k) = AuxR * *(pReg + _k ) - DR;


}

/* Função para o canto inferior direito*/
void canto_d_r( const int _j, const int _k,const int n, Fluxos f,  Multiplicadores m,  Fisicos fis)
{
	int lin = _j*n;
	register double AuxU, AuxL,DU,DL;	/* Auxiliares para cada lado das células */

	register double *shiReg = &fis.shi[lin],
	*betauReg = &fis.betau[lin],*betalReg = &fis.betal[lin],
	*pReg = &fis.p[lin];

	double *fReg = &fis.f[lin],
	*qdOldReg = &f.qd_old[lin],*qrOldReg = &f.qr_old[(_j-1)*n],
	*ldOldReg = &m.ld_old[lin],*lrOldReg = &m.lr_old[(_j-1)*n],
	*quReg = &f.qu[lin],*qlReg = &f.ql[lin];

   	AuxU = *(shiReg + _k )/(1+ *(betauReg+_k)**(shiReg + _k ));
   	AuxL = *(shiReg + _k )/(1+ *(betalReg+_k)**(shiReg + _k ));
   	DU = AuxU*(*(betauReg+_k)* *(qdOldReg + _k+1)+*(ldOldReg + _k+1));
   	DL = AuxL*(*(betalReg+_k)* *(qrOldReg +_k)+ *(lrOldReg +_k));
   	*(pReg +_k) = (*(fReg+_k) + DU + DL)/(AuxU + AuxL);
   	*(quReg+_k) = AuxU**(pReg +_k) - DU;
   	*(qlReg+_k) = AuxL**(pReg +_k) - DL;


}

/* Funcao para o canto superior dereito*/
void canto_u_r( const int _j, const int _k, const int n,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{

	register double AuxD,AuxL,DD,DL;	/* Auxiliares para cada lado das células */
	int lin = _j*n;

	register double *shiReg = &fis.shi[lin],
	*betadReg = &fis.betad[lin],*betalReg = &fis.betal[lin],
	*pReg = &fis.p[lin];

	double *fReg = &fis.f[lin],
	*qrOldReg = &f.qr_old[(_j-1)*n],*quOldReg = &f.qu_old[lin],
	*luOldReg = &m.lu_old[lin],*lrOldReg = &m.lr_old[(_j-1)*n],
	*qdReg = &f.qd[lin],*qlReg = &f.ql[lin];

	AuxD = *(shiReg + _k)/(1+*(betadReg + _k)**(shiReg + _k));
   	AuxL = *(shiReg + _k)/(1+*(betalReg + _k)**(shiReg + _k));
   	DD = AuxD*(*(betadReg + _k)**(quOldReg + _k-1)+*(luOldReg + _k-1));
   	DL = AuxL*(*(betalReg + _k)**(qrOldReg + _k)+*(lrOldReg + _k));
   	*(pReg + _k) = (*(fReg + _k) + DD + DL)/(AuxD + AuxL);
   	f.qd[lin + _k] = AuxD**(pReg + _k) - DD;
   	f.ql[lin + _k] = AuxL**(pReg + _k) - DL;



}

/* Função para a fronteira esquerda L */
void fronteira_l( const int _j, const int _k,const int n,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{
	int lin = _j*n;

	register double AuxU, AuxD, AuxR,DU, DD, DR; /* Auxiliares para cada lado das células */

	register double *shiReg = &fis.shi[lin],
	*betauReg = &fis.betau[lin],*betarReg = &fis.betar[lin],*betadReg = &fis.betad[lin],
	*pReg = &fis.p[lin];

	double *fReg = &fis.f[lin],
	*qdOldReg = &f.qd_old[lin],*qlOldReg = &f.ql_old[(_j+1)*n],*quOldReg = &f.qu_old[lin],
	*ldOldReg = &m.ld_old[lin],*luOldReg = &m.lu_old[lin],*llOldReg = &m.ll_old[(_j+1)*n],
	*quReg = &f.qu[lin],*qdReg = &f.qd[lin],*qrReg = &f.qr[lin];

	AuxU = *(shiReg + _k)/(1+*(betauReg + _k)**(shiReg + _k));
  	AuxD = *(shiReg + _k)/(1+*(betadReg + _k)**(shiReg + _k));
  	AuxR = *(shiReg + _k)/(1+*(betarReg + _k)**(shiReg + _k));
  	DU = AuxU*(*(betauReg + _k)**(qdOldReg + _k+1)+*(ldOldReg + _k+1));
  	DD = AuxD*(*(betadReg + _k)**(quOldReg + _k-1)+*(luOldReg + _k-1));
  	DR = AuxR*(*(betarReg + _k)**(qlOldReg + _k)+*(llOldReg + _k));
  	*(pReg + _k) = (*(fReg + _k) + DU + DR + DD)/(AuxU + AuxR + AuxD);
  	f.qu[lin + _k] = AuxU**(pReg + _k) - DU;
  	f.qd[lin + _k] = AuxD**(pReg + _k) - DD;
  	f.qr[lin + _k] = AuxR**(pReg + _k) - DR;



}

/* Função para a fronteira dereita R */
void fronteira_r( const int _j, const int _k,const int n,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{

	 int lin = _j * n;

	register double AuxU,AuxD, AuxL,DU, DD, DL;	/* Auxiliares para cada lado das celulas */

	register double *shiReg = &fis.shi[lin],
	*betauReg = &fis.betau[lin],*betadReg = &fis.betad[lin],*betalReg = &fis.betal[lin],
	*pReg = &fis.p[lin];

	double	*fReg = &fis.f[lin],
	*qdOldReg = &f.qd_old[lin],*qrOldReg = &f.qr_old[((_j-1)*n)],*quOldReg = &f.qu_old[lin],
	*ldOldReg = &m.ld_old[lin],*luOldReg = &m.lu_old[lin],*lrOldReg = &m.lr_old[((_j-1)*n)],
	*quReg = &f.qu[lin],*qdReg = &f.qd[lin],*qlReg = &f.ql[lin],*qrReg = &f.qr[lin];

	AuxU = *(shiReg + _k)/(1+*(betauReg + _k)**(shiReg + _k));
   	AuxD = *(shiReg + _k)/(1+*(betadReg + _k)**(shiReg + _k));
   	AuxL = *(shiReg + _k)/(1+*(betalReg + _k)**(shiReg + _k));
   	DU = AuxU*(*(betauReg + _k)**(qdOldReg + _k+1)+*(ldOldReg + _k+1));
   	DD = AuxD*(*(betadReg + _k)**(quOldReg + _k-1)+*(luOldReg + _k-1));
   	DL = AuxL*(*(betalReg + _k)**(qrOldReg + _k)+*(lrOldReg + _k));
   	*(pReg + _k) = (*(fReg + _k) + DU + DL + DD)/(AuxU + AuxL + AuxD);
   	*(quReg + _k) = AuxU**(pReg + _k) - DU;
   	*(qdReg + _k) = AuxD**(pReg + _k) - DD;
   	*(qlReg + _k) = AuxL**(pReg + _k) - DL;



}

/* Função para a fronteira superior U */
void fronteira_u( const int _j, const int _k,const int n,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{

	int lin = _j * n;

	register double AuxD, AuxR, AuxL, DD, DR, DL;	/* Auxiliares para cada lado das celulas */

	register double *shiReg = &fis.shi[lin],
	*betarReg = &fis.betar[lin],*betadReg = &fis.betad[lin],*betalReg = &fis.betal[lin],
	*pReg = &fis.p[lin];

	double *fReg = &fis.f[lin],
	*qrOldReg = &f.qr_old[(_j-1)*n],*qlOldReg = &f.ql_old[(_j+1)*n],*quOldReg = &f.qu_old[lin],
	*ldOldReg = &m.ld_old[lin],*luOldReg = &m.lu_old[lin],*lrOldReg = &m.lr_old[(_j-1)*n],*llOldReg = &m.ll_old[(_j+1)*n],
	*quReg = &f.qu[lin],*qdReg = &f.qd[lin],*qlReg = &f.ql[lin],*qrReg = &f.qr[lin];

	AuxL = *(shiReg + _k)/(1+*(betalReg + _k)**(shiReg + _k));
  	AuxR = *(shiReg + _k)/(1+*(betarReg + _k)**(shiReg + _k));
  	AuxD = *(shiReg + _k)/(1+*(betadReg + _k)**(shiReg + _k));
  	DL = AuxL*(*(betalReg + _k)**(qrOldReg+_k)+*(lrOldReg+_k));
  	DR = AuxR*(*(betarReg + _k)**(qlOldReg+_k)+*(llOldReg+_k));
  	DD = AuxD*(*(betadReg + _k)**(quOldReg+(_k-1))+*(luOldReg+(_k-1)));
  	*(pReg + _k) = (*(fReg+_k) + DD + DL + DR)/(AuxD + AuxL + AuxR);
  	*(qlReg+_k) = AuxL**(pReg + _k) - DL;
  	*(qrReg+_k) = AuxR**(pReg + _k) - DR;
  	*(qdReg+_k) = AuxD**(pReg + _k) - DD;


}

/* Função para a fronteira inferior D */
void fronteira_d( const int _j, const int _k, const int n, Fluxos f,  Multiplicadores m,  Fisicos fis)
{
	int lin = _j * n;
	register double AuxU,AuxR, AuxL,DU,DR, DL;	/* Auxiliares para cada lado das celulas */

	register double *shiReg = &fis.shi[lin],
	*betauReg = &fis.betau[lin],*betarReg = &fis.betar[lin],*betalReg = &fis.betal[lin],
	*pReg = &fis.p[lin];

	double	*fReg = &fis.f[lin],
	*qdOldReg = &f.qd_old[lin],*qrOldReg = &f.qr_old[(_j-1)*n],*qlOldReg = &f.ql_old[(_j+1)*n],
	*ldOldReg = &m.ld_old[lin],*lrOldReg = &m.lr_old[(_j-1)*n],*llOldReg = &m.ll_old[(_j+1)*n],
	*quReg = &f.qu[lin],*qlReg = &f.ql[lin],*qrReg = &f.qr[lin];

  	AuxL = *(shiReg + _k)/(1+*(betalReg + _k)**(shiReg + _k));
  	AuxR = *(shiReg + _k)/(1+*(betarReg + _k)**(shiReg + _k));
  	AuxU = *(shiReg + _k)/(1+*(betauReg + _k)**(shiReg + _k));
  	DL = AuxL*(*(betalReg + _k)**(qrOldReg+_k)+*(lrOldReg+_k));
  	DR = AuxR*(*(betarReg + _k)**(qlOldReg+_k)+*(llOldReg+_k));
  	DU = AuxU*(*(betauReg + _k)**(qdOldReg+_k+1)+*(ldOldReg+_k+1));
  	*(pReg + _k) = (*(fReg + _k) + DU + DL + DR)/(AuxU + AuxL + AuxR);
  	*(qlReg + _k) = AuxL**(pReg + _k) - DL;
  	*(qrReg + _k) = AuxR**(pReg + _k) - DR;
  	*(quReg + _k) = AuxU**(pReg + _k) - DU;


}


/* Função para os elementos interiores */
void internos( const int _j, const int _k, const int n, Fluxos f,  Multiplicadores m, Fisicos fis)
{
	int lin = _j * n;

	register double
	shi = fis.shi[lin+_k], *pReg = &fis.p[lin];

	register double AuxU, AuxD, AuxR, AuxL,DU, DD, DR, DL;	/* Auxiliares para cada lado das celulas */

	double *qdOldReg = &f.qd_old[lin],*qrOldReg = &f.qr_old[(_j-1)*n],*qlOldReg = &f.ql_old[(_j+1)*n],*quOldReg = &f.qu_old[lin],
	*quReg = &f.qu[lin],*qdReg = &f.qd[lin],*qlReg = &f.ql[lin],*qrReg = &f.qr[lin],
	*ldOldReg = &m.ld_old[lin],*luOldReg = &m.lu_old[lin],*lrOldReg = &m.lr_old[(_j-1)*n],*llOldReg = &m.ll_old[(_j+1)*n],
	*fReg = &fis.f[lin];

	double p,
	betau = fis.betau[lin+_k], betar = fis.betar[lin+_k], betad = fis.betad[lin+_k], betal = fis.betal[lin+_k];

  	AuxL = shi/(1+betal*shi);
  	AuxR = shi/(1+betar*shi);
  	AuxU = shi/(1+betau*shi);
  	AuxD = shi/(1+betad*shi);
  	DL = AuxL*(betal**(qrOldReg+_k)+*(lrOldReg+_k));
  	DR = AuxR*(betar**(qlOldReg+_k)+*(llOldReg+_k));
  	DU = AuxU*(betau**(qdOldReg+_k+1)+*(ldOldReg+_k+1));
  	DD = AuxD*(betad**(quOldReg+_k-1)+*(luOldReg+_k-1));

	*(pReg+_k) = (*(fReg+_k) + DU + DD + DR + DL)/(AuxU + AuxR + AuxD + AuxL);

	p = fis.p[lin+_k];

  	*(qlReg+_k) = AuxL*p - DL;
 	*(qrReg+_k) = AuxR*p - DR;
  	*(quReg+_k) = AuxU*p - DU;
 	*(qdReg+_k) = AuxD*p - DD;


 	//----------------------------------------------------------------------------
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


/*register double *shiReg = &fis.shi[lin],
	*betauReg = &fis.betau[lin],*betarReg = &fis.betar[lin],*betadReg = &fis.betad[lin],*betalReg = &fis.betal[lin],
	*qdOldReg = &f.qd_old[lin],*qrOldReg = &f.qr_old[lin],*qlOldReg = &f.ql_old[(_j+1)*n],*quOldReg = &f.qu_old[lin],
	*ldOldReg = &m.ld_old[lin],*luOldReg = &m.lu_old[lin],*lrOldReg = &m.lr_old[lin],*llOldReg = &m.ll_old[(_j+1)*n],
	*quReg = &f.qu[lin],*qdReg = &f.qd[lin],*qlReg = &f.ql[lin],*qrReg = &f.qr[lin],
	*pReg = &fis.p[lin],*fReg = &fis.f[lin];*/
