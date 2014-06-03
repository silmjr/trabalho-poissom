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
      double **qu;
      double **qu_old;
      double **qd;
      double **qd_old;
      double **qr;
      double **qr_old;
      double **ql;
      double **ql_old;
}Fluxos;

 typedef struct {
    /*
      *  Multiplicadores de Lagrange em cada um dos lados da célula espacial
    */        
      double **lu;
      double **lu_old;
      double **ld;
      double **ld_old;
      double **lr;
      double **lr_old;
      double **ll;
      double **ll_old;
}Multiplicadores;

 typedef struct  {
    /*
     *	Pressões atuais e antigas em cada uma das células
     */ 
     double **p;
     double **p_old;

    /*
     *	Betas da condição de Robin em cada um dos lados da célula espacial
     */        
     double **betau;
     double **betad;
     double **betar;
     double **betal;

       
     double **f;     /* fonte */
     double **perm;  /* permeabilidade das rochas */
     double **shi;	 /* variável auxiliar valor de shi */

  
}Fisicos;
       

/*
* Medida do tempo de processamento levando em consideração o tempo de relógio
*/
double walltime( double* ); 

double **cria_ponteiro(int,int);

/*
 *	Funcoes para cálculo do método de decomposicão de domínio
 */      
void canto_d_l(const int, const int,  Fluxos ,  Multiplicadores ,  Fisicos);
void canto_u_l(const int, const int,  Fluxos ,  Multiplicadores ,  Fisicos);
void canto_d_r(const int, const int,  Fluxos ,  Multiplicadores ,  Fisicos);
void canto_u_r(const int, const int,  Fluxos ,  Multiplicadores ,  Fisicos);
void fronteira_l(const int, const int,  Fluxos ,  Multiplicadores ,  Fisicos);
void fronteira_r(const int, const int,  Fluxos ,  Multiplicadores ,  Fisicos);
void fronteira_u(const int, const int,  Fluxos ,  Multiplicadores ,  Fisicos);
void fronteira_d(const int, const int,  Fluxos ,  Multiplicadores ,  Fisicos);
void internos(const int, const int,  Fluxos ,  Multiplicadores ,  Fisicos);






int main(int argc, char *argv[])
{
    double size=25600.00; /* dimensão da regiao */
    // Declaração das estruturas que substituiram as variaveis globais. 
    
     Fluxos f;
     Multiplicadores m;
     Fisicos fis;
	
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

	double **p_aux;

	if(argc != 2) 
	{
        printf("Erro na execucao!\nEx. ./[nome_executavel] [tamanho_da_malha+2]\n");
		  return (1);
	}
	else
	{
		N = atoi(argv[1]);
	}

	/*Cria ponteiro com linhas = N colunas = N*/
	f.qu = cria_ponteiro(N,N);
	f.qu_old = cria_ponteiro(N,N);
	f.qd = cria_ponteiro(N,N);
	f.qd_old = cria_ponteiro(N,N);
	f.qr = cria_ponteiro(N,N);
	f.qr_old = cria_ponteiro(N,N);
	f.ql = cria_ponteiro(N,N);
	f.ql_old = cria_ponteiro(N,N);

	m.lu = cria_ponteiro(N,N);
	m.lu_old = cria_ponteiro(N,N);
	m.ld = cria_ponteiro(N,N);
	m.ld_old = cria_ponteiro(N,N);
	m.lr = cria_ponteiro(N,N);
	m.lr_old = cria_ponteiro(N,N);
	m.ll = cria_ponteiro(N,N);
	m.ll_old = cria_ponteiro(N,N);

	fis.p = cria_ponteiro(N,N);
	fis.p_old = cria_ponteiro(N,N);
	fis.betau = cria_ponteiro(N,N);
	fis.betad = cria_ponteiro(N,N);
	fis.betar = cria_ponteiro(N,N);
	fis.betal = cria_ponteiro(N,N);
	fis.f = cria_ponteiro(N,N);
	fis.perm = cria_ponteiro(N,N);
	fis.shi = cria_ponteiro(N,N);

	startTime = walltime( &clockZero );
  	ticks1 = clock(); 

	/*
	* Inicialização das variaveis
	*/
	for(j=1; j<N-1; j++)
		for(k=1; k<N-1; k++)
		{
		 	fis.f[j][k]=0.0;
	      m.lu_old[j][k] = m.ld_old[j][k] = m.lr_old[j][k] = m.ll_old[j][k] = 0.0;
	      f.qu_old[j][k] = f.qd_old[j][k] = f.qr_old[j][k] = 0.0;
		 	f.ql_old[j][k]=0.0;
	      fis.p_old[j][k]=0.0;           
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
			fis.perm[j][k] = 1.0e-10;     
			fis.perm[j+i][k] = 1.0e-11;
		}     
		     
	/*
	 * Calcula os Beta da Condicao de Robin
  	 */ 

	Keff = (2*fis.perm[1][1]*fis.perm[1][2])/(fis.perm[1][2]+fis.perm[1][2]);
	fis.betau[1][1] = c*h/Keff;
	Keff = (2*fis.perm[1][1]*fis.perm[2][1])/(fis.perm[1][1]+fis.perm[2][1]);
	fis.betar[1][1] = c*h/Keff;
		 
	Keff = (2*fis.perm[1][n]*fis.perm[1][n-1])/(fis.perm[1][n]+fis.perm[1][n-1]);
	fis.betad[1][n] = c*h/Keff;
	Keff = (2*fis.perm[1][n]*fis.perm[2][n])/(fis.perm[1][n]+fis.perm[2][n]);
	fis.betar[1][n] = c*h/Keff;
	 
	Keff = (2*fis.perm[n][1]*fis.perm[n][2])/(fis.perm[n][1]+fis.perm[n][2]);
	fis.betau[n][1] = c*h/Keff;
	Keff = (2*fis.perm[n][1]*fis.perm[n-1][1])/(fis.perm[n][1]+fis.perm[n-1][1]);
	fis.betal[n][1] = c*h/Keff;
	
	Keff = (2*fis.perm[n][n]*fis.perm[n][n-1])/(fis.perm[n][n]+fis.perm[n][n-1]);
	fis.betad[n][n] = c*h/Keff;
	Keff = (2*fis.perm[n][n]*fis.perm[n-1][n])/(fis.perm[n][n]+fis.perm[n-1][n]);
	fis.betal[n][n] = c*h/Keff;
		
	for (i=2; i<n; i++)
	{
		Keff = (2*fis.perm[i][1]*fis.perm[i][2])/(fis.perm[i][1]+fis.perm[i][2]);
		fis.betau[i][1] = c*h/Keff;
		Keff = (2*fis.perm[i][1]*fis.perm[i-1][1])/(fis.perm[i][1]+fis.perm[i-1][1]);
		fis.betal[i][1] = c*h/Keff;
		Keff = (2*fis.perm[i][1]*fis.perm[i+1][1])/(fis.perm[i][1]+fis.perm[i+1][1]);
		fis.betar[i][1] = c*h/Keff;
		 			
		Keff = (2*fis.perm[i][n]*fis.perm[i][n-1])/(fis.perm[i][n]+fis.perm[i][n-1]);
		fis.betad[i][n] = c*h/Keff;
		Keff = (2*fis.perm[i][n]*fis.perm[i-1][n])/(fis.perm[i][n]+fis.perm[i-1][n]);
		fis.betal[i][n] = c*h/Keff;
		Keff = (2*fis.perm[i][n]*fis.perm[i+1][n])/(fis.perm[i][n]+fis.perm[i+1][n]);
		fis.betar[i][n] = c*h/Keff;
			
		Keff = (2*fis.perm[1][i]*fis.perm[1][i+1])/(fis.perm[1][i]+fis.perm[1][i+1]);
		fis.betau[1][i] = c*h/Keff;
		Keff = (2*fis.perm[1][i]*fis.perm[1][i-1])/(fis.perm[1][i]+fis.perm[1][i-1]);
		fis.betad[1][i] = c*h/Keff;
		Keff = (2*fis.perm[1][i]*fis.perm[2][i])/(fis.perm[1][i]+fis.perm[2][i]);
		fis.betar[1][i] = c*h/Keff;
		 		
		Keff = (2*fis.perm[n][i]*fis.perm[n][i+1])/(fis.perm[n][i]+fis.perm[n][i+1]);
		fis.betau[n][i] = c*h/Keff;
		Keff = (2*fis.perm[n][i]*fis.perm[n][i-1])/(fis.perm[n][i]+fis.perm[n][i-1]);
		fis.betad[n][i] = c*h/Keff;
		Keff = (2*fis.perm[n][i]*fis.perm[n-1][i])/(fis.perm[n][i]+fis.perm[n-1][i]);
		fis.betal[n][i] = c*h/Keff;   		
	}
		 

	for(j=2; j<n; j++)
		for(k=2; k<n; k++)
		{     
			Keff = (2*fis.perm[j][k]*fis.perm[j][k+1])/(fis.perm[j][k]+fis.perm[j][k+1]);
		   fis.betau[j][k] = c*h/Keff;
		   Keff = (2*fis.perm[j][k]*fis.perm[j][k-1])/(fis.perm[j][k]+fis.perm[j][k-1]);
		   fis.betad[j][k] = c*h/Keff;
		   Keff = (2*fis.perm[j][k]*fis.perm[j+1][k])/(fis.perm[j][k]+fis.perm[j+1][k]);
		   fis.betar[j][k] = c*h/Keff;
		   Keff = (2*fis.perm[j][k]*fis.perm[j-1][k])/(fis.perm[j][k]+fis.perm[j-1][k]);
		   fis.betal[j][k] = c*h/Keff;
	    }     
		 
	/*
	 * Inicializando valores da fonte
	 */  
	fis.f[1][1]=1.0e-7;
	fis.f[n][n]=-1.0e-7;
		 
	/*
	 * calculo de parâmetros que nao dependem das iterações
	 */
	
	aux = 1/h;
	
	for (j=1; j<=n; j++)
		for (k=1; k<=n; k++)
		{
	  		fis.shi[j][k]=2*fis.perm[j][k]*aux;
			fis.f[j][k]*=h;
		}

	/*
    * Ciclo até convergência do problema
    */
  	i = 1; 
    
	for (; ;)
	{	

	  	/*Cálculo da pressão e dos Fluxosem cada elemento */

		/*Canto inferior esquerdo [1][1]*/
      canto_d_l(1, 1,f,m,fis);

		/*Canto superior esquerdo [1][N-2]*/
		canto_u_l(1, n,f,m,fis);
		    
		/*Canto inferior direito [N-2][1]*/
		canto_d_r(n, 1,f,m,fis);

		/*Canto superior direito [N-2][N-2]*/
		canto_u_r(n, n,f,m,fis);		    

		/*Fronteira U [2...N-3][N-2]*/
		for (j=2; j<n; j++)
			fronteira_u(j, n,f,m,fis);
		    	
		/*Fronteira D [2...N-3][1]*/
		for (j=2; j<n; j++)
			fronteira_d(j, 1,f,m,fis);

		/*Fronteira R [N-2][2...N-3]*/
	   for (k=2; k<n; k++)		    
			fronteira_r(n, k,f,m,fis);      
			   
		/*Fronteira L [1][2...N-3]*/
		for (k=2; k<n; k++)
			fronteira_l(1, k,f,m,fis);
		   	
		/*Elementos internos [2..N-3][2..N-3]*/
		for (j=2; j<n; j++)
			for (k=2; k<n; k++)
			{
		  		AuxL = fis.shi[j][k]/(1+fis.betal[j][k]*fis.shi[j][k]);
	  			AuxR = fis.shi[j][k]/(1+fis.betar[j][k]*fis.shi[j][k]);
	  			AuxU = fis.shi[j][k]/(1+fis.betau[j][k]*fis.shi[j][k]);
				AuxD = fis.shi[j][k]/(1+fis.betad[j][k]*fis.shi[j][k]);
				DL = AuxL*(fis.betal[j][k]*f.qr_old[j-1][k]+m.lr_old[j-1][k]);
				DR = AuxR*(fis.betar[j][k]*f.ql_old[j+1][k]+m.ll_old[j+1][k]);
				DU = AuxU*(fis.betau[j][k]*f.qd_old[j][k+1]+m.ld_old[j][k+1]);
				DD = AuxD*(fis.betad[j][k]*f.qu_old[j][k-1]+m.lu_old[j][k-1]);
	  	
				fis.p[j][k] = (fis.f[j][k] + DU + DD + DR + DL)/(AuxU + AuxR + AuxD + AuxL);
				f.ql[j][k] = AuxL*fis.p[j][k] - DL;
				f.qr[j][k] = AuxR*fis.p[j][k] - DR;
				f.qu[j][k] = AuxU*fis.p[j][k] - DU;
				f.qd[j][k] = AuxD*fis.p[j][k] - DD;			
			}     	
					    	
		
			/* 
			 * Atualização dos multiplicadores de lagrange 
			 */
			
			M = 0.0;

			for (j=1; j<=n; j++)
				for (k=1; k<=n; k++)
			  	{
					m.lu[j][k] = fis.betau[j][k]*(f.qu[j][k] + f.qd_old[j][k+1]) + m.ld_old[j][k+1];
					m.ld[j][k] = fis.betad[j][k]*(f.qd[j][k] + f.qu_old[j][k-1]) + m.lu_old[j][k-1];
					m.lr[j][k] = fis.betar[j][k]*(f.qr[j][k] + f.ql_old[j+1][k]) + m.ll_old[j+1][k];
					m.ll[j][k] = fis.betal[j][k]*(f.ql[j][k] + f.qr_old[j-1][k]) + m.lr_old[j-1][k];
					M += fis.p[j][k];
				}                     


			M = M / (n*n);	 

			sum1 = 0.;
			sum2 = 0.; 
  
			/* Impondo a média zero na distriubição de pressões
			 * Início de cálculo de verificação de convergência
			 */
			for (j=1; j<=n; j++)
		   	for (k=1; k<=n; k++)
		  		{
					fis.p[j][k] -= M;
			  		m.lu[j][k] -= M;
			  		m.ld[j][k] -= M;
			  		m.lr[j][k] -= M;
			  		m.ll[j][k] -= M;

					aux = fis.p[j][k] - fis.p_old[j][k];
			  		sum1 += aux*aux;
			  		sum2 += fis.p[j][k]*fis.p[j][k];
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


/* Função para o canto inferior esquerdo*/	
void canto_d_l(const int _j, const int _k,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{

	double AuxU,AuxR,DU,DR;	/* Auxiliares para cada lado das células */

	AuxU = fis.shi[_j][_k]/(1+fis.betau[_j][_k]*fis.shi[_j][_k]);
   AuxR = fis.shi[_j][_k]/(1+fis.betar[_j][_k]*fis.shi[_j][_k]);
   DU = AuxU*(fis.betau[_j][_k]*f.qd_old[_j][_k+1]+m.ld_old[_j][_k+1]);
   DR = AuxR*(fis.betar[_j][_k]*f.ql_old[_j+1][_k]+m.ll_old[_j+1][_k]);
   fis.p[_j][_k] = (fis.f[_j][_k] + DU + DR)/(AuxU + AuxR);
   f.qu[_j][_k] = AuxU*fis.p[_j][_k] - DU;
   f.qr[_j][_k] = AuxR*fis.p[_j][_k] - DR;
	   
}

/* Função para o canto superior esquerdo*/	
void canto_u_l(const int _j, const int _k,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{

	double AuxD, AuxR,DD, DR;	/* Auxiliares para cada lado das células */

	AuxD = fis.shi[_j][_k]/(1+fis.betad[_j][_k]*fis.shi[_j][_k]);
   AuxR = fis.shi[_j][_k]/(1+fis.betar[_j][_k]*fis.shi[_j][_k]);
   DD = AuxD*(fis.betad[_j][_k]*f.qu_old[_j][_k-1]+m.lu_old[_j][_k-1]);
   DR = AuxR*(fis.betar[_j][_k]*f.ql_old[_j+1][_k]+m.ll_old[_j+1][_k]);
   fis.p[_j][_k] = (fis.f[_j][_k] + DD + DR)/(AuxD + AuxR);
   f.qd[_j][_k] = AuxD * fis.p[_j][_k] - DD;
   f.qr[_j][_k] = AuxR * fis.p[_j][_k] - DR;		

}

/* Função para o canto inferior direito*/	
void canto_d_r(const int _j, const int _k,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{

	double AuxU, AuxL,DU,DL;	/* Auxiliares para cada lado das células */

   AuxU = fis.shi[_j][_k]/(1+fis.betau[_j][_k]*fis.shi[_j][_k]);
   AuxL = fis.shi[_j][_k]/(1+fis.betal[_j][_k]*fis.shi[_j][_k]);
   DU = AuxU*(fis.betau[_j][_k]*f.qd_old[_j][_k+1]+m.ld_old[_j][_k+1]);
   DL = AuxL*(fis.betal[_j][_k]*f.qr_old[_j-1][_k]+m.lr_old[_j-1][_k]);
   fis.p[_j][_k] = (fis.f[_j][_k] + DU + DL)/(AuxU + AuxL);
   f.qu[_j][_k] = AuxU*fis.p[_j][_k] - DU;
   f.ql[_j][_k] = AuxL*fis.p[_j][_k] - DL;		

}

/* Funcao para o canto superior dereito*/		
void canto_u_r(const int _j, const int _k,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{

	double AuxD,AuxL,DD,DL;	/* Auxiliares para cada lado das células */

	AuxD = fis.shi[_j][_k]/(1+fis.betad[_j][_k]*fis.shi[_j][_k]);
   AuxL = fis.shi[_j][_k]/(1+fis.betal[_j][_k]*fis.shi[_j][_k]);
   DD = AuxD*(fis.betad[_j][_k]*f.qu_old[_j][_k-1]+m.lu_old[_j][_k-1]);
   DL = AuxL*(fis.betal[_j][_k]*f.qr_old[_j-1][_k]+m.lr_old[_j-1][_k]);
   fis.p[_j][_k] = (fis.f[_j][_k] + DD + DL)/(AuxD + AuxL);
   f.qd[_j][_k] = AuxD*fis.p[_j][_k] - DD;
   f.ql[_j][_k] = AuxL*fis.p[_j][_k] - DL;

}

/* Função para a fronteira esquerda L */		
void fronteira_l(const int _j, const int _k,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{

	double AuxU, AuxD, AuxR,DU, DD, DR; /* Auxiliares para cada lado das células */

	AuxU = fis.shi[_j][_k]/(1+fis.betau[_j][_k]*fis.shi[_j][_k]);
  	AuxD = fis.shi[_j][_k]/(1+fis.betad[_j][_k]*fis.shi[_j][_k]);
  	AuxR = fis.shi[_j][_k]/(1+fis.betar[_j][_k]*fis.shi[_j][_k]);
  	DU = AuxU*(fis.betau[_j][_k]*f.qd_old[_j][_k+1]+m.ld_old[_j][_k+1]);
  	DD = AuxD*(fis.betad[_j][_k]*f.qu_old[_j][_k-1]+m.lu_old[_j][_k-1]);
  	DR = AuxR*(fis.betar[_j][_k]*f.ql_old[_j+1][_k]+m.ll_old[_j+1][_k]);
  	fis.p[_j][_k] = (fis.f[_j][_k] + DU + DR + DD)/(AuxU + AuxR + AuxD);
  	f.qu[_j][_k] = AuxU*fis.p[_j][_k] - DU;
  	f.qd[_j][_k] = AuxD*fis.p[_j][_k] - DD;
  	f.qr[_j][_k] = AuxR*fis.p[_j][_k] - DR;		

}

/* Função para a fronteira dereita R */	
void fronteira_r(const int _j, const int _k,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{

	double AuxU,AuxD, AuxL,DU, DD, DL;	/* Auxiliares para cada lado das celulas */

	AuxU = fis.shi[_j][_k]/(1+fis.betau[_j][_k]*fis.shi[_j][_k]);
   AuxD = fis.shi[_j][_k]/(1+fis.betad[_j][_k]*fis.shi[_j][_k]);
   AuxL = fis.shi[_j][_k]/(1+fis.betal[_j][_k]*fis.shi[_j][_k]);
   DU = AuxU*(fis.betau[_j][_k]*f.qd_old[_j][_k+1]+m.ld_old[_j][_k+1]);
   DD = AuxD*(fis.betad[_j][_k]*f.qu_old[_j][_k-1]+m.lu_old[_j][_k-1]);
   DL = AuxL*(fis.betal[_j][_k]*f.qr_old[_j-1][_k]+m.lr_old[_j-1][_k]);
   fis.p[_j][_k] = (fis.f[_j][_k] + DU + DL + DD)/(AuxU + AuxL + AuxD);
   f.qu[_j][_k] = AuxU*fis.p[_j][_k] - DU;
   f.qd[_j][_k] = AuxD*fis.p[_j][_k] - DD;
   f.ql[_j][_k] = AuxL*fis.p[_j][_k] - DL;		

}

/* Função para a fronteira superior U */	
void fronteira_u(const int _j, const int _k,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{

	double AuxD, AuxR, AuxL, DD, DR, DL;	/* Auxiliares para cada lado das celulas */

	AuxL = fis.shi[_j][_k]/(1+fis.betal[_j][_k]*fis.shi[_j][_k]);
  	AuxR = fis.shi[_j][_k]/(1+fis.betar[_j][_k]*fis.shi[_j][_k]);
  	AuxD = fis.shi[_j][_k]/(1+fis.betad[_j][_k]*fis.shi[_j][_k]);
  	DL = AuxL*(fis.betal[_j][_k]*f.qr_old[_j-1][_k]+m.lr_old[_j-1][_k]);
  	DR = AuxR*(fis.betar[_j][_k]*f.ql_old[_j+1][_k]+m.ll_old[_j+1][_k]);
  	DD = AuxD*(fis.betad[_j][_k]*f.qu_old[_j][_k-1]+m.lu_old[_j][_k-1]);
  	fis.p[_j][_k] = (fis.f[_j][_k] + DD + DL + DR)/(AuxD + AuxL + AuxR);
  	f.ql[_j][_k] = AuxL*fis.p[_j][_k] - DL;
  	f.qr[_j][_k] = AuxR*fis.p[_j][_k] - DR;
  	f.qd[_j][_k] = AuxD*fis.p[_j][_k] - DD;		

}
		
/* Função para a fronteira inferior D */	
void fronteira_d(const int _j, const int _k,  Fluxos f,  Multiplicadores m,  Fisicos fis)
{
	double AuxU,AuxR, AuxL,DU,DR, DL;	/* Auxiliares para cada lado das celulas */

  	AuxL = fis.shi[_j][_k]/(1+fis.betal[_j][_k]*fis.shi[_j][_k]);
  	AuxR = fis.shi[_j][_k]/(1+fis.betar[_j][_k]*fis.shi[_j][_k]);
  	AuxU = fis.shi[_j][_k]/(1+fis.betau[_j][_k]*fis.shi[_j][_k]);
  	DL = AuxL*(fis.betal[_j][_k]*f.qr_old[_j-1][_k]+m.lr_old[_j-1][_k]);
  	DR = AuxR*(fis.betar[_j][_k]*f.ql_old[_j+1][_k]+m.ll_old[_j+1][_k]);
  	DU = AuxU*(fis.betau[_j][_k]*f.qd_old[_j][_k+1]+m.ld_old[_j][_k+1]);
  	fis.p[_j][_k] = (fis.f[_j][_k] + DU + DL + DR)/(AuxU + AuxL + AuxR);
  	f.ql[_j][_k] = AuxL*fis.p[_j][_k] - DL;
  	f.qr[_j][_k] = AuxR*fis.p[_j][_k] - DR;
  	f.qu[_j][_k] = AuxU*fis.p[_j][_k] - DU;		

}

/* Função para os elementos interiores */	
void internos(const int _j, const int _k,  Fluxos f,  Multiplicadores m, Fisicos fis)
{
	double AuxU, AuxD, AuxR, AuxL,DU, DD, DR, DL;	/* Auxiliares para cada lado das celulas */

  	AuxL = fis.shi[_j][_k]/(1+fis.betal[_j][_k]*fis.shi[_j][_k]);
  	AuxR = fis.shi[_j][_k]/(1+fis.betar[_j][_k]*fis.shi[_j][_k]);
  	AuxU = fis.shi[_j][_k]/(1+fis.betau[_j][_k]*fis.shi[_j][_k]);
  	AuxD = fis.shi[_j][_k]/(1+fis.betad[_j][_k]*fis.shi[_j][_k]);
  	DL = AuxL*(fis.betal[_j][_k]*f.qr_old[_j-1][_k]+m.lr_old[_j-1][_k]);
  	DR = AuxR*(fis.betar[_j][_k]*f.ql_old[_j+1][_k]+m.ll_old[_j+1][_k]);
  	DU = AuxU*(fis.betau[_j][_k]*f.qd_old[_j][_k+1]+m.ld_old[_j][_k+1]);
  	DD = AuxD*(fis.betad[_j][_k]*f.qu_old[_j][_k-1]+m.lu_old[_j][_k-1]);
	  	
	fis.p[_j][_k] = (fis.f[_j][_k] + DU + DD + DR + DL)/(AuxU + AuxR + AuxD + AuxL);
  	f.ql[_j][_k] = AuxL*fis.p[_j][_k] - DL;
 	f.qr[_j][_k] = AuxR*fis.p[_j][_k] - DR;
  	f.qu[_j][_k] = AuxU*fis.p[_j][_k] - DU;
 	f.qd[_j][_k] = AuxD*fis.p[_j][_k] - DD;		

}

/*Função que aloca as matrizes de tamanho linXcol*/
double** cria_ponteiro(int lin,int col)
{

	int i;
	double** p_aux;

	
	p_aux = (double**) malloc(lin*sizeof(double*));

	if(p_aux==NULL)
	{
		printf("Erro na alocacao de memoria\n");
		return NULL;
	}

	for(i=0;i<lin;i++)
	{
		p_aux[i] = (double*) malloc(col*sizeof(double));

		if(p_aux[i]==NULL)
		{
			printf("Erro na alocacao de memoria\n");
			return NULL;
		}
		
	}

	return p_aux;

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

