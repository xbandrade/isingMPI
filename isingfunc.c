/* Definições das funções */
#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h> 
#include <mpi.h>
#include "ising.h"
#include "define.h"

unsigned int seed = 5;     // Seed base do gerador de números aleatórios

struct somaEM sEM = {0, 0, 0, 0, 0, 0};
struct prop props;


int **alocaM(int m, int n) // Aloca dinamicamente uma matriz de int m x n
{
   int i, **M;
   M = malloc(m*sizeof(int *));
   for (i=0;i<m;i++)
      M[i] = malloc(n*sizeof(int));
   return M;
}


void randSpins(int **M) // Inicializa M com spins aleatórios -1 ou 1
{
  int i, j;
  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      M[i][j] = newRandInt(0)%2 ? -1 : 1;
}

int somaAdj(int **M, int x, int y) // Soma os spins adjacentes
{
	// Condições de contorno periódico na rede -> %N conecta o fim ao início
	return M[(x+1)%N][y%N] + M[(x-1)%N][y%N] + M[x%N][(y+1)%N] + M[x%N][(y-1)%N];
}

double deltaE(int **M, int x, int y, int sAdj)  // Calcula ΔE
{
	// Condições de contorno periódico na rede -> %N conecta o fim ao início
	return 2*J*M[x%N][y%N]*(sAdj);  
}

double eneTotal(int **M)  // Calcula energia média total
{
  int x, y;
  double en = 0;
  for(x=0;x<N;x++)
    for(y=0;y<N;y++)
      en -= M[x%N][y%N]*(M[(x+1)%N][y%N] + M[x%N][(y+1)%N]);
  return en*INV_TOTALN; 
}

double magTotal(int **M) // Calcula magnetização média total
{
  int x, y;
  double m = 0;
	for(x=0;x<N;x++)
		for(y=0;y<N;y++)
			m += M[x][y];
	return m*INV_TOTALN;
}

bool equilibrioM(double totalEne, double oldEne) // Verifica equilíbrio do sistema
{
	// Compara a energia total do sistema atual com a energia total anterior
	// Se o resultado for aceitável por mais de uma vez, o sistema está em equilíbrio
	if(fabs((oldEne-totalEne)/totalEne)<ERR_ENE) return true;
	else return false; 
}

double capTermica(struct somaEM *sEM, double temp) // Calcula capacidade térmica
{
	double E2M = sEM->somaE2/(double)sEM->cont; //<E²>
	double EM2 = ((sEM->somaE)*(sEM->somaE))/(((double)sEM->cont)*((double)sEM->cont)); //<E>²
	return (E2M-EM2)/(temp*temp);  // Cv = (<E²>-<E>²)/T²
}

double susceptMagnetica(struct somaEM *sEM, double temp) // Calcula susceptibilidade magnética
{ 
	double M2M = sEM->somaM2/(double)sEM->cont;                                 // <M²>
	double MM2 = (sEM->somaM/(double)sEM->cont)*(sEM->somaM/(double)sEM->cont); // <M>²
	return (M2M-MM2)/temp;                                                      // χ = (<M²>-<M>²)/T
}

double *calcFlipProb(double temp) // calcula previamente todas as exponenciais possíveis
{
   int k;
   double *p;
   p = malloc(9*sizeof(double));    // valores de spins adjacentes possíveis são 9: [-4, 4]
   for(k=0;k<9;k++)                 // spins -4 -> k=0, spins -3 -> k=1, spins -2 -> k=2, etc
      p[k] = exp((-2*J*(k-4)+B)/temp);
   return p;
}

void propsM(struct somaEM *sEM, struct prop *props, int **M, double temp) // Salva propriedades para cada temperatura
{
	props->energia            = eneTotal(M);
	props->magnetiz	          = magTotal(M); 
	props->capacidadeTermica  = capTermica(sEM, temp); 
	props->suscept            = susceptMagnetica(sEM, temp); 
	props->tempEq             = sEM->tempEq;
}

int newRandInt(int rank) // gera um número aleatorio inteiro [0, 32767]
{
	unsigned int k; // n>=0
	seed = (seed + 3*rank) * 1103515245 + 12345; // muda para cada rank
	k=((unsigned)(seed/65536)%32768);
	//return (double)k/32768;
   return k;
}

double newRandDouble(int rank)  // gera um número aleatorio double [0, 1]
{
	unsigned int k; // n>=0
	seed = (seed + 3*rank) * 1103515245 + 12345; // muda para cada rank
	k=((unsigned)(seed/65536)%32768);
	return (double)k/32768;
}

void monteCarlo(int **M, struct somaEM *sEM, double temp, int size, int rank) // Usa o Método de Monte Carlo até o sistema ficar equilibrado
{ 	
	int x = 0, y = 0, t = 0, k; 
	int sAdj;	                           // soma dos spins vizinhos
	double dE, dM;                         // delta E, delta M
	double ene = eneTotal(M);	           // energia do sistema inicial
	double mag = magTotal(M);              // magnetização do sistema inicial
	double oldEne = ene;	               // salva energia anterior para verificar equilíbrio
  	double *p;                             // p = exp((-2*J*(int)sAdj+B)/temp)
  	double probFlip;                       // probabilidade de flipar o spin
	bool eqM = false, doubleCheck = false; // Verifica equilibrio e faz um double check para confirmar

	p = calcFlipProb(temp);                // salva no array p os valores das exponenciais calculadas

	do{ 
		t++;  // contador
		x     = 1+newRandInt(rank)%N;
		y     = 1+newRandInt(rank)%N;
      	sAdj  = somaAdj(M, x, y);
		dE    = deltaE(M, x, y, sAdj)*INV_TOTALN;
		dM    = -2*M[x%N][y%N]*INV_TOTALN;
      
		//probFlip = M[x%N][y%N] == -1 ? 1/p[sAdj+4] : p[sAdj+4]; // checa se o spin é up ou down

		probFlip = p[sAdj+4];

		if(dE < 0 || newRandDouble(rank) < probFlip) // Verifica se o spin vai ser flipado
		{
			M[x%N][y%N] *= -1; 
			ene         += dE;					
			mag         += dM; 	
		}

		sEM->somaE2 += (ene*ene); // soma para calculo de <E>² e <E²>
		sEM->somaE	+= ene;
		sEM->somaM2 += (mag*mag); // soma para calculo de <M>² e <M²>
		sEM->somaM	+= mag;
		sEM->cont++;              // contador

		if(!(t%(N*N)) && !eqM) // faz verificação de equilíbrio a cada (N*N) passos de Monte Carlo
		{
			eqM = equilibrioM(ene, oldEne); 			 
			oldEne = ene;
			if(eqM && !doubleCheck) // força um double check pra confirmar equilibrio duas vezes seguidas
			{	
				doubleCheck = true; 			
				eqM = false;
			}
			else if(!eqM && doubleCheck)    // doubleCheck falhou -> não houve equilibrio 2x seguidas
				doubleCheck = false;
			else if(eqM && doubleCheck)	 // doubleCheck funcionou -> aceita o equilibrio do sistema
			{			
            sEM->tempEq = t;
            break; 		
			}
		}
	}while(t<MC_STEPS);
} 

void tempLoop(int **M, FILE *fp, FILE *fEvo, int size, int rank)
{ 
	int x, y, i, j;
	int start, end, step;
	int master = 0;
	int temp; 
	double t;

	step = size*TEMP_FINAL;
	end = (rank+1) * TEMP_FINAL;
	start = end + ((TEMP_INI/size)-TEMP_FINAL)*size; //

	for(temp=start;temp>=end;temp-=step)
	{ 
		t = 1.*temp/10; // Volta a temperatura para double
		monteCarlo(M, &sEM, t, size, rank); // Aplica o método de Monte Carlo para as temperaturas
		propsM(&sEM, &props, M, t);         // Atribui valores das propriedades na struct
		props.temperatura = t;	

		fprintf(fEvo, "%f\n", t);
		for(x=0;x<N;x++) // Salva configuração de spins da temperatura atual no arquivo
		{ 
			for(y=0;y<N;y++)
			fprintf(fEvo, M[x%N][y%N]==1 ? " 1 " : "-1 ");
			fprintf(fEvo, "\n");
		}
		fprintf(fEvo, "\n\n");

		//Salva as propriedades calculadas no arquivo
		if(props.tempEq==0) // não encontrou equilíbrio no limite de passos
		{
			fprintf(fp, "%2.1f\t\t%16.12f\t\t%11.12f\t\t%15.12f\t\t%17.12f\t\t%d\n", props.temperatura, props.energia, props.magnetiz,    
				props.suscept, props.capacidadeTermica, MC_STEPS);
			printf("%2.1f\t%16.8f\t%11.8f\t%15.8f\t\t%17.12f\t%14d\n", props.temperatura, props.energia, props.magnetiz, props.suscept,  
				props.capacidadeTermica, MC_STEPS);
		}
		else // encontrou equilíbrio dentro do limite de passos
		{
			fprintf(fp, "%2.1f\t\t%16.12f\t\t%11.12f\t\t%15.12f\t\t%17.12f\t\t%d\n", props.temperatura, props.energia, props.magnetiz,    
				props.suscept, props.capacidadeTermica, props.tempEq);
			printf("%2.1f\t%16.8f\t%11.8f\t%15.8f\t\t%17.12f\t%14d\n", props.temperatura, props.energia, props.magnetiz, props.suscept,  
				props.capacidadeTermica, props.tempEq);
		}
		sEM = (struct somaEM){0, 0, 0, 0, 0, 0}; // Reseta a struct e passa para a próxima iteração
	}
}

