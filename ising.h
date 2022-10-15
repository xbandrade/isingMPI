/* Header file para funções e estruturas */

/* STRUCTS */
struct somaEM  // Guarda somas de <E²>, <E>², <M²> e <M>²
{
	double somaE, somaE2;	  //soma de E e E²
	double somaM, somaM2;	  //soma de M e M²
	int cont;	    		 
	int tempEq;				  
};

struct prop // Guarda as propriedades calculadas
{ 		
	double temperatura, capacidadeTermica, energia, magnetiz, suscept; // 
	int tempEq;
}; 


/* FUNÇÕES */
int newRandInt(int rank); 								  // Retorna int [0, 32767]


double newRandDouble(int rank); 						  // Retorna double [0, 1]


void monteCarlo(int **M, struct somaEM *sEM, double temp, int size, int rank); // Usa o método de Monte Carlo até encontrar equilíbrio na rede


void tempLoop(int **M, FILE *fp, FILE *fEvo, int size, int rank); // Aplica MC para todas as temperaturas 


int **alocaM(int m, int n); 							  // Aloca matrix m x n dinamicamente


void randSpins(int **M); 								  // Inicializa M com spins -1 ou 1


int somaAdj(int **M, int x, int y); 					  // Soma dos spins adjacentes


double deltaE(int **M, int x, int y, int sAdj);  		  // Calcula ΔE


double eneTotal(int **M); 								  // Calcula energia total média


bool equilibrioM(double totalEne, double oldEne); 		  // Verifica equilíbrio do sistema


double capTermica(struct somaEM *sEM, double temp); 	  // Calcula capacidade térmica


double susceptMagnetica(struct somaEM *sEM, double temp); // Calcula susceptibilidade magnética


double *calcFlipProb(double temp); 						  // Calcula todas as exp. da distribuição de Boltzmann


void propsM(struct somaEM *sEM, struct prop *props, int **M, double temp); // Guarda propriedades EM para cada temperatura