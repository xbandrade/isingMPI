#include <stdio.h> 
#include <stdlib.h> 
#include <string.h>
#include <stdbool.h>
#include <time.h>
#include <math.h> 
#include <mpi.h>
#include "define.h"
#include "ising.h"

int main(int argc, char *argv[])
{
	FILE *fp0 = fopen("results/time.dat", "w");         // Salva o tempo gasto para rodar o programa
	FILE *fp1 = fopen("results/spinsIni.dat", "w");     // spins - configuração inicial
	FILE *fp2 = fopen("results/spinsFinal.dat", "w");   // spins - configuração final
	FILE *fp = fopen("results/propsFinal.dat", "w");    // salva propriedades calculadas
	FILE *fEvo = fopen("results/latticeEvo.dat", "w");  // salva a evolução da rede a cada temperatura

	int x, y, master = 0;
	int rank, size;
	int i, j;
	int **M;
	int Mini[N][N];
	M = alocaM(N, N); // aloca dinamicamente a matriz M
	double t1, t2;

    MPI_Status rstatus; // Status com informações da rotina
    MPI_Init(&argc, &argv); // Inicia ambiente MPI
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	
	t1 = MPI_Wtime(); // Grava o tempo inicial
	if(rank == master)
	{	
		randSpins(M);
		for(i=0;i<N;i++)
			for(j=0;j<N;j++)
				Mini[i][j] = M[i][j]; // Salva matriz inicial em uma matriz estática para bcast

		// grava o tempo inicial
		printf("\n- Parallel 2D Ising Model - Lattice size: %dx%d - Processes: %d -\n\n", N, N, size);

		for(x=0;x<N;x++) // Salva configuração inicial de spins no arquivo
		{ 
			for(y=0;y<N;y++)
				fprintf(fp1, M[x%N][y%N]==1 ? " 1 " : "-1 ");
			fprintf(fp1, "\n");
		}
		if(fp == NULL || fEvo == NULL)
		{ 
			fprintf(stderr, "%s", "File I/O Error\n");
			exit(1); 
		}
		printf("%s\t%s\t\t%s\t\t%s\t\t%s\t%s\n", "Temperature", "Energy", "Magnetization", "Heat Capacity", "Mag. Susceptibility", "Steps");
		fprintf(fp, "%s\t\t%s\t\t\t%s\t\t\t%s\t\t%s\t%s\n", "Temperature", "Energy", "Magnetization", "Heat Capacity", "Mag. Susceptibility", "Steps");	
	}
	// Envia a matriz iniciada para todos
    MPI_Bcast(&Mini, N*N, MPI_INT, master, MPI_COMM_WORLD);
	if(rank!=master) // Todos exceto o mestre salvam matriz inicial em uma matrix dinâmica
		for(i=0;i<N;i++)
			for(j=0;j<N;j++)
				M[i][j] = Mini[i][j];

	MPI_Barrier(MPI_COMM_WORLD); // Espera o processo mestre terminar suas ações antes de todos entrarem no loop
	tempLoop(M, fp, fEvo, size, rank); // Aplica o método para cada passo de temperatura a partir da matriz M
	MPI_Barrier(MPI_COMM_WORLD); // Espera todos os processos sairem do loop

	if(rank == master)
	{
		t2 = MPI_Wtime(); // Grava o tempo final
		printf("\nElapsed time: %f seconds\n", t2 - t1);
		fprintf(fp0, "\nElapsed time: %f seconds\n", t2 - t1);
		fclose(fp);	
		fclose(fp0);
		fclose(fp1);	
		fclose(fp2);	
		fclose(fEvo);	// Fecha os arquivos
		printf("\nDone!\n");
	}

    MPI_Finalize(); // Finaliza o ambiente MPI
	return 0;
} 
