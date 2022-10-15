/* Define parâmetros para a simulação */

#define N 100                 // Tamanho do lado da rede quadrada 

#define TOTALN (N*N)          // Número total de spins

#define INV_TOTALN 1./(N*N)   // Total inverso de spins -> 1/TOTALN = 1*INV_TOTALN

#define ERR_ENE 1e-5          // Erro de energia aceitável no cálculo do equilíbrio do sistema

#define B 0                   // Campo magnético externo

#define J 1                   // Energia de interação

#define TEMP_INI 48           // Temperatura inicial * 10 << para usar int no loop

#define TEMP_FINAL 1          // Temperatura final * 10 << para usar int no loop

#define TEMP_STEP 1           // Passo da temperatura * 10 << para usar int no loop

#define MAX_SIZE 48           // TempIni/TempStep

#define MC_STEPS 2100000000   // Limite de passos para encontrar equilíbrio
