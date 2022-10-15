## Modelo de Ising 2D

Modelo de Ising 2D paralelizado com MPI - HPC2

➡️ make clean; make; mpirun -np (númeroDeProcessos) ./ising.x

#### *ARQUIVOS*

> - isingfunc.c: definição das funções<br/>
> - define.h: definição dos parâmetros<br/>
> - spinsIni.dat e spinsFinal.dat: configuração inicial e final dos spins<br/>
> - latticeEvo.dat: evolução da configuração de spins a cada passo<br/>
> - propsFinal.dat: propriedades magnéticas do sistema para cada temperatura<br/>
> - time.dat: tempo de execução do programa
