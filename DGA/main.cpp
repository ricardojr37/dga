#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>

#define MAXd 250 // valor máximo de d  - número de doadores
#define MAXr 250 // valor máximo de r  - número de receptores
#define MAXh 200   // valor máximo do número de hospitais
#define NPop 30  // valor máximo do número de individuos da população


float sol[MAXr], sopt[MAXr], custo_opt;
float  mut_opt, pop_opt,  p[MAXr][MAXr], pc[MAXr][MAXr], set[MAXr];
float cross_opt;

char *recep[120] ={"ctsR60_16.txt", "ctsR60_1_60D.txt", "ctsR100_1_100D.txt"};
char *doad[120] ={"ctsD32_16.txt", "ctsD60_1_60R.txt", "ctsD100_1_100R.txt"};
char *hosps[120] ={"H_CH27_7.txt", "H_CH.txt", "H_CH (3).txt"};
int id_doador = -1, id_hospital =-1;

struct sreceptor {
  int tipo_sanguineo;
  float localizacao_x;
  float localizacao_y;
  int prioridade;
};
struct sdoador {
  int tipo_sanguineo;
  float localizacao_x;
  float localizacao_y;
};
struct shospital{
  int custo;
  float localizacao_x;
  float localizacao_y;
};
int numReceptor=0, numDoador = 0, numHospital=0, melhorH=0;
sreceptor receptor[MAXr];
sdoador doador[MAXd];
shospital hospital[MAXh];
long double dist=0, somatorio[NPop], pesoW[MAXr];



bool Existe(int valores[],int tam, int valor){
    for(int i = 0;i<tam;i++){
        if(valores[i]==valor)
            return true;
    }
    return false;
}
int Repetido(int vetor[],int tam){//Retorna a posicao do primeiro valor repetido no vetor
  int numeros_existentes[tam];
  numeros_existentes[0]=vetor[0];
  int k=1;
  for(int i=1;i<tam;i++){
    for(int j=0;j<k;j++){
      if(vetor[i]==numeros_existentes[j]){
        return i;
      }
    }
    numeros_existentes[k]=vetor[i];
    k++;
  }
  return -1;
}

void GeraAleatorios(int numeros[],int quantNumeros,int Limite){
    int v;
    int i = 1;
    while(i<=quantNumeros){
      v = rand() % Limite;
      if(v==0){
        while(v==0){
          v = rand() % Limite;
        }
      }
      while(Existe(numeros,i,v)){
        v = rand() % Limite;
        if(v==0){
        while(v==0){
          v = rand() % Limite;
        }
        }
      }
      numeros[i] = v;
      //printf("%d ", numeros[i]);
      i++;
    }
    //printf("\n") ;
}



/* ======================================================================
				   timeused
   ====================================================================== */
/* This timing routine is based on the ANSI-C procedure "clock", which
 * has a resolution of 1000000 ticks per second. This however implies
 * that we pass the limit of a long integer after only 4295 seconds.
 * The following routine attempts to correct such situations by adding
 * the constant ULONG_MAX to the counter whenever wraparound can be
 * detected. But the user is advised to use a timing routine like "times"
 * (which however is not ANSI-C standard) for measuring longer time
 * periods.                                                               */

void timeused(double *time)
{
  static double tstart, tend, tprev;

  if (time == NULL) {
    clock(); /* one extra call to initialize clock */
    tstart = tprev = clock();
  } else {
    tend = clock();
    if (tend < tprev) tstart -= ULONG_MAX; /* wraparound occured */
    tprev = tend;
    *time = (tend-tstart) / CLOCKS_PER_SEC; /* convert to seconds */
  }
}
/* ========================================================================= */
/* ================   Calcular distancia   ================================= */
double distancia(float x1, float x2, float y1, float y2){
  double distance=sqrt((pow(x1-x2,2)+pow(y1-y2,2)));
  return distance;
}


/* ========================================================================= */
/* ================   Entrada de dados  ==================================== */

void input(int h){
   int i, j, k, l;
  //Pegar os valores dos receptores
  FILE *file;
  file = fopen(recep[h],"r");
  fscanf(file, "%d ", &numReceptor);
  //printf("Problema com n=%d receptores\n", numReceptor);
  for (j = 1; j <= numReceptor; j++) {
    fscanf(file, "%f", &receptor[j].localizacao_x);
    fscanf(file, "%f ", &receptor[j].localizacao_y);
    fscanf(file, "%d ", &receptor[j].tipo_sanguineo);
    fscanf(file, "%d ", &receptor[j].prioridade);
  }
  fclose(file);
  //Pegar os valores dos doadores
  file = fopen(doad[h],"r");
  fscanf(file, "%d ", &numDoador);
  //printf("\nProblema com n=%d doadores\n", numDoador);
  for (j = 0; j < numDoador; j++) {
    fscanf(file, "%f", &doador[j].localizacao_x);
    fscanf(file, "%f ", &doador[j].localizacao_y);
    fscanf(file, "%d ", &doador[j].tipo_sanguineo);
  }
  fclose(file);
  //Pegar os valores dos hospitais
  file = fopen(hosps[h],"r");
  fscanf(file, "%d ", &numHospital);
  //printf("\nProblema com n=%d hospitais\n", numHospital);
  for (j = 0; j < numHospital; j++) {
    fscanf(file, "%f", &hospital[j].localizacao_x);
    fscanf(file, "%f ", &hospital[j].localizacao_y);
    fscanf(file, "%d ", &hospital[j].custo);
  }
  fclose(file);

  //Calcular o PESO W
 double maxRH = 0, maxDH = 0, maxCH = 0, soma = 0, maior = 0, aux_copia=0, maior_indice=0, aux2, resultado = 0;
  for(i=0;i<MAXr;i++)pesoW[i]=0;
  //printf("%lf", maior);
  for(int hosp=0;hosp<numHospital;hosp++){
    for(j = 1;j<=numReceptor;j++){
      if(hosp==0)pesoW[receptor[j].prioridade] +=1;
      if(receptor[j].prioridade>maior_indice)maior_indice=receptor[j].prioridade;
      for(i = 0;i<numDoador;i++){
        soma =distancia(doador[i].localizacao_x, hospital[hosp].localizacao_x, doador[i].localizacao_y, hospital[hosp].localizacao_y)+distancia(receptor[j].localizacao_x, hospital[hosp].localizacao_x, receptor[j].localizacao_y, hospital[hosp].localizacao_y)+hospital[hosp].custo;
        if(maior == 0.0) maior = soma;
        if(soma > maior) maior = soma;
      }
    }
  }
  for(i = maior_indice;i>0;i--){
    if(i==maior_indice){
      aux_copia=pesoW[i];
      pesoW[i] = maior;
    }
    else{
      aux2=pesoW[i];
      pesoW[i]=(aux_copia+1)*pesoW[i+1];
      aux_copia=aux2;
    }
  }
}

int copia_idD, copia_idH;

/* ========================================================================= */
/* ================   Funcao Custo g    ==================================== */
void g(float s2[MAXr]){
  int i, j;
  float tg;
  double soma=0, menor = 0;
  long double resultado = 0;
  int aux_doador[MAXd];
  id_doador = -1;
  id_hospital =-1;
  for(i = 0;i<numDoador;i++) aux_doador[i] = doador[i].tipo_sanguineo;
  for(int r=1;r<=numReceptor;r++){
    //printf("%d ", (int)s2[r]);
    id_doador = -1;
    id_hospital =-1;
    copia_idH = id_hospital;
    copia_idD = id_doador;
    menor = 0;
    //soma=0;
    for(int i = 0;i<numDoador;i++){
    //printf("%d ", aux_doador[i]);
    if(receptor[(int)s2[r]].tipo_sanguineo == 1){
      if(aux_doador[i] == 1){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          //printf("%lf-", soma);
        }
      }
    }
    if(receptor[(int)s2[r]].tipo_sanguineo == 2){
      if(aux_doador[i] == 2){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)s2[r]].tipo_sanguineo == 3){
      if(aux_doador[i] == 3){
        soma = 0;
       for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)s2[r]].tipo_sanguineo == 4){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }

        }
      }
    }
    if(receptor[(int)s2[r]].tipo_sanguineo == 1){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)s2[r]].tipo_sanguineo == 2){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)s2[r]].tipo_sanguineo == 3){
       if(aux_doador[i] == 4 || aux_doador[i] == 1|| aux_doador[i] == 2){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;

          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    //printf("%d\n", menor);
    }
    //printf("-> %lf \n", soma);
    if(menor!=0){
    //copia_idH = id_hospital;
    //copia_idD = id_doador;
    aux_doador[id_doador] = 0;
    resultado += pesoW[receptor[(int)s2[r]].prioridade] - menor;
    //return menor;
    }
  }
  tg=resultado;
  if (tg>custo_opt) {
    custo_opt=tg; sopt[0]=tg;
    for (j = 1; j <= numReceptor; j++) sopt[j]=s2[j];
   }
}
/* ========================================================================= */
/* ================   Função para saber a eficiencia =============== */
int eficiencia(int s2[MAXr]){
  int i, j;
  float tg;
  int contador=0;
  double soma=0, menor = 0;
  long double resultado = 0;
  int aux_doador[MAXd];
  id_doador = -1;
  id_hospital =-1;
  for(i = 0;i<numDoador;i++) aux_doador[i] = doador[i].tipo_sanguineo;
  for(int r=1;r<=numReceptor;r++){
    //printf("%d ", (int)s2[r]);
    id_doador = -1;
    id_hospital =-1;
    copia_idH = id_hospital;
    copia_idD = id_doador;
    menor = 0;
    //soma=0;
    for(int i = 0;i<numDoador;i++){
      //printf("%d ", aux_doador[i]);
      if(receptor[(int)s2[r]].tipo_sanguineo == 1){
      if(aux_doador[i] == 1){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          //printf("%lf-", soma);
        }
      }
      }
      if(receptor[(int)s2[r]].tipo_sanguineo == 2){
      if(aux_doador[i] == 2){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
      }
      if(receptor[(int)s2[r]].tipo_sanguineo == 3){
      if(aux_doador[i] == 3){
        soma = 0;
       for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
      }
      if(receptor[(int)s2[r]].tipo_sanguineo == 4){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }

        }
      }
      }
      if(receptor[(int)s2[r]].tipo_sanguineo == 1){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
      }
      if(receptor[(int)s2[r]].tipo_sanguineo == 2){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
      }
      if(receptor[(int)s2[r]].tipo_sanguineo == 3){
       if(aux_doador[i] == 4 || aux_doador[i] == 1|| aux_doador[i] == 2){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;

          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
      }
      //printf("%d\n", menor);
      }
      //printf("-> %lf \n", soma);
      if(menor!=0){
      //copia_idH = id_hospital;
      //copia_idD = id_doador;
      aux_doador[id_doador] = 0;
      contador++;
      //return menor;
      }
  }
  return contador;
}
/* ========================================================================= */
/* ================   Funcao Custo g para população inicial  =============== */

void g_p(float s2[MAXr]){
  int i, j, k, w, nigual;
  float tg;
  double soma=0, menor = 0;
  long double resultado = 0;
  int aux_doador[MAXd];
  for(int i = 0;i<numDoador;i++) aux_doador[i] = doador[i].tipo_sanguineo;
  for(int r=1;r<=numReceptor;r++){
    //printf("%d ", (int)s2[r]);
    id_doador = -1;
    id_hospital =-1;
    copia_idH = id_hospital;
    copia_idD = id_doador;
    menor = 0;
    //soma=0;
    for(int i = 0;i<numDoador;i++){
    //printf("%d ", aux_doador[i]);
    if(receptor[(int)s2[r]].tipo_sanguineo == 1){
      if(aux_doador[i] == 1){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          //printf("%lf-", soma);
        }
      }
    }
    if(receptor[(int)s2[r]].tipo_sanguineo == 2){
      if(aux_doador[i] == 2){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)s2[r]].tipo_sanguineo == 3){
      if(aux_doador[i] == 3){
        soma = 0;
       for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)s2[r]].tipo_sanguineo == 4){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }

        }
      }
    }
    if(receptor[(int)s2[r]].tipo_sanguineo == 1){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)s2[r]].tipo_sanguineo == 2){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)s2[r]].tipo_sanguineo == 3){
       if(aux_doador[i] == 4 || aux_doador[i] == 1|| aux_doador[i] == 2){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)s2[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)s2[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;

          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    //printf("%d\n", menor);
    }
    //printf("-> %lf \n", soma);
    if(menor!=0){
      //copia_idH = id_hospital;
      //copia_idD = id_doador;
      aux_doador[id_doador] = 0;
      resultado += pesoW[receptor[(int)s2[r]].prioridade] - menor;
      //return menor;
    }
  }
  tg=resultado;
  //printf("%Ld\n", tg);
  if (tg>custo_opt) {
    custo_opt=tg; sopt[0]=tg;
    for (j = 1; j <= numReceptor; j++) sopt[j]=s2[j];
  }
  //printf("\n");
  //inserir a solução s2 em P, se ela é melhor que a pior solução de P
  i=0; k=-1;
  while (i<NPop) {
    if (tg>p[i][0]) {
      w=i-1; k=i; i=NPop;
      while (w!=-1) {
        if (tg==p[w][0]) {
          nigual=0;
          for (j=1; j<=numReceptor; j++) if (p[w][j]!=s2[j]) nigual=nigual+1 ;
          if (nigual<(numReceptor/2)) {
            w=0; k=-1;
          }
          w=w-1;
        }
        else {
          w=-1;
        }
      } // end while
    } // end if custo<p[i][0] ...
       i=i+1;
  } // end while
  // insere ss em P na posição k
  if (k!=-1) {
    // Baixa todas as linhas da Matriz P a partir da linha k até a linha NPop-1
    for (j=NPop-2; j>=k; j--) {
      w=j+1;
      for (i=0; i<=numReceptor; i++) p[w][i]= p[j][i];
    }
    // inseri o vetor ss na linha k da matriz P
    for (i = 1; i <= numReceptor; i++) p[k][i]=s2[i]; p[k][0]=tg;
  } // fim do if (k!=-1)

}


/* ========================================================================= */
/* ================   Funcao Custo g para população do crosover  =========== */

void g_pc(long int ss[MAXr]){
  int i, j, k, w, nigual;
  float tg;
  double soma=0, menor = 0;
  long double resultado = 0;
  int aux_doador[MAXd];
  for(int i = 0;i<numDoador;i++) aux_doador[i] = doador[i].tipo_sanguineo;
  for(int r=1;r<=numReceptor;r++){
    //printf("%d ", (int)ss[r]);
    id_doador = -1;
    id_hospital =-1;
    copia_idH = id_hospital;
    copia_idD = id_doador;
    menor = 0;
    //soma=0;
    for(int i = 0;i<numDoador;i++){
    //printf("%d ", aux_doador[i]);
    if(receptor[(int)ss[r]].tipo_sanguineo == 1){
      if(aux_doador[i] == 1){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          //printf("%lf-", soma);
        }
      }
    }
    if(receptor[(int)ss[r]].tipo_sanguineo == 2){
      if(aux_doador[i] == 2){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)ss[r]].tipo_sanguineo == 3){
      if(aux_doador[i] == 3){
        soma = 0;
       for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)ss[r]].tipo_sanguineo == 4){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }

        }
      }
    }
    if(receptor[(int)ss[r]].tipo_sanguineo == 1){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)ss[r]].tipo_sanguineo == 2){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)ss[r]].tipo_sanguineo == 3){
       if(aux_doador[i] == 4 || aux_doador[i] == 1|| aux_doador[i] == 2){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;

          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    //printf("%d\n", menor);
    }
    //printf("-> %lf \n", soma);
    if(menor!=0){
    //copia_idH = id_hospital;
    //copia_idD = id_doador;
    aux_doador[id_doador] = 0;
    resultado += pesoW[receptor[(int)ss[r]].prioridade] - menor;
    //return menor;
  }
  //printf("%Lf\n", resultado);
  }
  tg=resultado;
   //printf("%ld\n", tg);
   //printf("%d\n", custo_opt);
  if (tg>custo_opt) {
    custo_opt=tg; sopt[0]=tg;
    for (j = 1; j <= numReceptor; j++) sopt[j]=ss[j];
  }

  //inserir a solução ss em P, se ela é melhor que a pior solução de P
  i=0; k=-1;
  while (i<NPop) {
    if (tg>pc[i][0]) {
      w=i-1; k=i; i=NPop;
      while (w!=-1) {
        if (tg==pc[w][0]) {
          nigual=0;
          for (j=1; j<=numReceptor; j++) if (pc[w][j]!=ss[j]) nigual=nigual+1;
          if (nigual<(numReceptor/2)) { w=0; k=-1; }
          w=w-1;
        }
        else {
          w=-1;
        }
      } // end while
    } // end if tg<p[i][0] ...
    i=i+1;
  } // end while
  // insere ss em P na posição k
  if (k!=-1) {
    // Baixa todas as linhas da Matriz P a partir da linha k até a linha NPop-1
    for (j=NPop-2; j>=k; j--) {
      w=j+1;
      for (i=0; i<=numReceptor; i++) pc[w][i]= pc[j][i];
    }
    // inseri o vetor ss na linha k da matriz P
    for (i = 1; i <= numReceptor; i++) pc[k][i]=ss[i]; pc[k][0]=tg;

  } // fim do if (k!=-1)
}

/* ====================================================================== */
/* ===============      Procedure População Inicial ===================== */

void PopulacaoInicial(){
    int i, j, k, v, ii, start, set[MAXr], s[MAXr], min;

    custo_opt=0;
    for (i=0; i<NPop; i++) p[i][0]=0;

    for(j = 0;j<NPop;j++){
    int numeros[numReceptor];
    GeraAleatorios(numeros,numReceptor,numReceptor+1);
      //printf("%d- ", j);
      for(i = 1;i<=numReceptor;i++){
        sol[i] = numeros[i];
        //printf("%d ", numeros[i]);
        //if(sol[i]==0) printf("AAAAAAAAAAAAAAAAA");
      }
      //printf("\n");
      g_p(sol);
    }
    //mostrar a popinicial
    /*
      for(j = 0;j<NPop;j++){
        printf("%d-> ", j);
        for(i = 0;i<=numReceptor;i++){
          printf("%d ", p[j][i]);
        }
        printf("\n");
      }
    */
}

/* ====================================================================== */
/* =====   Croosover com 1 Ponto de corte  ============================== */

void Cruzamento1P(){
   int i, j, k, v1, v2, c1;
   long int o1[MAXr], o2[MAXr], r1[MAXr], r2[MAXr];


   c1= numReceptor/2;
   //printf("  c1=%d \n", c1);

   //cruzamento de todos com todos
    for (i=0; i<=(NPop-2); i++) {
     for (j=i+1; j<=(NPop-1); j++) {

       // Coloca a 1a parte de cada cromossomo em cada um dos filhos
       for (k=1; k<=numReceptor; k++) { r1[k]=0; r2[k]=0; }
       for (k=1; k<=c1; k++) {
           o1[k]=p[i][k]; o2[k]=p[j][k];
           r1[ (int)p[i][k] ]= 1; r2[ (int)p[j][k] ]= 1;
       }
       // Preenche o restante do cromossomo de cada filho
       k=1; v1=c1+1; v2=c1+1;
       while (k<=numReceptor) {
         if (r1[ (int)p[j][k] ]==0) {
            o1[v1]=p[j][k];
            v1=v1+1;
         }
         if (r2[ (int)p[i][k] ]==0) {
            o2[v2]=p[i][k];
            v2=v2+1;
         }
         //if (i==0) printf("i=%d j=%d v1=%d v2=%d k=%d \n", i, j, v1, v2, k);
         k=k+1;
       }  // end while
       /*
       if ((i==0) && (j==1)) {
         printf("\n Um tipo de cruzamento \n");
         printf("P1 -");
         for (k=1; k<=numReceptor; k++) printf("%lu ", p[i][k]);
         printf("\n");
         printf("P2 -");
         for (k=1; k<=numReceptor; k++) printf("%lu ", p[j][k]);
         printf("\n");
         printf("O1 -");
         for (k=1; k<=numReceptor; k++) printf("%lu ", o1[k]);
         printf("\n");
         printf("O2 -");
         for (k=1; k<=numReceptor; k++) printf("%lu ", o2[k]);
         printf("\n");
       }
       */

       g_pc(o1); g_pc(o2);


       // Coloca a 2a parte de cada cromossomo em cada um dos filhos
       for (k=1; k<=numReceptor; k++) { r1[k]=0; r2[k]=0; }
       for (k=c1+1; k<=numReceptor; k++) {
           o1[k]=p[i][k]; o2[k]=p[j][k];
           r1[ (int)p[i][k] ]= 1; r2[ (int)p[j][k] ]= 1;
       }
       // Preenche o restante do cromossomo de cada filho
       k=1; v1=1; v2=1;
       while (k<=numReceptor) {
         if (r1[ (int)p[j][k] ]==0) {
            o1[v1]=p[j][k];
            v1=v1+1;
         }
         if (r2[ (int)p[i][k] ]==0) {
            o2[v2]=p[i][k];
            v2=v2+1;
         }
         //if (i==0) printf("i=%d j=%d v1=%d v2=%d k=%d \n", i, j, v1, v2, k);
         k=k+1;
       }  // end while
       g_pc(o1); g_pc(o2);

       //printf("i=%d j=%d custoOpt=%d \n", i, j, custo_opt);
     }  // end for j
   }   // end for i
}


/* ====================================================================== */
/* =====    Crossover com 2 pontos de corte    ========================== */

void Cruzamento2P(){
   int i, j, k, v1, v2, c1, c2;
   long int o1[MAXr], o2[MAXr], r1[MAXr], r2[MAXr];

   c1= numReceptor/3; c2= 2*c1; c1=c1 + 1; c2= c2 + 1;
   //printf("c1=%d c2=%d \n", c1, c2);

   //cruzamento de todos com todos
   for (i=0; i<=(NPop-2); i++) {
     for (j=i+1; j<=(NPop-1); j++) {
       // Coloca a 2a parte de cada cromossomo em cada um dos filhos
       for (k=1; k<=numReceptor; k++) { r1[k]=0; r2[k]=0; }
       for (k=c1; k<=c2; k++) {
           o1[k]=p[i][k]; o2[k]=p[j][k];
           r1[ (int)p[i][k] ]= 1; r2[ (int)p[j][k] ]= 1;
           //printf("%d ",p[j][k]);
       }
       // Preenche o restante do cromossomo de cada filho
       k=1; v1=1; v2=1;
       while (k<=numReceptor) {
         if (r1[ (int)p[j][k] ]==0) {
            o1[v1]=p[j][k];
            v1=v1+1;
            if (v1==c1) v1=c2+1;
         }
         if (r2[ (int)p[i][k] ]==0) {
            o2[v2]=p[i][k];
            v2=v2+1;
            if (v2==c1) v2=c2+1;
         }
         //if (i==0) printf("i=%d j=%d v1=%d v2=%d k=%d \n", i, j, v1, v2, k);
         k=k+1;
       }  // end while
       /*
        if ((i==0) && (j==1)) {
         printf("\n Um tipo de cruzamento \n");
         printf("P1 ->");
         for (k=1; k<=numReceptor; k++) printf("%ld ", p[i][k]);
         printf("\n");
         printf("P2 ->");
         for (k=1; k<=numReceptor; k++) printf("%ld ", p[j][k]);
         printf("\n");
         printf("O1 ->");
         for (k=1; k<=numReceptor; k++) printf("%ld ", o1[k]);
         printf("\n");
         printf("O2 ->");
         for (k=1; k<=numReceptor; k++) printf("%ld ", o2[k]);
         printf("\n");
       }
       */ 
       g_pc(o1); g_pc(o2);


       // Coloca a 1a e 2a parte de cada cromossomo em cada um dos filhos
       for (k=1; k<=numReceptor; k++) { r1[k]=0; r2[k]=0; }
       for (k=1; k<=c2; k++) {
           o1[k]=p[i][k]; o2[k]=p[j][k];
           r1[(int) p[i][k] ]= 1; r2[(int) p[j][k] ]= 1;
       }
       // Preenche o restante do cromossomo de cada filho
       k=1; v1=c2+1; v2=c2+1;
       while (k<=numReceptor) {
         if (r1[ (int)p[j][k] ]==0) {
            o1[v1]=p[j][k];
            v1=v1+1;
         }
         if (r2[ (int)p[i][k] ]==0) {
            o2[v2]=p[i][k];
            v2=v2+1;
         }
         //if (i==0) printf("i=%d j=%d v1=%d v2=%d k=%d \n", i, j, v1, v2, k);
         k=k+1;
       }  // end while
       g_pc(o1); g_pc(o2);


       // Coloca a 1a e 3a parte de cada cromossomo em cada um dos filhos
       for (k=1; k<=numReceptor; k++) { r1[k]=0; r2[k]=0; }
       for (k=1; k<c1; k++) {
           o1[k]=p[i][k]; o2[k]=p[j][k];
           r1[ (int)p[i][k] ]= 1; r2[(int) p[j][k] ]= 1;
       }
       for (k=c2+1; k<=numReceptor; k++) {
           o1[k]=p[i][k]; o2[k]=p[j][k];
           r1[ (int)p[i][k] ]= 1; r2[ (int)p[j][k] ]= 1;
       }
       // Preenche o restante do cromossomo de cada filho
       k=1; v1=c1; v2=c1;
       while (k<=numReceptor) {
         if (r1[ (int)p[j][k] ]==0) {
            o1[v1]=p[j][k];
            v1=v1+1;
         }
         if (r2[ (int)p[i][k] ]==0) {
            o2[v2]=p[i][k];
            v2=v2+1;
         }
         //if (i==0) printf("i=%d j=%d v1=%d v2=%d k=%d \n", i, j, v1, v2, k);
         k=k+1;
       }  // end while
       g_pc(o1); g_pc(o2);


       // Coloca a 2a e 3a parte de cada cromossomo em cada um dos filhos
       for (k=1; k<=numReceptor; k++) { r1[k]=0; r2[k]=0; }
       for (k=c1; k<=numReceptor; k++) {
           o1[k]=p[i][k]; o2[k]=p[j][k];
           r1[ (int)p[i][k] ]= 1; r2[ (int)p[j][k] ]= 1;
       }
       // Preenche o restante do cromossomo de cada filho
       k=1; v1=1; v2=1;
       while (k<=numReceptor) {
         if (r1[ (int)p[j][k] ]==0) {
            o1[v1]=p[j][k];
            v1=v1+1;
            //if (v1==c1) v1=c2+1;
         }
         if (r2[ (int)p[i][k] ]==0) {
            o2[v2]=p[i][k];
            v2=v2+1;
            //if (v2==c1) v2=c2+1;
         }
         //if (i==0) printf("i=%d j=%d v1=%d v2=%d k=%d \n", i, j, v1, v2, k);
         k=k+1;
       }  // end while
       g_pc(o1); g_pc(o2);


       // Coloca a 1a parte de cada cromossomo em cada um dos filhos
       for (k=1; k<=numReceptor; k++) { r1[k]=0; r2[k]=0; }
       for (k=1; k<c1; k++) {
           o1[k]=p[i][k]; o2[k]=p[j][k];
           r1[ (int)p[i][k] ]= 1; r2[ (int)p[j][k] ]= 1;
       }
       // Preenche o restante do cromossomo de cada filho
       k=1; v1=c1; v2=c1;
       while (k<=numReceptor) {
         if (r1[ (int)p[j][k] ]==0) {
            o1[v1]=p[j][k];
            v1=v1+1;
            //if (v1==c1) v1=c2+1;
         }
         if (r2[ (int)p[i][k] ]==0) {
            o2[v2]=p[i][k];
            v2=v2+1;
            //if (v2==c1) v2=c2+1;
         }
         //if (i==0) printf("i=%d j=%d v1=%d v2=%d k=%d \n", i, j, v1, v2, k);
         k=k+1;
       }  // end while
       g_pc(o1); g_pc(o2);


       // Coloca a 3a parte de cada cromossomo em cada um dos filhos
       for (k=1; k<=numReceptor; k++) { r1[k]=0; r2[k]=0; }
       for (k=c2+1; k<=numReceptor; k++) {
           o1[k]=p[i][k]; o2[k]=p[j][k];
           r1[ (int)p[i][k] ]= 1; r2[ (int)p[j][k] ]= 1;
       }
       // Preenche o restante do cromossomo de cada filho
       k=1; v1=1; v2=1;
       while (k<=numReceptor) {
         if (r1[ (int)p[j][k] ]==0) {
            o1[v1]=p[j][k];
            v1=v1+1;
            //if (v1==c1) v1=c2+1;
         }
         if (r2[ (int)p[i][k] ]==0) {
            o2[v2]=p[i][k];
            v2=v2+1;
            //if (v2==c1) v2=c2+1;
         }
         //if (i==0) printf("i=%d j=%d v1=%d v2=%d k=%d \n", i, j, v1, v2, k);
         k=k+1;
       }  // end while
       g_pc(o1); g_pc(o2);

       //printf("i=%d j=%d custoOpt=%d \n", i, j, custo_opt);
     }  // end for j
   }  // end for i
}


/* ====================================================================== */
/* =====    Crossover Partially Mapped    =============================== */

void CruzamentoPM(){
   int i, j, k, v, c1, c2;
   long int o1[MAXr], o2[MAXr], r1[MAXr], r2[MAXr];

   c1= numReceptor/3; c2= 2*c1; c1=c1 + 1; c2= c2 + 1;
   //printf("  c1=%d c2=%d \n", c1, c2);

   //cruzamento de todos os individuos distintos da população 2 a 2
   for (i=0; i<=(NPop-2); i++) {
     for (j=i+1; j<=(NPop-1); j++) {

       // Coloca a 2a parte de cada cromossomo em cada um dos filhos
       for (k=1; k<=numReceptor; k++) { r1[k]=0; r2[k]=0; }
       for (k=c1; k<=c2; k++) {
           o1[k]=p[i][k]; r1[(int) p[i][k] ]= 1;
           o2[k]=p[j][k]; r2[ (int)p[j][k] ]= 1;
       }
       // Preenche o restante do cromossomo de cada filho com o 1o e 3o bloco do outro pai
       k=1;
       while (k<=numReceptor) {
            o1[k]=p[j][k]; r1[ (int)p[j][k] ]= r1[(int) p[j][k] ] + 1;
            o2[k]=p[i][k]; r2[ (int)p[i][k] ]= r2[ (int)p[i][k] ] + 1;
            k=k+1;
            if (k==c1) k=c2+1;
       }  // end while
       // Faz o ajuste para gerar uma rota viável
       k=c1;
       while (k<=c2) {
         if (r1[ (int)p[j][k] ]==0) {
            v=1;
            while (v<=numReceptor) {
               if (r1[ o1[v] ]==2) { o1[v]=p[j][k]; v=numReceptor; }
               v=v+1;
               if (v==c1) v=c2+1;
            }
         }
         if (r2[ (int)p[i][k] ]==0) {
            v=1;
            while (v<=numReceptor) {
               if (r2[ o2[v] ]==2) { o2[v]=p[i][k]; v=numReceptor; }
               v=v+1;
               if (v==c1) v=c2+1;
            }
         }
         k=k+1;
       }  // end while
       /*
       if ((i==0) && (j==1)) {    //exibir ou não exibir
         printf("Um tipo de cruzamento \n");
         printf("DIM-");
         for (k=1; k<=numReceptor; k++) printf("%2d ",k);
         printf("\n");
         printf("P1 -");
         for (k=1; k<=numReceptor; k++) printf("%2d ", p[i][k]);
         printf("\n");
         printf("O1 -");
         for (k=1; k<=numReceptor; k++) printf("%2d ", o1[k]);
         printf("\n");
         printf("R1 -");
         for (k=1; k<=numReceptor; k++) printf("%2d ", r1[k]);
         printf("\n\n");
         printf("R2 -");
         for (k=1; k<=numReceptor; k++) printf("%2d ", r2[k]);
         printf("\n");
         printf("O2 -");
         for (k=1; k<=numReceptor; k++) printf("%2d ", o2[k]);
         printf("\n");
         printf("P2 -", j);
         for (k=1; k<=numReceptor; k++) printf("%2d ", p[j][k]);
         printf("\n");
       } //exibir ou não exibir
       */
       g_pc(o1); g_pc(o2);

       //printf("i=%d j=%d custoOpt=%d \n", i, j, custo_opt);
     }  // end for j
   }  // end for i
}


/* ====================================================================== */
/* =====    Crossover 2 Blocks   ======================================== */

void Cruzamento2B(){
   int i, j, k, v, w, c1, c2;
   long int o1[MAXr], r1[MAXr], b3[MAXr], b1[MAXr], b2[MAXr];

   c1= numReceptor/3; c2= 2*c1; c1=c1 + 1; c2= c2 + 1;
   //printf("  c1=%d c2=%d \n", c1, c2);

   //cruzamento de todos os individuos distintos da população 2 a 2
   for (i=0; i<=(NPop-2); i++) {
     for (j=i+1; j<=(NPop-1); j++) {
       for (k=1; k<=numReceptor; k++) r1[k]=0;

       // Coloca o bloco do meio do pai1 em b1
       v=1;
       for (k=c1; k<=c2; k++) { b1[v]=p[i][k]; r1[(int) p[i][k] ]= 1; v=v+1;}

       // Preenche b2 com a ordem dada de pai2 com os elementos que não estão em b1.
       v=1;
       for (k=1; k<=numReceptor; k++) {
          if (r1[ (int)p[j][k] ]==0) { b2[v]= p[j][k]; v=v+1; }
       }

       // Preenche b3 com a ordem inversa dada em b2, b3= (b2<--).
       w=1;
       while (w<v) { b3[w]= b2[v-w]; w=w+1; }

       // Gera o 1o filho O1=b1 b2;
       v=1; for (k=c1; k<=c2; k++) { o1[v]=b1[v]; v=v+1; }
       w=1; for (k=v; k<=numReceptor; k++) { o1[k]=b2[w]; w=w+1; }
       g_pc(o1);

       // Gera o 2o filho O1=b1 (b2<--);
       w=1; for (k=v; k<=numReceptor; k++) { o1[k]=b3[w]; w=w+1; }
       g_pc(o1);

       // Gera o 3o filho O1=b2 b1;
       for (k=1; k<w; k++) o1[k]=b2[k];
       v=1; for (k=w; k<=numReceptor; k++) { o1[k]=b1[v]; v=v+1; }
       g_pc(o1);

       // Gera o 4o filho O1=(b2<--) b1;
       for (k=1; k<w; k++) o1[k]=b3[k];
       g_pc(o1);

       // Coloca o bloco do meio do pai2 em b1
       for (k=1; k<=numReceptor; k++) r1[k]=0;
       v=1;
       for (k=c1; k<=c2; k++) { b1[v]=p[j][k]; r1[(int) p[j][k] ]= 1; v=v+1;}

       // Preenche b2 com a ordem dada de pai2 com os elementos que não estão em b1.
       v=1;
       for (k=1; k<=numReceptor; k++) {
          if (r1[ (int)p[i][k] ]==0) { b2[v]= p[i][k]; v=v+1; }
       }  // end while

       // Preenche b3 com a ordem inversa dada em b2, b3= (b2<--).
       w=1;
       while (w<v) { b3[w]= b2[v-w]; w=w+1; }

       // Gera o 1o filho O1=b1 b2;
       v=1; for (k=c1; k<=c2; k++) { o1[v]=b1[v]; v=v+1; }
       w=1; for (k=v; k<=numReceptor; k++) { o1[k]=b2[w]; w=w+1; }
       g_pc(o1);

       // Gera o 2o filho O1=b1 (b2<--);
       w=1; for (k=v; k<=numReceptor; k++) { o1[k]=b3[w]; w=w+1; }
       g_pc(o1);

       // Gera o 3o filho O1=b2 b1;
       for (k=1; k<w; k++) o1[k]=b2[k];
       v=1; for (k=w; k<=numReceptor; k++) { o1[k]=b1[v]; v=v+1; }
       g_pc(o1);

       // Gera o 4o filho O1=(b2<--) b1;
       for (k=1; k<w; k++) o1[k]=b3[k];
       g_pc(o1);
       //printf("i=%d j=%d custoOpt=%d \n", i, j, custo_opt);
     }  // end for j
   }  // end for i
}

/* ========================================================================= */
/* ================   Funcao Custo g para população do crosover  =========== */

void g_pcMUT(long int ss[MAXr]){
  int i, j, k, w, nigual;
  float tg;
  double soma=0, menor = 0;
  long double resultado = 0;
  int aux_doador[MAXd];
  for(int i = 0;i<numDoador;i++) aux_doador[i] = doador[i].tipo_sanguineo;
  for(int r=1;r<=numReceptor;r++){
    //printf("%d ", (int)ss[r]);
    id_doador = -1;
    id_hospital =-1;
    copia_idH = id_hospital;
    copia_idD = id_doador;
    menor = 0;
    //soma=0;
    for(int i = 0;i<numDoador;i++){
    //printf("%d ", aux_doador[i]);
    if(receptor[(int)ss[r]].tipo_sanguineo == 1){
      if(aux_doador[i] == 1){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          //printf("%lf-", soma);
        }
      }
    }
    if(receptor[(int)ss[r]].tipo_sanguineo == 2){
      if(aux_doador[i] == 2){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)ss[r]].tipo_sanguineo == 3){
      if(aux_doador[i] == 3){
        soma = 0;
       for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)ss[r]].tipo_sanguineo == 4){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }

        }
      }
    }
    if(receptor[(int)ss[r]].tipo_sanguineo == 1){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)ss[r]].tipo_sanguineo == 2){
      if(aux_doador[i] == 4){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    if(receptor[(int)ss[r]].tipo_sanguineo == 3){
       if(aux_doador[i] == 4 || aux_doador[i] == 1|| aux_doador[i] == 2){
        soma = 0;
        for(int j = 0; j<numHospital;j++){
          soma = distancia(receptor[(int)ss[r]].localizacao_x, hospital[j].localizacao_x, receptor[(int)ss[r]].localizacao_y, hospital[j].localizacao_y)+distancia(doador[i].localizacao_x, hospital[j].localizacao_x, doador[i].localizacao_y, hospital[j].localizacao_y)+hospital[j].custo;
          if(menor == 0){
            menor = soma;
            id_doador = i;
            id_hospital=j;

          }
          if(soma < menor){
            menor = soma;
            id_doador = i;
            id_hospital=j;
          }
        }
      }
    }
    //printf("%d\n", menor);
    }
    //printf("-> %lf \n", soma);
    if(menor!=0){
    //copia_idH = id_hospital;
    //copia_idD = id_doador;
    aux_doador[id_doador] = 0;
    resultado += pesoW[receptor[(int)ss[r]].prioridade] - menor;
    //return menor;
    }
    //printf("%Lf\n", resultado);
  }
  tg=resultado;
  if (tg>mut_opt) mut_opt=tg;
  if (tg>custo_opt) {
    custo_opt=tg; sopt[0]=tg;
    for (j = 1; j <= numReceptor; j++) sopt[j]=ss[j];
  }
  //inserir a solução ss em P, se ela é melhor que a pior solução de P
  i=0; k=-1;
  while (i<NPop) {
    if (tg>pc[i][0]) {
      w=i-1; k=i; i=NPop;
      while (w!=-1) {
        if (tg==pc[w][0]) {
          nigual=0;
          for (j=1; j<=numReceptor; j++) if (pc[w][j]!=ss[j]) nigual=nigual+1;
          if (nigual<(numReceptor/2)){ w=0; k=-1; }
            w=w-1;
          }
          else {
            w=-1;
          }
      } // end while
    } // end if tg<p[i][0] ...
    i=i+1;
  } // end while
  // insere ss em P na posição k
  if (k!=-1) {
    // Baixa todas as linhas da Matriz P a partir da linha k até a linha NPop-1
    for (j=(NPop-2); j>=k; j--) {
      w=j+1;
      for (i=0; i<=numReceptor; i++) pc[w][i]= pc[j][i];
    }
    // inseri o vetor ss na linha k da matriz P
    for (i = 1; i <= numReceptor; i++) pc[k][i]=ss[i]; pc[k][0]=tg;
  } // fim do if (k!=-1)
}

/* ====================================================================== */
/* ===============      Procedure Mutacao - CT regenerando PC  ========== */

void CT(long int s[MAXr]){
  int i, j, a, k, k1, k2, kk;
  long int sr[MAXr], set[MAXr];

  k=numReceptor/2;
  for (kk=1; kk<=k; kk++) { // Executa a celulas-tronco para cada valor de k=1,2,..., n/2.
    for (i=1; i<=numReceptor; i++) { set[i]=0; sr[i]=0; }
    k1=0;
    for (i=1; i<=numReceptor; i++) {
       if (set[ s[i] ]==0) {
          j=1; k1=k1+1;
          while (j<=numReceptor) {
            if (s[i]==sopt[j]) { a=j; j=numReceptor; }
            j=j+1;
          } // end while
          sr[k1]=s[i];
          set[s[i]]=set[s[i]]+1;
          j=a+1; k2=1;
          while ((j<=numReceptor) && (k2<=kk)) {
            if (set[ (int)sopt[j] ]==0) {
              k1=k1+1;
              sr[k1]=sopt[j];
              k2=k2+1;
              set[ (int)sopt[j] ]= set[ (int)sopt[j] ] +1;
            }
            j=j+1;
          }
       } //end if
    } // end for i
    g_pcMUT(sr);
  } // for kk
}


/* ====================================================================== */
/* ===============      Procedure AG   ================================== */

void AGct(){
    int i, j, k, bl, imedio, media1, media2, i_opt, z_opt, stop;
    long int solP[MAXr], m1p[3], vm1[3], vm2[3], custo;

    i=1; imedio=1; stop=0; bl=1;
    vm1[0]=0; vm1[1]=1; vm2[0]=1; vm2[1]=0;
    cross_opt=0; mut_opt=0;
    PopulacaoInicial();
    pop_opt=custo_opt;
    //printf("%d", pop_opt);
    media1=0; media2=1; i_opt=1;

    //while ((imedio<=2) && (i_opt<=10)) {
    while (stop==0) {
       //z_opt=custo_opt;

       // Aplicação do operador crossover
       custo= custo_opt;
       for (j=0; j<=(NPop-1); j++) pc[j][0]=0;
       Cruzamento1P();  // 1o tipo de crossover - com 1 ponto de corte
       //Cruzamento2P();  // 2o tipo de crossover - com 2 pontos de corte
       //CruzamentoPM();  // 3o tipo de crossover - com Partially Mapped
       //Cruzamento2B();  // 4o tipo de crossover - com 2 blocos

       if (custo!=custo_opt) {
         if (custo_opt>=cross_opt) cross_opt=custo_opt;
       }

       // O procedimento CT faz o papel do operador mutação em todas as NPop soluções da população atual
       media1=0; custo= custo_opt;
       for (j=0; j<=(NPop-1); j++) {
         for (k=1; k<=numReceptor; k++) solP[k]=p[j][k];
         media1=media1 + p[j][0];
         CT(solP);
       }
       media1=media1/NPop;

       // Verifica quantas repetições sem mudança na solução otima
       //if (z_opt==custo_opt) i_opt=i_opt+1; else i_opt=1;

       // Este procedimento faz a atualização da populaçao, os melhores NPop filhos.
       media2=0;
       for (j=0; j<=(NPop-1); j++) {
           media2=media2 + pc[j][0];
           for (k=0; k<=numReceptor; k++) {
               p[j][k]=pc[j][k];
           }
       }
       media2=media2/NPop;

       // Atualizar vetores vm1 e vm2
       for (k=0; k<2; k++) { vm1[k]=vm1[k+1]; vm2[k]=vm2[k+1]; }
       vm1[2]=media1; vm2[2]=media2;
       // Verificar 2a condição de parada
       k=0; j=0;
       while (k<3) {
         if (vm1[k]==vm2[k]) { j=j+1; }
         k=k+1;
       }
       if (vm1[0]==vm1[2]) j=2;
       if (j==2) imedio=imedio+1;

       //printf("\ni=%3d imedio=%d  -  m1=%d  m2=%d z=%d ", i, imedio, media1, media2, custo_opt);
       //printf("\n      vm1=%3d vm2=%d, vm1=%d  vm2=%d, vm1=%d vm2=%d", vm1[0], vm2[0], vm1[1], vm2[1], vm1[2], vm2[2]);
       if (imedio==2) stop=1;

       i=i+1;
    }
    printf("  %d \n", i-1);
    printf("\n\nMelhor desempenho da Populacao Inicial : %.0f \n", pop_opt);
    printf("Melhor desempenho do Operador Crossover: %.0f \n", cross_opt);
    printf("Melhor desempenho do Operador Mutacao  : %.0f \n", mut_opt);
}


/* ====================================================================== */
/* ===============      Procedure Apresenta Resultados    =============== */

void Resultado(int z)
{
    int i;
    int resultado = 0;
    int vetor_aux[MAXr];

    //printf("  %d \n", custo_opt);
    if (z==1) {
      printf("Sequencia  = ");
      printf("%.0f",sopt[0]);
      for (i = 1; i <= numReceptor; i++) {
         printf(" %.0f", sopt[i]);
         vetor_aux[i] =sopt[i];
      }
      resultado=eficiencia(vetor_aux);
      printf("\nNUMERO DE TRANSPLANTES: %d", resultado);
      printf("\n");
    }
}

/* ====================================================================== */
/* ===============      Procedure MAIN       ============================ */

int main(){
    int i, j;  double dev, time;
    FILE *trace;
    trace=fopen("reeves_DGA_KTP_cross2B_NPop=100_1011.txt", "a");

    for(i = 0;i<2;i++){
      custo_opt=0;
      input(i);
      timeused(NULL);
      AGct();
      timeused(&time);
      //dev=100*(custo_opt - zotimo[i])/(double)(zotimo[i]);
      //printf("\n%9s - %3d %.0f %.2lf seg.", recep[i], NPop, custo_opt, time);
      //printf("\n%9s - %.0f  %.0f  %.0f - z=%.0f %.2lf seg.", recep[i], pop_opt, cross_opt, mut_opt, custo_opt, time);
      Resultado(1);
      //fprintf(trace,"%s %d %.0f  %.0f  %.0f  %.0f  %.2lf\n", recep[i], NPop, pop_opt, cross_opt, mut_opt, custo_opt, time);
    }

/*
    custo_opt=0;
    input(0);
    timeused(NULL);
    AGct();
    timeused(&time);
    //dev=100*(custo_opt - zotimo[i])/(double)(zotimo[i]);
    printf("\n%9s - %3d %d %.2lf seg.", recep[0], NPop, custo_opt, time);
    printf("\n%9s - %5d  %ld  %d - z=%d %.2lf seg.", recep[0], pop_opt, cross_opt, mut_opt, custo_opt, time);
    Resultado(1);
    fprintf(trace,"%s %d %d  %ld  %d  %d  %.2lf\n", recep[0], NPop, pop_opt, cross_opt, mut_opt, custo_opt, time);
*/
    fclose(trace);
    printf("\n\nEntre com um num. inteiro para sair: ");
    scanf ("%d", &i);
    return 0;
}
