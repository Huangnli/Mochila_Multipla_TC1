/* mochila_multipla.c
codigo exemplo de uso do glpk para MIP (mixed integer programming)
By: Edna A. Hoshino

codigo das heuristicas
By: Adriano Rodrigues Alves
By: Julio Huang
By: Luís Fernando Leite França
By: Patrick Escorsi Silva (F)

Problema da mochila multipla: 
dados um conjunto de itens I={1,2,...,n} e um conjunto de mochilas K={1,2,..., k}
em que cada item i tem um peso pi e um valor vi
e cada mochila k tem uma capacidade Ck,
encontrar um subconjunto dos itens que podem ser transportados nas mochilas de modo
a maximizar o valor transportado. Ou seja, a soma do valor dos itens escolhidos
deve ser maxima e a soma dos pesos dos itens escolhidos em cada mochila nao pode
ultrapassar a sua capacidade.

Modelo de PLI para mochila multipla:

max v1(x11 + x12 + ... x1k) + v2(x21+x22+...+x2k) + ... + vn (xn1+xn2+...+xnk)
s.a.:
   p1x11 + p2x21 + ... pnxn1 <= C1
   p1x12 + p2x22 + ... pnxn2 <= C2
   ...
   p1x1k + p2x2k + ... pnxnk <= Ck
   
   x11 + x12 + ... + x1k <= 1
   x21 + x22 + ... + x2k <= 1
   ...
   xn1 + xn2 + ... + xnk <= 1



variaveis: xij binarias = 1 se o item i eh transportado na mochila j

*/

#include <stdio.h>
#include <stdlib.h>
#include <glpk.h>
#include <time.h>
#include <string.h>

#define EPSILON 0.000001

#ifdef DEBUG
#define PRINTF(...) printf(__VA_ARGS__)
#else
#define PRINTF(...)
#endif

typedef struct
{
  int num;      /* numero do item */
  double valor; /* valor do item */
  int peso;     /* peso do item */
  int index;    /* mochila ao qual o item pertence */
} Titem;

typedef struct
{
  int n;       /* total de itens */
  Titem *item; /* conjunto dos itens */
  int k;       /* total de mochilas */
  int *C;      /* capacidade das mochilas */
} Tinstance;

int carga_lp(glp_prob **lp, Tinstance I);
int carga_instancia(char *filename, Tinstance *I);
void free_instancia(Tinstance I);
int RandomInteger(int low, int high);
int comparador(const void *valor1, const void *valor2);
int comparador_num(const void *num1, const void *num2);
double guloso(Tinstance I);
double random_heuristica(Tinstance I);
double guloso_melhorada(Tinstance I);
void troca(Titem *a, Titem *b);
double heuristica(Tinstance I, int tipo);
double otimiza_PLI(Tinstance I, int tipo, double *x);
void gerar_arquivo_sol(char *filename, double z, Tinstance I);
void gerar_arquivo_out(char *filename, int tipo, double z, double tempo);
double destroy_rins(Tinstance I, double z, double *x);
double repair_rins(Tinstance I, double z, double *x);

/* carrega o modelo de PLI nas estruturas do GLPK */
int carga_lp(glp_prob **lp, Tinstance I)
{
  int *ia, *ja, nrows, ncols, i, k, row, col, nz;
  double *ar;
  char name[80]; // nome da restricao

  nrows = I.k + I.n; // 1 restricao de capacidade para cada mochila + 1 para cada item
  ncols = I.n * I.k;

  // Aloca matriz de coeficientes
  ia = (int *)malloc(sizeof(int) * (I.n * I.k * 2 + 1));
  ja = (int *)malloc(sizeof(int) * (I.n * I.k * 2 + 1));
  ar = (double *)malloc(sizeof(double) * (I.n * I.k * 2 + 1));

  // Cria problema de PL
  *lp = glp_create_prob();
  glp_set_prob_name(*lp, "mochila_multipla");
  glp_set_obj_dir(*lp, GLP_MAX);

  // Configura restricoes
  glp_add_rows(*lp, nrows);

  // criar uma restricao de capacidade para cada mochila
  row = 1;
  for (k = 0; k < I.k; k++)
  {
    sprintf(name, "capacidade_Mochila_%d", row); /* nome das restricoes */
    glp_set_row_name(*lp, row, name);
    glp_set_row_bnds(*lp, row, GLP_UP, 0.0, I.C[row - 1]);
    row++;
  }

  // criar uma restricao de unicidade para cada item
  for (i = 0; i < I.n; i++)
  {
    sprintf(name, "unicidade_%d", row); /* nome das restricoes */
    glp_set_row_name(*lp, row, name);
    glp_set_row_bnds(*lp, row, GLP_UP, 0.0, 1.0);
    row++;
  }
  // Configura variaveis
  glp_add_cols(*lp, ncols);

  col = 1;
  for (k = 1; k <= I.k; k++)
  {
    for (i = 0; i < I.n; i++)
    {
      name[0] = '\0';
      sprintf(name, "x%d_%d", i + 1, k); /* as variaveis referem-se `as variaveis xi_k para cada item i e cada mochila k */
      glp_set_col_name(*lp, col, name);
      glp_set_col_bnds(*lp, col, GLP_DB, 0.0, 1.0);
      glp_set_obj_coef(*lp, col, I.item[i].valor);
      glp_set_col_kind(*lp, col, GLP_BV); // especifica que a variaval xik eh binaria
      col++;
    }
  }

  // Configura matriz de coeficientes ...
  nz = 1;
  //coeficientes para as restricoes de cada item
  for (k = 1; k <= I.k; k++) //para cada mochila
  {
    for (i = 1; i <= I.n; i++) //para cada item
    {
      // restr de capacidade
      ia[nz] = k;                  // linha (indice da restricao)
      ja[nz] = (k - 1) * I.n + i;  // coluna (indice da variavel = xik)
      ar[nz] = I.item[i - 1].peso; // coeficiente da matriz de coeficientes na linha e coluna
      nz++;
    }
  }

  // coeficientes para as restricoes de unicidade
  for (i = 1; i <= I.n; i++) //para cada item
  {
    for (k = 1; k <= I.k; k++) //para cada mochila
    {
      // restr de unicidade
      ia[nz] = I.k + i;
      ja[nz] = (k - 1) * I.n + i;
      ar[nz] = 1.0;
      nz++;
    }
  }

  // Carrega PL
  glp_load_matrix(*lp, nz - 1, ia, ja, ar);

  // libera memoria
  free(ia);
  free(ja);
  free(ar);
  return 1;
}

/* carrega os dados da instancia de entrada */
int carga_instancia(char *filename, Tinstance *I)
{
  FILE *fin;
  int i, capacidade, item, peso;
  double valor;

  fin = fopen(filename, "r");
  if (!fin)
  {
    printf("\nProblema na abertura do arquivo %s\n", filename);
    return 0;
  }

  fscanf(fin, "%d %d", &(I->n), &(I->k));

  // aloca memória
  (*I).C = (int *)malloc(sizeof(int) * ((*I).k));
  (*I).item = (Titem *)malloc(sizeof(Titem) * ((*I).n));

  for (i = 0; i < (*I).k; i++)
  {
    fscanf(fin, "%d", &capacidade);
    (*I).C[i] = capacidade;
  }

  for (i = 0; i < (*I).n; i++)
  {
    fscanf(fin, "%d %d %lf", &item, &peso, &valor);
    if (item < 1 || item > (*I).n)
    {
      fclose(fin);
      return 0;
    }
    (*I).item[i].num = item;
    (*I).item[i].peso = peso;
    (*I).item[i].valor = valor;
  }

#ifdef DEBUG
  printf("n=%d k=%d\n", (*I).n, (*I).k);
  for (i = 0; i < (*I).k; i++)
  {
    printf("C[%d]=%d\n", i + 1, (*I).C[i]);
  }
  for (i = 0; i < (*I).n; i++)
  {
    printf("p[%d]=%d e v[%d]=%lf\n", (*I).item[i].num, (*I).item[i].peso, (*I).item[i].num, (*I).item[i].valor);
  }
#endif
  fclose(fin);
  return 1;
}

/* libera memoria alocada pelo programa para guardar a instancia */
void free_instancia(Tinstance I)
{
  free(I.item);
  free(I.C);
}

/* sorteia um numero aleatorio entre [low,high] */
int RandomInteger(int low, int high)
{
  int k;
  double d;

  d = (double)rand() / ((double)RAND_MAX + 1);
  k = d * (high - low + 1);
  return low + k;
}

/* resolve o problema de PLI usando o GLPK */
double otimiza_PLI(Tinstance I, int tipo, double *x)
{
  glp_prob *lp;
  double z, valor;
  glp_smcp param_lp;
  glp_iocp param_ilp;
  int status, i, k;

  // desabilita saidas do GLPK no terminal
  glp_term_out(GLP_OFF);

  // carga do lp
  carga_lp(&lp, I);

  // configura simplex
  glp_init_smcp(&param_lp);
  param_lp.msg_lev = GLP_MSG_ON;

  // configura optimizer
  glp_init_iocp(&param_ilp);
  param_ilp.msg_lev = GLP_MSG_ALL;
  param_ilp.tm_lim = 1000; // tempo limite do solver de PLI
  param_ilp.out_frq = 100;

  // Executa Solver de PL
  glp_simplex(lp, &param_lp); // resolve o problema relaxado
  if (tipo == 2)
  {
    glp_intopt(lp, &param_ilp); // resolve o problema inteiro
  }

  if (tipo == 2)
  {
    status = glp_mip_status(lp);
    PRINTF("\nstatus=%d\n", status);
  }
  // Recupera solucao
  if (tipo == 1)
    z = glp_get_obj_val(lp);
  else
    z = glp_mip_obj_val(lp);

  for (k = 0; k < I.k; k++)
  {
    for (i = 0; i < I.n; i++)
    {
      if (tipo == 1)
        valor = glp_get_col_prim(lp, k * I.n + i + 1); // recupera o valor da variavel xik relaxado (continuo)
      else
        valor = glp_mip_col_val(lp, k * I.n + i + 1); // recupera o valor da variavel xik
      if (valor > EPSILON)
        PRINTF("x%d_%d = %.2lf\n", I.item[i].num, k + 1, valor);
      x[k * I.n + i] = valor;
    }
  }

#ifdef DEBUG
  // Grava solucao e PL
  PRINTF("\n---LP gravado em mochila.lp e solucao em mochila.sol");
  glp_write_lp(lp, NULL, "mochila.lp");
  if (tipo == 1)
    glp_print_sol(lp, "mochila.sol");
  else
    glp_print_mip(lp, "mochila.sol");
#endif
  // Destroi problema
  glp_delete_prob(lp);
  return z;
}

// Função auxiliar de comparacao para o qsort
int comparador(const void *valor1, const void *valor2)
{
  if ((*(Titem *)valor1).valor > (*(Titem *)valor2).valor)
  {
    return -1;
  }
  else if ((*(Titem *)valor1).valor == (*(Titem *)valor2).valor)
  {
    return 0;
  }
  else
  {
    return 1;
  }
}

// Função auxiliar de comparacao para o qsort
int comparador_num(const void *num1, const void *num2)
{
  if ((*(Titem *)num1).num > (*(Titem *)num2).num)
  {
    return 1;
  }
  else if ((*(Titem *)num1).num == (*(Titem *)num2).num)
  {
    return 0;
  }
  else
  {
    return -1;
  }
}

//Primeira heuristica implementada pelo grupo
double guloso(Tinstance I)
{
  double z = 0.0; // Melhor resposta;
  int j;          // Indice da mochila

  qsort(I.item, I.n, sizeof(Titem), comparador); //Ordenacao da lista de itens pelo valor de cada item

  // Inicializa os index do item com valor zero
  for (int i = 0; i < I.n; i++)
    I.item[i].index = 0;

  // Percorre toda lista de itens
  for (int i = 0; i < I.n; i++)
  {
    j = 0;
    while (j < I.k) // Tenta colocar o item em alguma mochila
    {
      if (I.item[i].peso <= I.C[j]) // Verifica se o peso do item não
      {                             // ultrapassa a capacidade da mochila
        I.item[i].index = j + 1;
        I.C[j] -= I.item[i].peso;
        z += I.item[i].valor;
        j = I.k;
      }
      else
        j++;
    }
  }

  return z;
}

// Função que troca dois itens da lista de itens de lugar
void troca(Titem *a, Titem *b)
{
  Titem aux;

  aux.num = (*a).num;
  aux.valor = (*a).valor;
  aux.peso = (*a).peso;
  aux.index = (*a).index;

  (*a).num = (*b).num;
  (*a).valor = (*b).valor;
  (*a).peso = (*b).peso;
  (*a).index = (*b).index;

  (*b).num = aux.num;
  (*b).valor = aux.valor;
  (*b).peso = aux.peso;
  (*b).index = aux.index;
}

//Segunda heuristica implementada pelo grupo
double random_heuristica(Tinstance I)
{
  double z = 0.0; // Melhor resposta
  int i;          // Item escolhido aleatoriamente
  int j;          // Indice da mochila
  int n = I.n;    // Numero de itens restantes na lista de itens

  for (int k = 0; k < I.n; k++)
    I.item[k].index = 0;

  srand(time(NULL));

  while (n > 0)
  {
    i = RandomInteger(0, n - 1); // Escolhe um item aleatoriamente
    j = 0;

    while (j < I.k) // Verifica em qual mochila colocar o item
    {
      if (I.item[i].peso <= I.C[j])
      {
        I.item[i].index = j + 1;
        I.C[j] -= I.item[i].peso;
        z += I.item[i].valor;
        j = I.k;
      }
      else
        j++;
    }

    troca(&I.item[i], &I.item[n - 1]); // Troca o item escolhido pelo ultimo item
    n--;                               // Tira o ultimo item da lista para ele nao ser escolhido novamente
  }

  return z;
}

double destroy_rins(Tinstance I, double z, double *x)
{
  for (int i = 0; i < I.n; i++)
  {
    if (I.item[i].index != 0)
    {
      // printf("x: %f\n", x[(I.item[i].index - 1) * I.n + i]);
      if (x[(I.item[i].index - 1) * I.n + i] < 1.0)
      {                                               // item vai ser tirado da mochila
        I.C[(I.item[i].index - 1)] += I.item[i].peso; // ajusta o peso preenchido na mochila
        z -= I.item[i].valor;                         // ajusta o valor da solução
        printf("mochila: %d item: %d\n", I.item[i].index, I.item[i].num);
        I.item[i].index = 0; // retira o item da mochila
      }
    }
  }
  return z;
}

void destroy(Tinstance I)
{
  int maior_peso = 0;
  int maior_num;

  for (int i = 0; i < I.k; i++)
  {
    for (int j = 0; j < I.n; j++)
    {
      if (I.item[j].index == i + 1)
      { // Verifica se o item pertence a mochila i
        if (I.item[j].peso > maior_peso)
        { // Pega o maior peso dentre os itens da mochila i
          maior_peso = I.item[j].peso;
          maior_num = j;
        }
      }
    }
    I.C[i] += maior_peso;
    I.item[maior_num].index = 0;
  }
}

double repair_rins(Tinstance I, double z, double *x)
{
  int j;
  for (int i = 0; i < I.n; i++)
  {
    j = I.k - 1;
    while (j >= 0) // Tenta colocar o item em alguma mochila
    {
      if (I.item[i].peso <= I.C[j] && I.item[i].index == 0) // Verifica se o peso do item não
      {                                                     // ultrapassa a capacidade da mochila e se o item ainda não
        I.item[i].index = j + 1;                            // foi escolhido
        I.C[j] -= I.item[i].peso;
        z += I.item[i].valor;
        j = I.k;
        // printf("oiii %lf\n", z2);
      }
      else
        j--;
    }
  }
  return z;
}

//Heuristica guloso melhorada
double guloso_melhorada(Tinstance I)
{
  double z1 = 0.0; // Melhor resposta do relaxado
  double z2 = 0.0; // Melhor resposta do guloso
  double soma = 0.0;
  double *x, *x2;
  int tipo = 1;
  int j;
  Tinstance I_PLI;

  // aloca memória
  (I_PLI).C = (int *)malloc(sizeof(int) * ((I).k));
  (I_PLI).item = (Titem *)malloc(sizeof(Titem) * ((I).n));

  x = (double *)malloc(sizeof(double) * (I.n * I.k));
  z1 = otimiza_PLI(I, tipo, x);
  z2 = guloso(I);
  qsort(I.item, I.n, sizeof(Titem), comparador_num);

  /*for (int k = 0; k < I.k; k++) { // imprime os itens pegos pela heuristica relaxada
    for (i = 0; i < I.n; i++){
        printf("%lf\n", x[k*I.n + i]);
    }
  }*/

  z2 = destroy_rins(I, z2, x);

  // Percorre a lista de itens
  // repair:
  z2 = repair_rins(I, z2, x);

  // destroi novamente
  z2 = destroy_rins(I, z2, x);

  //destroy(I);

  //for (int i = 0; i < I.n; i++)
  //{
  //  if (I.item[i].index != 0)
  //  {
  //    printf("\n%d %f\n", I.item[i].num, I.item[i].valor);
  //  }
  //}

  j = 0;
  for (int i = 0; i < I.n; i++)
  {
    if (I.item[i].index == 0)
    {                                    // item nao foi pego
      I_PLI.item[j].num = I.item[i].num; // I_PLI copia o item nao pego
      I_PLI.item[j].valor = I.item[i].valor;
      I_PLI.item[j].peso = I.item[i].peso;
      I_PLI.item[j].index = I.item[i].index;
      j++;
    }
  }

  I_PLI.n = j + 1; //Total de itens que ainda pode ser levado
  I_PLI.k = I.k;   //Total de mochila
  for (int i = 0; i < I.k; i++)
  {
    I_PLI.C[i] = I.C[i];
    //printf("\n%d %d\n", I_PLI.C[i], I.C[i]);
  }

  x2 = (double *)malloc(sizeof(double) * (I_PLI.n * I_PLI.k));
  z1 = otimiza_PLI(I_PLI, 2, x2);

  for (int i = 0; i < I_PLI.n * I_PLI.k; i++)
  {
    if (x2[i] == 1)
    {
      printf("\nx2= %d\n", i % I_PLI.n + 1);
    }
  }

  printf("\nz1: %f\n", z1);

  for (int i = 0; i < I_PLI.n; i++)
  {
    for (int j = 0; j < I_PLI.k; j++)
    {
      if (x2[j * I_PLI.n + i] == 1.0) // item levado
      {
        I.item[I_PLI.item[i].num - 1].index = (j + 1);
        // printf("\nnum: %d %f %d %f\n", I_PLI.item[i].num, I_PLI.item[i].valor, I.item[I_PLI.item[i].num - 1].num, I.item[I_PLI.item[i].num - 1].valor);
        break;
      }
    }
  }

  // for (int i = 0; i < I_PLI.n * I_PLI.k; i++)
  // {
  //   if (x2[i] == 1.0) // item levado
  //   {
  //     printf("\ni: %d\n", i%I_PLI.n);
  //   }
  // }

  for (int i = 0; i < I.n; i++)
  {
    if (I.item[i].index != 0)
    {
      soma += I.item[i].valor;
      // printf("\nnum2:%d %f\n", I.item[i].num, I.item[i].valor);
      // printf("\nnum2: %d\n", I.item[i].num);
    }
  }

  // libera memoria alocada
  free_instancia(I_PLI);
  free(x2);
  free(x);
  return soma;
}

/* heuristica a ser implementada */
double heuristica(Tinstance I, int tipo)
{
  double z = 0.0;

  if (tipo == 3)
  {
    z = guloso(I);
  }
  else if (tipo == 4)
  {
    z = random_heuristica(I);
  }
  else
  {
    z = guloso_melhorada(I);
  }

  return z;
}

void gerar_arquivo_sol(char *filename, double z, Tinstance I)
{
  FILE *arquivo_saida;
  char nomeArquivo[64];
  int soma;

  strcpy(nomeArquivo, filename);
  strcat(nomeArquivo, ".sol");

  arquivo_saida = fopen(nomeArquivo, "w");

  fprintf(arquivo_saida, "%.0lf %d\n", z, I.k);

  for (int j = 1; j <= I.k; j++)
  {
    soma = 0;
    for (int i = 0; i < I.n; i++)
    {
      if (I.item[i].index == j)
        soma += 1;
    }
    fprintf(arquivo_saida, "mochila %d %d\n", j, soma);
    if (soma != 0)
    {
      for (int i = 0; i < I.n; i++)
      {
        if (I.item[i].index == j)
          fprintf(arquivo_saida, "%d ", I.item[i].num);
      }
    }
    fprintf(arquivo_saida, "\n");
  }

  fclose(arquivo_saida);
}

void gerar_arquivo_out(char *filename, int tipo, double z, double tempo)
{
  FILE *arquivo_saida;
  char nomeArquivo[100];
  const char *gerador;
  int status;
  char UB[10];

  const char *implementacao1 = "-1";
  const char *implementacao2 = "-2";

  const char *heuristica1 = "-1";
  const char *heuristica2 = "-2";

  strcpy(nomeArquivo, filename);

  if (tipo < 3)
  {
    strcat(nomeArquivo, implementacao1);
    strcat(nomeArquivo, "-0");
    if (tipo == 1)
      gerador = "1:relaxação";
    else
      gerador = "2:branch-and-bound";
  }
  else
  {
    strcat(nomeArquivo, implementacao2);
    if (tipo == 3)
    {
      strcat(nomeArquivo, heuristica1);
      gerador = "3:heuristica gulosa";
    }
    else
    {
      strcat(nomeArquivo, heuristica2);
      gerador = "4:heuristica aleatória";
    }
    status = 10;
  }

  strcat(nomeArquivo, ".out");

  arquivo_saida = fopen(nomeArquivo, "w");

  if (tipo < 3)
    sprintf(UB, "%.0lf", z);
  else
    sprintf(UB, " ");

  fprintf(arquivo_saida, "%s;%s;%lf;%.0lf;%s;%d", filename, gerador, tempo, z, UB, status);

  fclose(arquivo_saida);
}

/* programa principal */
int main(int argc, char **argv)
{
  double z, *x, tempo;
  clock_t antes, agora;
  int tipo;
  Tinstance I;

  // checa linha de comando
  if (argc < 3)
  {
    printf("\nSintaxe: mochila <instancia.txt> <tipo>\n\t<tipo>: 1 = relaxacao linear, 2 = solucao inteira\n");
    exit(1);
  }

  // ler a entrada
  if (!carga_instancia(argv[1], &I))
  {
    printf("\nProblema na carga da instância: %s", argv[1]);
    exit(1);
  }

  tipo = atoi(argv[2]);
  if (tipo < 1 || tipo > 5)
  {
    printf("Tipo invalido\nUse: tipo=1 (relaxacao linear), 2 (solucao inteira), 3 (heuristica gulosa), 4 (heuristica aleatoria), 5 (heuristica gulosa melhorada)\n");
    exit(1);
  }

  antes = clock();
  if (tipo < 3)
  {
    // aloca memoria para a solucao
    x = (double *)malloc(sizeof(double) * (I.n * I.k));
    z = otimiza_PLI(I, tipo, x);
    free(x);
  }
  else
  {
    // heuristica
    z = heuristica(I, tipo);
  }
  agora = clock();

  PRINTF("Valor da solucao: %lf\tTempo gasto=%lf\n", z, ((double)agora - antes) / CLOCKS_PER_SEC);

  printf("%s;%d;%d;%d;%.0lf;%lf\n", argv[1], tipo, I.n, I.k, z, ((double)agora - antes) / CLOCKS_PER_SEC);

  tempo = ((double)agora - antes) / CLOCKS_PER_SEC;

  if (tipo > 2)
  {
    gerar_arquivo_sol(argv[1], z, I);
    gerar_arquivo_out(argv[1], tipo, z, tempo);
  }

  // libera memoria alocada
  free_instancia(I);

  return 0;
}

/* eof */
