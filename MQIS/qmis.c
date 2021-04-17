#include<stdio.h>                       /* MAXIMAL QUASI-IN-STAR */
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<errno.h>
#define CBUF 100000
#define STR_IN "maximal_in_stars.txt"
//#define STR_IN "OUTPUT.txt"
/*
*         Graph input details
*                   e.g.      1         2         -->  vertex 2 incident into vertex 1 -->                   [ 1 <-------  2 ]
*
*/
struct node
{
          int vertex;
          struct node *next;
} *curr, *o_curr, *prev, *o_prev, *i_strt, *o_strt, *now, *combst , *combpr;
int ns=0;

void getGraph(int, int *, int **, char *);
void writeSet(struct node*);
void genComb(int, int);
void readCombNWriteSet(char *, struct node *, struct node *, struct node *, int, int);
void freeSet(struct node*);
int getSetSize(struct node*);

void getEdgeGraph(int, int *, int *,int **, int **,char *);


int main(void)
{
          int order, i, vertex, deg, in_cnt, out_cnt, out_lim, size_set, max_size = 0 , temp_size;
          float gma, f;
          char fname[20], temp[10];
          FILE *fcomb, *finfo, *dict;

          printf("Enter order of the Graph: ");
          scanf("%d", &order);
          printf("\nEnter value for Gamma: ");
          scanf("%f", &gma);
          int *IN_LIST[order], *OUT_LIST[order];
          int IN_DEG[order], OUT_DEG[order];
          for(i=0; i<order; i++)
          {
                    IN_LIST[i] = OUT_LIST[i] = NULL;
                    IN_DEG[i] = OUT_DEG[i] = 0;
          }


          /// ---------------------------------------------------------------------------------------------------------------


          getEdgeGraph(order, IN_DEG, OUT_DEG, IN_LIST, OUT_LIST, "Graph.txt");


          for(vertex=0; vertex<order; vertex++)             /// Star Processing
          {
                    if(IN_DEG[vertex] < 1)
                    {
                              if(OUT_DEG[vertex] > 1)           /// Cannot be side vertex
                                        IN_DEG[vertex] = -1;

                    }
                    else
                    {
                              if(gma == 0.0)
                              {
                                        now = (struct node *) malloc(sizeof (struct node));                   /// Adding center vertex to set
                                        now->vertex = vertex + 1;
                                        now->next = NULL;
                                        i_strt = prev = now;

                                        curr = (struct node*) OUT_LIST[vertex];
                                        out_cnt = 0;

                                        while(curr)
                                        {
                                                  if(IN_DEG[curr->vertex - 1] == -1 || (IN_DEG[curr->vertex - 1]+OUT_DEG[curr->vertex - 1]) != 1)
                                                  {
                                                            curr = curr->next;
                                                            continue;
                                                  }
                                                  now = (struct node *) malloc(sizeof (struct node));         /// Adding side node to set
                                                  now->vertex = curr->vertex;
                                                  now->next = NULL;
                                                  prev->next = now;
                                                  prev = now;

                                                  IN_DEG[curr->vertex - 1] = -1;
                                                  curr = curr->next;            /// Next pointer
                                                  out_cnt++;
                                        }

                                        if(out_cnt >1)
                                        {

                                                  temp_size = getSetSize(i_strt);
                                                            if(temp_size > max_size)
                                                                      max_size = temp_size;

                                                  writeSet(i_strt);             /// Write Set ...
                                                  freeSet(i_strt);
                                        }

                                        free(prev);
                                        free(i_strt);
                                        free(now);

                              }
                              else      /// gma not = 0.0
                              {
                                        now = (struct node *) malloc(sizeof (struct node));                   /// Adding center vertex to set
                                        now->vertex = vertex + 1;
                                        now->next = NULL;
                                        i_strt = prev = now;

                                        curr = (struct node*) IN_LIST[vertex];
                                        in_cnt = 0;
                                        while(curr)                   /// In degree processing
                                        {
                                                  if(IN_DEG[curr->vertex - 1] == -1 || (IN_DEG[curr->vertex - 1]+OUT_DEG[curr->vertex - 1]) != 1)
                                                  {
                                                            curr = curr->next;
                                                            continue;
                                                  }

                                                  now = (struct node *) malloc(sizeof (struct node));         /// Adding side node to set
                                                  now->vertex = curr->vertex;
                                                  now->next = NULL;
                                                  prev->next = now;
                                                  prev = now;

                                                  IN_DEG[curr->vertex - 1] = -1;
                                                  curr = curr->next;            /// Next pointer
                                                  in_cnt++;
                                        }

                                        out_cnt = 0;
                                        o_curr = (struct node*)OUT_LIST[vertex];
                                        o_strt  = o_curr;

                                        while(o_curr)                 /// Out degree processing
                                        {
                                                  if(IN_DEG[o_curr->vertex - 1] == -1 || (IN_DEG[o_curr->vertex - 1]+OUT_DEG[o_curr->vertex - 1]) != 1)
                                                  {
                                                            o_curr = o_curr->next;
                                                            continue;
                                                  }

                                                  now = (struct node *) malloc(sizeof (struct node));

                                                  if(out_cnt==0)
                                                            o_prev = now;

                                                  now->vertex = o_curr->vertex;
                                                  now->next = NULL;
                                                  o_prev->next = now;
                                                  o_prev = now;
                                                  IN_DEG[o_curr->vertex - 1] = -1;
                                                  o_curr = o_curr->next;
                                                  out_cnt++;
                                        }
                                                                                                              /// Quasi in-star processing

                                        f = (float)in_cnt / gma;
                                        out_lim = floor((double)f);
                                        out_lim = out_lim - in_cnt;             /// Permissible out degrees

                                        if(out_lim >= out_cnt)
                                        {
                                                  if(in_cnt+out_cnt > 1) //CHNG
                                                  {
                                                            prev->next = o_strt;

                                                            temp_size = getSetSize(i_strt);
                                                            if(temp_size > max_size)
                                                                      max_size = temp_size;

                                                            writeSet(i_strt);             /// Write Set

                                                            freeSet(i_strt);
                                                  }
                                        }
                                        else
                                        {
                                                  strcpy(fname,"");             /// Permutation
                                                  strcpy(temp, "");
                                                  strcpy(fname, "d");
                                                  sprintf(temp, "%d%d.txt", out_cnt, out_lim);
                                                  strcat(fname, temp);

                                                  dict = fopen(fname, "r");

                                                  if(dict == NULL)
                                                  {
                                                            fclose(dict);
                                                            genComb(out_cnt, out_lim);
                                                  }


                                                  readCombNWriteSet(fname, i_strt, o_strt, prev, out_cnt, out_lim);

                                                  fclose(dict);
                                        }

                                        free(prev);
                                        free(curr);
                                        free(i_strt);
                                        free(now);
                                        free(o_curr);
                                        free(o_prev);
                              }

                    }   /// else end
          }

          finfo = fopen("LOG_output.txt", "w");
          if(finfo == NULL)
                    printf("\nERROR!  creating output log.");
          fprintf(finfo, "Graph order: %d", order);
          fprintf(finfo, "\nGamma: %0.2f", gma);
          fprintf(finfo, "\nNo of in-star found: %d", ns);
          fprintf(finfo, "\nMax  in-star size: %d", max_size);
          fclose(finfo);

          return(1);
}         /// main end ...

/* ********************************* FUNCTION DEFINITIONS ***************************************** */

void getGraph(int order, int *arr_deg, int **arr_list, char *fname)
{
          FILE *fp;
          char buf[CBUF], *ltok;
          int i = 0, deg, vertex;

          fp = fopen(fname, "r");
          if(fp == NULL)
          {
                    printf("\nERROR!!! in opening file.");
                    exit(0);
          }

          while(fgets(buf, CBUF, fp) != NULL)
          {
                    deg=0;
                    ltok =(char*) strtok(buf, "\n\t");
                    while (ltok != NULL)
                    {
                              if(deg != 0)
                              {
                                        vertex = atoi(ltok);
                                        curr = (struct node *) malloc(sizeof (struct node));
                                        curr->next = NULL;
                                        if(deg == 1)
                                                  *arr_list =(int*) curr;
                                        else
                                                  prev->next = curr;
                                        curr->vertex = vertex;
                                        prev = curr;
                              }
                              deg++;
                              ltok =(char*) strtok (NULL, "\t\n");
                    }
                    *arr_deg = deg - 1;
                    arr_deg ++;
                    arr_list++;
                    i++;
          }
          fclose(fp);
}

void writeSet(struct node *strt)
{

          FILE *fp;
          struct node *temp;

          fp = fopen(STR_IN, "a+");
          if(fp == NULL)
          {
                    printf("\n\n\t\t :( ERROR IN WRITING SET...");
                    printf("\nError no: %d", errno);
                   exit(0);
          }

          while(strt)                             /// Writing Set
          {
                    fprintf(fp, "%d\t", strt->vertex);
                    strt = strt->next;
          }

          fprintf(fp, "\n");
          fclose(fp);
          ns++;                                   /// # of stars counter
}

void genComb(int n, int k)
{
          int i, c, j, x;
          char fname[20], temp[10];
          FILE *d;
          printf("\n\t\t\t\t\t\tComb: n = %d\tk = %d", n, k);
          strcpy(fname, "d");
          sprintf(temp, "%d%d.txt", n, k);
          strcat(fname, temp);
          d = fopen(fname, "w");
          if(d == NULL)
          {
                    printf("\nERROR! cannot create dictionary.");
                    exit(0);
          }

          for (i=0; i<(1<<n); i++)
          {
                    for (j=0,c=0; j<32; j++)
                              if (i & (1<<j))
                                        c++;
                    if (c == k)
                    {
                              x = 0;
                              for (j=0;j<32; j++)
                              {
                                        if (i & (1<<j))
                                        {
                                                  x++;
                                                  if(x == k)
                                                            fprintf (d, "%i", j);
                                                  else
                                                            fprintf (d, "%i\t", j);
                                        }
                              }
                              fprintf (d, "\n");
                    }
        }
        fclose(d);
        printf("\tfin...");
}

void readCombNWriteSet(char *fname, struct node *ii_strt, struct node *oo_strt, struct node *endptr, int n, int r)
{
          int *arr, i = 0;
          char buf[CBUF], *ltok;
          FILE *fp;

          arr = (int *) malloc(sizeof(int)*n);

          while(oo_strt)
          {
                    arr[i] = oo_strt->vertex;

                    oo_strt = oo_strt->next;
                    i++;
          }

          fp = fopen(fname, "r");
          if(fp == NULL)
          {
                    printf("\nERROR!!! while reading dictionary.");
                    exit(0);
          }

          while(fgets(buf, CBUF, fp) != NULL)
          {
                    ltok =(char*) strtok(buf, "\n\t");
                    combst = NULL;
                    while (ltok != NULL)
                    {
                              now = (struct node *) malloc(sizeof(struct node));
                              now->vertex = arr[atoi(ltok)];
                              now->next = NULL;
                              if(combst == NULL)
                                        combst = now;
                              else
                                        combpr->next = now;
                              combpr = now;

                              ltok =(char*) strtok (NULL, "\t\n");
                    }
                    endptr->next = combst;
                    writeSet(ii_strt);
                    freeSet(combst);
          }
          free(arr);
          freeSet(ii_strt);
          fclose(fp);
}

void freeSet(struct node *f)
{
          struct node *temp;
          while(f)                             /// Freeing memory.
          {
                    temp = f;
                    f = f->next;
                    free(temp);
          }
}

int getSetSize(struct node *x)
{
          int i=0;
          while(x)
          {
                    x  = x->next;
                    i++;
          }
          return i;
}

/* NEW Function definition*/
void getEdgeGraph(int order, int *arr_in_deg, int *arr_out_deg, int **arr_in_list, int **arr_out_list, char *fname)
{
          FILE *fp;
          char buf[CBUF], *ltok;
          int i = 0, deg, miRNA, Gene;
          struct node *tmp, *prv;

          fp = fopen(fname, "r");
          if(fp == NULL)
          {
                    printf("\nERROR!!! in opening file.");
                    exit(0);
          }

          while(fgets(buf, CBUF, fp) != NULL)
          {
                    ltok =(char*) strtok(buf, "\n\t");
                    miRNA =  atoi(ltok);
                    ltok =(char*) strtok(NULL, "\t\n");
                    Gene = atoi(ltok);

                    if(arr_in_list[miRNA-1] == NULL)                  /// In-List Processing
                    {
                              curr = (struct node *) malloc(sizeof (struct node));
                              arr_in_list[miRNA-1] = (int*) curr;

                              curr->vertex = Gene;
                              curr->next = NULL;
                              arr_in_deg[miRNA-1]++;
                    }
                    else
                    {
                              tmp = (struct node*) arr_in_list[miRNA-1];
                              while(tmp)
                              {
                                        prv = tmp;
                                        tmp = tmp->next;
                              }
                              curr = (struct node *) malloc(sizeof (struct node));

                              curr->next = NULL;
                              curr->vertex = Gene;
                              prv->next=curr;
                              arr_in_deg[miRNA-1]++;
                    }

                    if(arr_out_list[Gene-1] == NULL)                  /// Out-List Processing
                    {
                              curr = (struct node *) malloc(sizeof (struct node));
                              arr_out_list[Gene-1] = (int*) curr;

                              curr->vertex = miRNA;
                              curr->next = NULL;
                              arr_out_deg[Gene-1]++;
                    }
                    else
                    {
                              tmp = (struct node*) arr_out_list[Gene-1];
                              while(tmp)
                              {
                                        prv = tmp;
                                        tmp = tmp->next;
                              }
                              curr = (struct node *) malloc(sizeof (struct node));

                              curr->next = NULL;
                              curr->vertex = miRNA;
                              prv->next=curr;
                              arr_out_deg[Gene-1]++;
                    }
          }
          fclose(fp);
}
