/************************************************************************************************
*Implementation of reconfiguration of distribution system using particle swarm optimization 	*
*written by Vasudevan.B (vasu@ee.iitkgp.ernet.in)						*
*Last Modified on 27.05.15									* 
/************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define NBus 69
#define NPART 05 
#define MITER  02
#define c1 0.3
#define c2 0.45
#define wmax .9
#define wmin .4
#define Hi 100
#define Low 0
#define ntie 1
#define gbest 1.25

//Structure pointer
typedef struct
{
	int fb;
	int tb;
	int status;
}Feederdata;

typedef struct
{
	int bus_no;
	double Pg;
	double Qg;
	double Pl;
	double Ql;
}Busdata;

typedef struct
{
	double fit,up_fit;
	double v_int,v_fin;		//v_int, v_fin are the inital velocity and final velocity.
	int *p, *updated_p;
	double x;
	double avg_ls;		//average load served
}abc;

//Function declaration
void Adjacency_matrix_formation();
void read_generation_data();
double fitness_function();
void Risk_factor_calculation();
void pso();
double rand_number();
int random_bus();


//Pointer declaration
Feederdata *fd; Busdata *bd;
double *generation,*load_weight,*gen_node_value;
int *gen_node,**adj;



//File Pointer
FILE *fpdt,*fgen,*fadm,*fwei,*test;


//Variable declaration
int lweight_choice,No_gen,i,j,k,l,ref_mat[3]; int iter = 1;
double risk_factor,load_served,Tot_load_served,sum_risk_factor,pbest,new_gbest;



int main(int argc,char *argv[])
{
  
  int i,j,**graph1;
  printf("\n********************************************************************************\n");
  printf("\t\tProgram to determine optimal Tie switch using PSO\n");
  printf("********************************************************************************\n");
  Adjacency_matrix_formation();	//Initial calcualtion of adjacency matrix formation
  read_generation_data();
  
  graph1 = malloc((NBus+1) * sizeof(int *));
  if(graph1 == NULL)
  {
    printf("out of memory\n");
    exit(1);
  }
  for(i = 1; i <= (NBus+1); i++)
  {
    graph1[i] = malloc((NBus+1) * sizeof(int));
    if(graph1[i] == NULL)
    {
	printf("out of memory\n");
	exit(1);
    }
  }

  for (i = 1; i <= NBus; i++)		//Copying the adjacency matrix to graph for further processing 
  {
    for (j = 1; j <= NBus; j++)
      {
	graph1[i][j] = adj[i][j];
      }
  }
  
  risk_factor = 0;
  test = fopen("Output/Testing.pdt", "w");
  if (test == NULL)
  {
     printf("\nError in creating output file\n");
  }
  /*printf("\n-----------------------------------------");
  printf("\nBase Case without any tie switch\n");
  printf("-----------------------------------------\n");*/
  fitness_function();

  for(i = 1; i <= (NBus+1); i++)
  {
     free(graph1[i]);
     graph1[i] == NULL;
  }
  free(graph1);
  graph1 == NULL;	
  
  //pso();
  fclose(test);
   
	  

  //Freeing memory
  free(fd); free(bd);   
  for(i = 1; i <= (NBus+1); i++)
  {
     free(adj[i]);
     adj[i] == NULL;
  }
  free(adj);
  free(generation); free(gen_node_value); free(load_weight);
  free(gen_node); 
}

void pso()
{
   abc particle[NPART + 1];   
   FILE *fpso;
   int oneswitch[3], twoswitch[5], threeswitch[7], fourswitch[9], fiveswitch[11],p,a;
   int part_no,**graph1;
   double tempa, sum_fit,w,xx, yy, new_gbest,t1, t2;
   w = 0; xx = 0; yy = 0;

   for (i = 1; i <= NPART; i++)
   {
	particle[i].p = (int*)malloc(NBus+1 * sizeof(int));
	particle[i].updated_p = (int*)malloc(NBus+1 * sizeof(int));
   }

   fpso = fopen("Output/PSO_Output.pdt", "w");
   if (fpso == NULL)
   {
	printf("\nError in creating output file\n");
   }

  /***************************************************************************
   *  Step 1 : creating initial particle using random number		     *
   ***************************************************************************/
   fprintf(fpso, "\n*******************************************************************");
   fprintf(fpso, "\nInitial particles\t velocity(v_int)\n");
   fprintf(fpso, "*******************************************************************\n\n");
   
   for (i = 1; i <= NPART; i++)	//Random bus number and its velocity 
   {
	for (j = 1; j <= (2 * ntie); j++)
	{
		particle[i].p[j] = random_bus();
	}
	//Condition to check whether both from and to buses are to be same.
	for (k = 1; k <= ((2 * ntie) - 1);)
	{
		if (particle[i].p[k] == particle[i].p[k + 1])
		{
	label2:
			particle[i].p[k + 1] = random_bus();
			if (particle[i].p[k + 1] == particle[i].p[k])
			{
				goto label2;
			}
		}
		k = k + 2;
	}
	//Printing purpose only
	for (j = 1; j <= (2 * ntie); j++)
	{
		fprintf(fpso, "%d\t", particle[i].p[j]);
	}
	//Initial Velocity
	particle[i].v_int = rand_number();
	fprintf(fpso, "\t%lf", particle[i].v_int);
	fprintf(fpso, "\n");
   }
  /***************************************************************************
   *  		Step 2 : Run objective fucntion and find pbest		     *
   ***************************************************************************/
   do
   {
      if (ntie == 1)		//One switch combintions
      {
	fprintf(fpso, "\n\n**************************\nIteration = %d\n**************************\n\n", iter);
	fprintf(fpso, "------------------------------------------\nParticle\t\t\tAverage Load Served\n------------------------------------------\n");
	for (p = 1; p <= NPART; p++)	//NPART
	{
	   //Memory Allocation
   	   graph1 = malloc((NBus+1) * sizeof(int *));
           if(graph1 == NULL)
           {
              printf("out of memory\n");
              exit(1);
           }
           for(i = 1; i < (NBus+1); i++)
           {
             graph1[i] = malloc((NBus+1) * sizeof(int));
             if(graph1[i] == NULL)
             {
 	        printf("out of memory\n");
 	        exit(1);
             }
           }

           for (i = 1; i <= NBus; i++)		//Copying the adjacency matrix to graph for further processing 
	   {
		for (j = 1; j <= NBus; j++)
		{
	           graph1[i][j] = adj[i][j];
		}
	   }
           oneswitch[0] = particle[p].p[1]; oneswitch[1] = particle[p].p[2];
        
      	  graph1[oneswitch[0]][oneswitch[1]] = 1;
	  graph1[oneswitch[1]][oneswitch[0]] = 1;

	  fitness_function();
	  fprintf(fpso, "\nparticle[%2d]\t\t\t%lf", p, (Tot_load_served / (NBus - 2)));
	  particle[p].fit = (Tot_load_served / (NBus - 2));
	  //Freeing Memory
   	  for(i = 1; i <= (NBus+1); i++)
          {
      	      free(graph1[i]);
              graph1[i] == NULL;
          }
   	  free(graph1); graph1 == NULL;
	}
      }

        if (ntie == 2)		//Two switch combintions
        {
	  fprintf(fpso, "\n\n**************************\nIteration = %d\n**************************\n\n", iter);
	  fprintf(fpso, "------------------------------------------\nParticle\t\t\tAverage Load Served\n------------------------------------------\n");
	  for (p = 1; p <= NPART; p++)	//NPART
	  {
	     //Memory Allocation
   	     graph1 = malloc((NBus+1) * sizeof(int *));
             if(graph1 == NULL)
             {
                printf("out of memory\n");
                exit(1);
             }
             for(i = 1; i < (NBus+1); i++)
             {
               graph1[i] = malloc((NBus+1) * sizeof(int));
               if(graph1[i] == NULL)
               {
 	          printf("out of memory\n");
 	          exit(1);
               }
             }

             for (i = 1; i <= NBus; i++)		//Copying the adjacency matrix to graph for further processing 
	     {
		  for (j = 1; j <= NBus; j++)
		  {
	             graph1[i][j] = adj[i][j];
		  }  
	     }
             twoswitch[0] = particle[p].p[1]; twoswitch[1] = particle[p].p[2];
             twoswitch[2] = particle[p].p[3]; twoswitch[3] = particle[p].p[4];
        
      	     graph1[twoswitch[0]][twoswitch[1]] = 1;
	     graph1[twoswitch[1]][twoswitch[0]] = 1;
             graph1[twoswitch[2]][twoswitch[3]] = 1;
	     graph1[twoswitch[3]][twoswitch[2]] = 1;

	     fitness_function();
	     fprintf(fpso, "\nparticle[%2d]\t\t\t%lf", p, (Tot_load_served / (NBus - 2)));
	     particle[p].fit = (Tot_load_served / (NBus - 2));
	    //Freeing Memory
   	    for(i = 1; i <= (NBus+1); i++)
            {
      	        free(graph1[i]);
                graph1[i] == NULL;
            }
   	    free(graph1); graph1 == NULL;
	  }
        }
       if (ntie == 3)		//Three switch combintions
       {
	 fprintf(fpso, "\n\n**************************\nIteration = %d\n**************************\n\n", iter);
	 fprintf(fpso, "------------------------------------------\nParticle\t\t\tAverage Load Served\n------------------------------------------\n");
	 for (p = 1; p <= NPART; p++)	//NPART
	 {
	     //Memory Allocation
   	     graph1 = malloc((NBus+1) * sizeof(int *));
             if(graph1 == NULL)
             {
                printf("out of memory\n");
              exit(1);
           }
           for(i = 1; i < (NBus+1); i++)
           {
             graph1[i] = malloc((NBus+1) * sizeof(int));
             if(graph1[i] == NULL)
             {
 	        printf("out of memory\n");
 	        exit(1);
             }
           }

           for (i = 1; i <= NBus; i++)		//Copying the adjacency matrix to graph for further processing 
	   {
		for (j = 1; j <= NBus; j++)
		{
	           graph1[i][j] = adj[i][j];
		}
	   }
           threeswitch[0] = particle[p].p[1]; threeswitch[1] = particle[p].p[2];
           threeswitch[2] = particle[p].p[3]; threeswitch[3] = particle[p].p[4];
           threeswitch[4] = particle[p].p[5]; threeswitch[5] = particle[p].p[6];
        
      	  graph1[threeswitch[0]][threeswitch[1]] = 1;
	  graph1[threeswitch[1]][threeswitch[0]] = 1;
          graph1[threeswitch[2]][threeswitch[3]] = 1;
	  graph1[threeswitch[3]][threeswitch[2]] = 1;
          graph1[threeswitch[4]][threeswitch[5]] = 1;
	  graph1[threeswitch[5]][threeswitch[4]] = 1;


	  fitness_function();
	  fprintf(fpso, "\nparticle[%2d]\t\t\t%lf", p, (Tot_load_served / (NBus - 2)));
	  particle[p].fit = (Tot_load_served / (NBus - 2));
	  //Freeing Memory
   	  for(i = 1; i <= (NBus+1); i++)
          {
      	      free(graph1[i]);
              graph1[i] == NULL;
          }
   	  free(graph1); graph1 == NULL;
	}
      }

      if (ntie == 4)		//Four switch combintions
      {
	fprintf(fpso, "\n\n**************************\nIteration = %d\n**************************\n\n", iter);
	fprintf(fpso, "------------------------------------------\nParticle\t\t\tAverage Load Served\n------------------------------------------\n");
	for (p = 1; p <= NPART; p++)	//NPART
	{
	   //Memory Allocation
   	   graph1 = malloc((NBus+1) * sizeof(int *));
           if(graph1 == NULL)
           {
              printf("out of memory\n");
              exit(1);
           }
           for(i = 1; i < (NBus+1); i++)
           {
             graph1[i] = malloc((NBus+1) * sizeof(int));
             if(graph1[i] == NULL)
             {
 	        printf("out of memory\n");
 	        exit(1);
             }
           }

           for (i = 1; i <= NBus; i++)		//Copying the adjacency matrix to graph for further processing 
	   {
		for (j = 1; j <= NBus; j++)
		{
	           graph1[i][j] = adj[i][j];
		}
	   }
           fourswitch[0] = particle[p].p[1]; fourswitch[1] = particle[p].p[2];
           fourswitch[2] = particle[p].p[3]; fourswitch[3] = particle[p].p[4];
           fourswitch[4] = particle[p].p[5]; fourswitch[5] = particle[p].p[6];
           fourswitch[6] = particle[p].p[7]; fourswitch[7] = particle[p].p[8];
        
      	  graph1[fourswitch[0]][fourswitch[1]] = 1;
	  graph1[fourswitch[1]][fourswitch[0]] = 1;
          graph1[fourswitch[2]][fourswitch[3]] = 1;
	  graph1[fourswitch[3]][fourswitch[2]] = 1;
          graph1[fourswitch[4]][fourswitch[5]] = 1;
	  graph1[fourswitch[5]][fourswitch[4]] = 1;
          graph1[fourswitch[6]][fourswitch[7]] = 1;
	  graph1[fourswitch[7]][fourswitch[6]] = 1;


	  fitness_function();
	  fprintf(fpso, "\nparticle[%2d]\t\t\t%lf", p, (Tot_load_served / (NBus - 2)));
	  particle[p].fit = (Tot_load_served / (NBus - 2));
	  //Freeing Memory
   	  for(i = 1; i <= (NBus+1); i++)
          {
      	      free(graph1[i]);
              graph1[i] == NULL;
          }
   	  free(graph1); graph1 == NULL;
	}
      }

      if (ntie == 5)		//Five switch combintions
      {
	fprintf(fpso, "\n\n**************************\nIteration = %d\n**************************\n\n", iter);
	fprintf(fpso, "------------------------------------------\nParticle\t\t\tAverage Load Served\n------------------------------------------\n");
	for (p = 1; p <= NPART; p++)	//NPART
	{
	   //Memory Allocation
   	   graph1 = malloc((NBus+1) * sizeof(int *));
           if(graph1 == NULL)
           {
              printf("out of memory\n");
              exit(1);
           }
           for(i = 1; i < (NBus+1); i++)
           {
             graph1[i] = malloc((NBus+1) * sizeof(int));
             if(graph1[i] == NULL)
             {
 	        printf("out of memory\n");
 	        exit(1);
             }
           }

           for (i = 1; i <= NBus; i++)		//Copying the adjacency matrix to graph for further processing 
	   {
		for (j = 1; j <= NBus; j++)
		{
	           graph1[i][j] = adj[i][j];
		}
	   }
           fiveswitch[0] = particle[p].p[1]; fiveswitch[1] = particle[p].p[2];
           fiveswitch[2] = particle[p].p[3]; fiveswitch[3] = particle[p].p[4];
           fiveswitch[4] = particle[p].p[5]; fiveswitch[5] = particle[p].p[6];
           fiveswitch[6] = particle[p].p[7]; fiveswitch[7] = particle[p].p[8];
           fiveswitch[8] = particle[p].p[9]; fiveswitch[9] = particle[p].p[10];
        
      	  graph1[fiveswitch[0]][fiveswitch[1]] = 1;
	  graph1[fiveswitch[1]][fiveswitch[0]] = 1;
          graph1[fiveswitch[2]][fiveswitch[3]] = 1;
	  graph1[fiveswitch[3]][fiveswitch[2]] = 1;
          graph1[fiveswitch[4]][fiveswitch[5]] = 1;
	  graph1[fiveswitch[5]][fiveswitch[4]] = 1;
          graph1[fiveswitch[6]][fiveswitch[7]] = 1;
	  graph1[fiveswitch[7]][fiveswitch[6]] = 1;
          graph1[fiveswitch[8]][fiveswitch[9]] = 1;
	  graph1[fiveswitch[9]][fiveswitch[8]] = 1;


	  fitness_function();
	  fprintf(fpso, "\nparticle[%2d]\t\t\t%lf", p, (Tot_load_served / (NBus - 2)));
	  particle[p].fit = (Tot_load_served / (NBus - 2));
	  //Freeing Memory
   	  for(i = 1; i <= (NBus+1); i++)
          {
      	      free(graph1[i]);
              graph1[i] == NULL;
          }
   	  free(graph1); graph1 == NULL;
	}
      }
   
      //Determine the pbest and gbest for optimization
      tempa = 0.0; sum_fit = 0.0;  part_no = 0;		//Initialize part_no as zero
      tempa = particle[1].fit;

      //Finding best particle
      for (i = 1; i <= NPART; i++)	//NPART
      {
	  if (tempa < particle[i].fit)
	  {
	     tempa = particle[i].fit;
	     part_no = i;
	  }
      }

     pbest = particle[part_no].fit;
     if (iter == 1)
     {
	 new_gbest = gbest;
     }
     if (pbest > new_gbest)
	 new_gbest = pbest;

     fprintf(fpso, "\n\nBest particle number = %d\t\tPbest = %lf\t\tgbest = %lf", part_no, pbest, new_gbest);
     //fprintf(fpso, "\n\n%lf", new_gbest);

  /***************************************************************************
   *  			Step 3 : Update Velocity			     *
   ***************************************************************************/
   fprintf(fpso, "\n\n------------------------------------\n");
   fprintf(fpso, "Updated velocity");
   fprintf(fpso, "\n------------------------------------");
   w = wmax - (wmax - wmin)*iter / MITER;   //  w is dynamic weight factor
   for (p = 1; p <= NPART; p++)
   {
	xx = rand_number();
	yy = rand_number();
	particle[p].v_fin = (particle[p].v_int * w) + (c1*xx*(pbest - particle[p].fit)) + (c2*yy*(new_gbest - particle[p].fit));
	fprintf(fpso, "\n%lf", particle[p].v_fin);
   }

   //Assigning to local search space
   for (i = 1; i <= NPART; i++)	//NPART
   {
	particle[i].x = particle[i].fit;
   }

   //Update particle with new velocity
   fprintf(fpso, "\n\n------------------------------------\n");
   fprintf(fpso, "New particle in search space");
   fprintf(fpso, "\n------------------------------------");
   for (p = 1; p <= NPART; p++)
   {
	particle[p].x = particle[p].x + particle[p].v_fin;
	fprintf(fpso, "\n%lf", particle[p].x);
   }

  /***************************************************************************
   *  	Step 4 : Mapping updated particle to new particle		     *
   ***************************************************************************/
   //Generating new pair of particles based on updated velocity
   //fprintf(fpso, "\n\n**************************\nInside condition check\n**************************\n\n");
   for (a = 1; a <= NPART; a++)	//Random bus number and its velocity 
   {
   label4:
	//Random bus pair creation
	for (j = 1; j <= (2 * ntie); j++)
	{
            particle[a].updated_p[j] = random_bus();
	}

	//Condition to check whether both from and to buses are to be same.
	for (k = 1; k <= ((2 * ntie) - 1);)
	{
   	    if (particle[a].updated_p[k] == particle[a].updated_p[k + 1])
	    {
	     label3:
		particle[a].updated_p[k + 1] = random_bus();
		if (particle[a].updated_p[k + 1] == particle[a].updated_p[k])
		{
		   goto label3;
		}
	    }
	    k = k + 2;
	}
        //Printing purpose only
	//fprintf(fpso, "\n");
	for (j = 1; j <= (2 * ntie); j++)
	{
	    particle[a].p[j] = random_bus();
	    //fprintf(fpso, "%d\t", particle[a].p[j]);
	}
	//printf("\n");
        
        //Generated particle (random bus pairs) checking over it 
        if (ntie == 1)		//One switch combintions
        {
	  fprintf(fpso, "\n\n**************************\nIteration = %d\n**************************\n\n", iter);
	  fprintf(fpso, "------------------------------------------\nParticle\t\t\tAverage Load Served\n------------------------------------------\n");
	  for (p = 1; p <= NPART; p++)	//NPART
	  {
	     //Memory Allocation
   	     graph1 = malloc((NBus+1) * sizeof(int *));
             if(graph1 == NULL)
             {
                printf("out of memory\n");
                exit(1);
             }
             for(i = 1; i < (NBus+1); i++)
             {
               graph1[i] = malloc((NBus+1) * sizeof(int));
               if(graph1[i] == NULL)
               {
 	          printf("out of memory\n");
 	          exit(1);
               }
             }

             for (i = 1; i <= NBus; i++)		//Copying the adjacency matrix to graph for further processing 
	     {
		  for (j = 1; j <= NBus; j++)
		  {
	             graph1[i][j] = adj[i][j];
	  	  }
	     }
             oneswitch[0] = particle[p].p[1]; oneswitch[1] = particle[p].p[2];
        
      	    graph1[oneswitch[0]][oneswitch[1]] = 1;
	    graph1[oneswitch[1]][oneswitch[0]] = 1;

	    fitness_function();
	    fprintf(fpso, "\nparticle[%2d]\t\t\t%lf", p, (Tot_load_served / (NBus - 2)));
	    particle[p].fit = (Tot_load_served / (NBus - 2));
	    //Freeing Memory
   	    for(i = 1; i <= (NBus+1); i++)
            {
      	        free(graph1[i]);
                graph1[i] == NULL;
            }
   	    free(graph1); graph1 == NULL;
	  }
        }
       
        if (ntie == 2)		//Two switch combintions
        {
	  fprintf(fpso, "\n\n**************************\nIteration = %d\n**************************\n\n", iter);
	  fprintf(fpso, "------------------------------------------\nParticle\t\t\tAverage Load Served\n------------------------------------------\n");
	  for (p = 1; p <= NPART; p++)	//NPART
	  {
	     //Memory Allocation
   	     graph1 = malloc((NBus+1) * sizeof(int *));
             if(graph1 == NULL)
             {
                printf("out of memory\n");
                exit(1);
             }
             for(i = 1; i < (NBus+1); i++)
             {
               graph1[i] = malloc((NBus+1) * sizeof(int));
               if(graph1[i] == NULL)
               {
 	          printf("out of memory\n");
 	          exit(1);
               }
             }

             for (i = 1; i <= NBus; i++)		//Copying the adjacency matrix to graph for further processing 
	     {
		  for (j = 1; j <= NBus; j++)
		  {
	             graph1[i][j] = adj[i][j];
		  }
	     }
             twoswitch[0] = particle[p].p[1]; twoswitch[1] = particle[p].p[2];
             twoswitch[2] = particle[p].p[3]; twoswitch[3] = particle[p].p[4];
        
            graph1[twoswitch[0]][twoswitch[1]] = 1;
	    graph1[twoswitch[1]][twoswitch[0]] = 1;
            graph1[twoswitch[2]][twoswitch[3]] = 1;
	    graph1[twoswitch[3]][twoswitch[2]] = 1;

	    fitness_function();
	    fprintf(fpso, "\nparticle[%2d]\t\t\t%lf", p, (Tot_load_served / (NBus - 2)));
	    particle[p].fit = (Tot_load_served / (NBus - 2));
	    //Freeing Memory
   	    for(i = 1; i <= (NBus+1); i++)
            {
      	       free(graph1[i]);
               graph1[i] == NULL;
            }
   	    free(graph1); graph1 == NULL;
	  }
        }
        
        if (ntie == 3)		//Three switch combintions
        {
	  fprintf(fpso, "\n\n**************************\nIteration = %d\n**************************\n\n", iter);
	  fprintf(fpso, "------------------------------------------\nParticle\t\t\tAverage Load Served\n------------------------------------------\n");
	  for (p = 1; p <= NPART; p++)	//NPART
	  {
	     //Memory Allocation
   	     graph1 = malloc((NBus+1) * sizeof(int *));
             if(graph1 == NULL)
             {
                printf("out of memory\n");
                exit(1);
             }
             for(i = 1; i < (NBus+1); i++)
             {
               graph1[i] = malloc((NBus+1) * sizeof(int));
               if(graph1[i] == NULL)
               {
 	          printf("out of memory\n");
 	          exit(1);
               }
             }

             for (i = 1; i <= NBus; i++)		//Copying the adjacency matrix to graph for further processing 
	     {
		  for (j = 1; j <= NBus; j++)
		  {
	             graph1[i][j] = adj[i][j];
		  }
	     }
             threeswitch[0] = particle[p].p[1]; threeswitch[1] = particle[p].p[2];
             threeswitch[2] = particle[p].p[3]; threeswitch[3] = particle[p].p[4];
             threeswitch[4] = particle[p].p[5]; threeswitch[5] = particle[p].p[6];
        
      	    graph1[threeswitch[0]][threeswitch[1]] = 1;
	    graph1[threeswitch[1]][threeswitch[0]] = 1;
            graph1[threeswitch[2]][threeswitch[3]] = 1;
	    graph1[threeswitch[3]][threeswitch[2]] = 1;
            graph1[threeswitch[4]][threeswitch[5]] = 1;
	    graph1[threeswitch[5]][threeswitch[4]] = 1;


	    fitness_function();
	    fprintf(fpso, "\nparticle[%2d]\t\t\t%lf", p, (Tot_load_served / (NBus - 2)));
	    particle[p].fit = (Tot_load_served / (NBus - 2));
	    //Freeing Memory
   	    for(i = 1; i <= (NBus+1); i++)
            {
      	        free(graph1[i]);
                graph1[i] == NULL;
            }
   	    free(graph1); graph1 == NULL;
	  }
        }

        if (ntie == 4)		//Four switch combintions
        {
	  fprintf(fpso, "\n\n**************************\nIteration = %d\n**************************\n\n", iter);
	  fprintf(fpso, "------------------------------------------\nParticle\t\t\tAverage Load Served\n------------------------------------------\n");
	  for (p = 1; p <= NPART; p++)	//NPART
	  {
	     //Memory Allocation
   	     graph1 = malloc((NBus+1) * sizeof(int *));
             if(graph1 == NULL)
             {
                printf("out of memory\n");
                exit(1);
             }
             for(i = 1; i < (NBus+1); i++)
             {
               graph1[i] = malloc((NBus+1) * sizeof(int));
               if(graph1[i] == NULL)
               {
 	          printf("out of memory\n");
 	          exit(1);
               }
             }

             for (i = 1; i <= NBus; i++)		//Copying the adjacency matrix to graph for further processing 
	     {
		  for (j = 1; j <= NBus; j++)
		  {
	             graph1[i][j] = adj[i][j];
		  }
	     }
             fourswitch[0] = particle[p].p[1]; fourswitch[1] = particle[p].p[2];
             fourswitch[2] = particle[p].p[3]; fourswitch[3] = particle[p].p[4];
             fourswitch[4] = particle[p].p[5]; fourswitch[5] = particle[p].p[6];
             fourswitch[6] = particle[p].p[7]; fourswitch[7] = particle[p].p[8];
        
      	    graph1[fourswitch[0]][fourswitch[1]] = 1;
	    graph1[fourswitch[1]][fourswitch[0]] = 1;
            graph1[fourswitch[2]][fourswitch[3]] = 1;
	    graph1[fourswitch[3]][fourswitch[2]] = 1;
            graph1[fourswitch[4]][fourswitch[5]] = 1;
	    graph1[fourswitch[5]][fourswitch[4]] = 1;
            graph1[fourswitch[6]][fourswitch[7]] = 1;
	    graph1[fourswitch[7]][fourswitch[6]] = 1;


	    fitness_function();
	    fprintf(fpso, "\nparticle[%2d]\t\t\t%lf", p, (Tot_load_served / (NBus - 2)));
	    particle[p].fit = (Tot_load_served / (NBus - 2));
	    //Freeing Memory
   	    for(i = 1; i <= (NBus+1); i++)
            {
      	        free(graph1[i]);
                graph1[i] == NULL;
            }
   	    free(graph1); graph1 == NULL;
	  }
        }
        
        if (ntie == 5)		//Five switch combintions
        {
	  fprintf(fpso, "\n\n**************************\nIteration = %d\n**************************\n\n", iter);
	  fprintf(fpso, "------------------------------------------\nParticle\t\t\tAverage Load Served\n------------------------------------------\n");
	  for (p = 1; p <= NPART; p++)	//NPART
	  {
	     //Memory Allocation
   	     graph1 = malloc((NBus+1) * sizeof(int *));
             if(graph1 == NULL)
             {
                printf("out of memory\n");
                exit(1);
             }
             for(i = 1; i < (NBus+1); i++)
             {
               graph1[i] = malloc((NBus+1) * sizeof(int));
               if(graph1[i] == NULL)
               {
 	          printf("out of memory\n");
 	          exit(1);
               }
             }

             for (i = 1; i <= NBus; i++)		//Copying the adjacency matrix to graph for further processing 
	     {
		  for (j = 1; j <= NBus; j++)
		  {
	             graph1[i][j] = adj[i][j];
		  }
	     }
             fiveswitch[0] = particle[p].p[1]; fiveswitch[1] = particle[p].p[2];
             fiveswitch[2] = particle[p].p[3]; fiveswitch[3] = particle[p].p[4];
             fiveswitch[4] = particle[p].p[5]; fiveswitch[5] = particle[p].p[6];
           fiveswitch[6] = particle[p].p[7]; fiveswitch[7] = particle[p].p[8];
           fiveswitch[8] = particle[p].p[9]; fiveswitch[9] = particle[p].p[10];
        
      	  graph1[fiveswitch[0]][fiveswitch[1]] = 1;
	  graph1[fiveswitch[1]][fiveswitch[0]] = 1;
          graph1[fiveswitch[2]][fiveswitch[3]] = 1;
	  graph1[fiveswitch[3]][fiveswitch[2]] = 1;
          graph1[fiveswitch[4]][fiveswitch[5]] = 1;
	  graph1[fiveswitch[5]][fiveswitch[4]] = 1;
          graph1[fiveswitch[6]][fiveswitch[7]] = 1;
	  graph1[fiveswitch[7]][fiveswitch[6]] = 1;
          graph1[fiveswitch[8]][fiveswitch[9]] = 1;
	  graph1[fiveswitch[9]][fiveswitch[8]] = 1;


	  fitness_function();
	  fprintf(fpso, "\nparticle[%2d]\t\t\t%lf", p, (Tot_load_served / (NBus - 2)));
	  particle[p].fit = (Tot_load_served / (NBus - 2));
	  //Freeing Memory
   	  for(i = 1; i <= (NBus+1); i++)
          {
      	      free(graph1[i]);
              graph1[i] == NULL;
          }
   	  free(graph1); graph1 == NULL;
	}
      }
      
      t1 = 0.0; t2 = 0.0;
      t1 = particle[a].fit;
      t2 = particle[a].fit + (particle[a].fit *0.1);
      
      //fprintf(fpso, "\n\t\t\t\t\t\t\t%lf\t\t%lf", t1, particle[a].up_fit);
      if (particle[a].fit > particle[a].up_fit)
	  goto label4;
      else if (particle[a].fit == particle[a].up_fit)
      {
	particle[a].p = particle[a].updated_p;
      }
      else
      {
	  particle[a].p = particle[a].updated_p;
      }

      //Termination condition have to include
       
     	      
    }
    //printf("\nIter = %d",iter);
    iter = iter+1;
   }while(iter <= MITER); 	      

   
   
   
   	     
   fclose(fpso);
   //Pso function final memory clearing stage
   for (i = 1; i <= NPART; i++)
   {
      free(particle[i].p);
      free(particle[i].updated_p);
   }

}

double rand_number()
{
	int x1;
	double no1;
	x1 = rand() % (Hi - Low + 1) + 0;		//high = 100, low =0
	no1 = x1*.01;
	return no1;
}

int random_bus()
{
	int x2;
label1:
	x2 = rand() % (NBus+1);
	if (x2 <= NBus && x2 > 1)
		return x2;
	else
	{
		goto label1;
	}
}

void Risk_factor_calculation()
{
   int **graph,*Is1_Node, *Is2_Node,*x, *y, *temp, *temp1,*LAF;
   int node_is1, node_is2,No_Island,t1, t2, t3, temp2,add = 0;
   double ls,temp_rf,I_gen, I_tot_gen,I_load, I_tot_load;
   node_is1 = 0; node_is2 = 0; 	No_Island = 0;
  
   //Memory Allocation
   graph = malloc((NBus+1) * sizeof(int *));
   if(graph == NULL)
   {
     printf("out of memory\n");
     exit(1);
   }
   for(i = 1; i < (NBus+1); i++)
   {
     graph[i] = malloc((NBus+1) * sizeof(int));
     if(graph[i] == NULL)
     {
 	printf("out of memory\n");
 	exit(1);
     }
   }

   for (i = 1; i <= NBus; i++)		
   {
     for (j = 1; j <= NBus; j++)
       {
	 graph[i][j] = adj[i][j];
       }
   }
  
   graph[ref_mat[0]][ref_mat[1]] = 0;	//graph
   graph[ref_mat[1]][ref_mat[0]] = 0;

  /************************************************************
   *   Determination of nodes in each every island            * 
   ************************************************************/

   x = (int *)malloc((NBus+1) * sizeof(int));
   y = (int *)malloc((NBus+1) * sizeof(int));
   temp = (int *)malloc((NBus+1) * sizeof(int));
   temp1 = (int *)malloc((NBus+1) * sizeof(int));
   Is1_Node = (int *)malloc((NBus+1) * sizeof(int));
   Is2_Node = (int *)malloc((NBus+1) * sizeof(int));
   LAF = (int *)malloc((NBus+1) * sizeof(int));
     
   x[1] = 1;
   for (i = 1; i <= NBus; i++)
   {
      y[i] = x[i];
      //printf("\n%d", y[i]);
   }

 doagain:
   for (i = 1; i <= NBus; i++)
   {
      y[i] = x[i];
   }
   
  /**********************************
   *	Matrix Multiplication	     *	
   **********************************/
   for (i = 1; i <= NBus; i++)
   {
      for (j = 1; j <= NBus; j++)
	{
            add = graph[i][j] * x[j] + add;
	}
      temp[i] = add;
      add = 0;
   }
  
  /**********************************
   *	Matrix Addition 	    *	
   **********************************/
   for (i = 1; i <= NBus; i++)
   {
	add = 0;
	temp1[i] = temp[i] + x[i];
	if (temp1[i] != 0)
	{
		temp1[i] = temp1[i] / temp1[i];
		x[i] = temp1[i];
	}
	else
	{
		x[i] = 0;
	}
   }
   
   t1 = 0; t2 = 0; t3 = 0;
   for (i = 1; i <= NBus; i++)
   {
	if (x[i] == y[i] && x[i] == 1)
	{
		t1 = t1++;
	}

	if (x[i] == y[i] && x[i] == 0)
	{
		t2 = t2++;
	}
	if (x[i] != y[i])
		t3 = t3++;
   }

   if (t1 == NBus)
   {
	//fprintf(test, "\nAll nodes are connected\n");
	for (j = 1; j <= NBus; j++)
	{
		LAF[j] = 1;
	}
	No_Island = 1;
   }
   if ((t1 + t2) != NBus)
   {
	goto doagain;
   }
   if ((t2 != 0) && ((t1 + t2) == NBus))
   {
	node_is1 = t1;
	node_is2 = t2;
	for (i = 1, k = 1, l = 1; i <= NBus; i++)
 	{
		if (x[i] == 1)
		{
			Is1_Node[k] = i;
			k++;
		}
		else
		{
			Is2_Node[l] = i;
			l++;
		}
	}
	//Printing it in to a testing file 
	fprintf(test, "\n\nNode in Island-1 are as follows\n\n");
	for (i = 1; i <= node_is1; i++)
	    fprintf(test, "%d\t", Is1_Node[i]);
	fprintf(test, "\n\nNode in Island-2 are as follows\n\n");
	for (i = 1; i <= node_is2; i++)
	    fprintf(test, "%d\t", Is2_Node[i]);
   }
   
   /************************************************************
   *   		Formation of Load Affordability Factor         * 
   *************************************************************/
   temp2 = 0; I_tot_gen = 0; I_gen = 0; I_load = 0; I_tot_load = 0;
   
   //LAF matrix for island 1 
   for (j = 1; j <= node_is1; j++)		//Making the nodes in island are to be 1.
   {
	temp2 = Is1_Node[j];
	LAF[temp2] = 1;
   }
   
   if (node_is2 != 0)		//Island formed with any number of nodes
   {
	for (j = 1; j <= node_is2; j++)		//going to check generation and load at any second island alone.
	{
	    temp2 = Is2_Node[j];
	    I_gen = generation[temp2];
	    I_tot_gen = I_gen + I_tot_gen;
	    I_load = bd[temp2].Pl;
	    I_tot_load = I_tot_load + I_load;
	    I_gen = 0.0; I_load = 0.0;
	}
	//printf("\nI_tot_gen = %lf\tI_tot_load = %lf", I_tot_gen, I_tot_load);

	if ((I_tot_gen - I_tot_load) < 0)
	{
	    for (j = 1; j <= node_is2; j++)		//going to check generation and load at any second island alone.
	    {
		temp2 = Is2_Node[j];
		LAF[temp2] = 0;
	    }
	}
	else
	{
	    for (j = 1; j <= node_is2; j++)	 //Load at the island can be able to surve the load	
	    {
		temp2 = Is2_Node[j];
		LAF[temp2] = 1;
	    }
	}
   }
   
   //Prinitng it into a testing file 
   /*fprintf(test, "\n");
   for (i = 1; i <= NBus; i++)
   {
	fprintf(test, "%d\t", LAF[i]);
   }*/
   
   temp_rf = 0; ls = 0; load_served = 0;
   
   //Risk factor calculation 
   for (i = 1; i <= NBus; i++)	
   {
	temp_rf = load_weight[i] * LAF[i];
	risk_factor = temp_rf + risk_factor;
	temp_rf = 0;
	ls = (bd[i].Pl * LAF[i]) + ls;
   }
   risk_factor = 1 - risk_factor;
   load_served = ls;
   //printf("\nRisk factor = %lf\tload_served = %lf",risk_factor,load_served);
   //fprintf(test,"\nRisk factor = %lf", risk_factor);

   //Freeing Memory
   for(i = 1; i <= (NBus+1); i++)
   {
      free(graph[i]);
      graph[i] == NULL;
   }
   free(graph); graph == NULL;
   free(x); free(y); free(temp); free(temp1);
   free(Is1_Node); free(Is2_Node); free(LAF);
}


double fitness_function()
{
  int nc;
  sum_risk_factor = 0.0;
  Tot_load_served = 0.0;
  //fprintf(test, "\nIteration = %d\n", k);
  for (nc = 2; nc <= (NBus-1); nc++)	//considering the grid connected faults alone
  {
     ref_mat[0] = fd[nc].fb;		//From bus and to bus to be opened
     ref_mat[1] = fd[nc].tb;
     Risk_factor_calculation();
     sum_risk_factor = sum_risk_factor + risk_factor;
     Tot_load_served = Tot_load_served + load_served;
     risk_factor = 0;
     /*if (nc==(NBus-1))
     {
       	printf("\n\n-------------------------------------------------------------------------------\n");
	printf("Sum of Risk_factor = %lf\tAverage_load_served = %lf", sum_risk_factor / (NBus - 2), Tot_load_served/ (NBus - 2));
	printf("\n-------------------------------------------------------------------------------\n");
     }*/
  }
  return(Tot_load_served);
}


void read_generation_data()
{
  int aa;
  fgen = fopen("Input/gen_data69.pdt", "r");
  if (fgen == NULL)
  {
    printf("\nError in opening generator data file\nPlease check");
    exit(1);
  }
  fscanf(fgen, "%d", &No_gen);

  gen_node = (int*)malloc((10) *sizeof(int));
  generation = (double*)malloc((NBus+1) * sizeof(double));
  gen_node_value = (double*)malloc((NBus+1) * sizeof(double));
  if(generation == NULL || gen_node == NULL || gen_node_value == NULL)
    printf("\nError - read_generation_data function\n");
  for (aa = 1; aa <= No_gen; aa++)
  {
    fscanf(fgen, "%d %lf", &gen_node[aa], &gen_node_value[aa]);
    generation[gen_node[aa]] = gen_node_value[aa];
  }
  //printing the genration value and node at which it present.
  /*printf("\n\nDGs value and its connected bus\n");
  for (aa = 1; aa <= NBus; aa++)
  {
    printf("Bus[%d] = %lf\n", aa, generation[aa]);
  }*/
  fclose(fgen);
  fgen == NULL;
}

void Adjacency_matrix_formation()
{
  int  j, l,i,cnt;
  double inv_cnt;			//inverse of cnt variable;

  //Memory allocation
  fd = (Feederdata*)malloc((NBus+1) * sizeof(Feederdata)); 
  bd = (Busdata*)malloc((NBus+1) * sizeof(Busdata));	
  if(fd == NULL)
     printf("\nError in memory Allocation\n\n");

  
  //Reading pdt file 
  fpdt = fopen("Input/PdTemp69.pdt", "r");
  if (fpdt == NULL)
  {
   printf("The file PdTemp.pdt was not opened\n");
   exit(1);
  }
		
 else
 {
  for (i = 1; i <= (NBus-1); i++)
  {
   fscanf(fpdt, "%d %d %d", &fd[i].fb, &fd[i].tb, &fd[i].status);
   //printf("\n%d %d %d", fd[i].fb, fd[i].tb, fd[i].status);
  }
  for (i = 1; i <= NBus; i++)
  {
   fscanf(fpdt, "%d %lf %lf %lf %lf", &bd[i].bus_no, &bd[i].Pg, &bd[i].Qg, &bd[i].Pl, &bd[i].Ql);
   //printf("\n%d\t%lf\t%lf\t%lf\t%lf", bd[i].bus_no, bd[i].Pg, bd[i].Qg, bd[i].Pl, bd[i].Ql);
  }
 }
 fclose(fpdt);
 fpdt == NULL;
 //user preferences
 printf("\n-------------------------\nload weightage\n-------------------------\n");
 printf("\nPress 0 for Equal load weightage\nPress 1 for unequal load weightage:\t");
 scanf("%d", &lweight_choice);

 load_weight = (double *)malloc((NBus+1) * sizeof(double)); 
 if (lweight_choice == 0)	//Equal load weightage
 {
   cnt = 0; //cnt - count variable
   inv_cnt = 0;
   for (i = 1; i <= NBus; i++)
   {
    if (bd[i].Pl > 0)
      {
	cnt++;
      }	
   }
   inv_cnt = 0.02083;	//(1/14) for 15 bus system (1/48) for 69 bus 
  
   for (i = 1; i <= NBus; i++)
  {
   if (bd[i].Pl > 0)
   {
    load_weight[i] = inv_cnt;
    //printf("\nLoad_weighhtage[%d] = %lf", i, load_weight[i]);
   }
   else
   {
    load_weight[i] = 0;
    //printf("\nLoad_weighhtage[%d] = %lf", i, load_weight[i]);
   }
  }
 }
 else						//Unequal load weightage
 {
   fwei = fopen("Input/Weight69.pdt", "r");
   if (fwei == NULL)
   {
     printf("The file weightage.pdt was not opened\n");
     exit(1);
   }
   else
   {
     //printf("\nLoad Weightage are as follows\n");
     for (i = 1; i <= NBus; i++)
     {
	fscanf(fwei, "%lf", &load_weight[i]);
	//printf("\nload[%d] = %lf",i,load_weight[i]);
     }
     fclose(fwei);
     fwei == NULL;
   }
  }

  adj = malloc((NBus+1) * sizeof(int *));
  if(adj == NULL)
  {
    printf("out of memory\n");
    exit(1);
  }
  for(i = 1; i <= (NBus+1); i++)
  {
    adj[i] = malloc((NBus+1) * sizeof(int));
    if(adj[i] == NULL)
    {
	printf("out of memory\n");
	exit(1);
    }
  }

  //Creating the adjacency matrix for the given network
  for (i = 1; i <= NBus; i++)
  {
    for (j = 1; j <= NBus; j++)
    {
      if (i == j)
	adj[i][j] = 0;		//Making the diagonal element to be 0
      else			//Upper and Lower diagonal element
      {
	for (l = 1; l<=(NBus-1); l++)
	{
	  if (i == fd[l].fb && j == fd[l].tb && fd[l].status == 1)
	  {
	    adj[i][j] = 1;
	  }
	  if (i == fd[l].tb && j == fd[l].fb && fd[l].status == 1)
	  {
	    adj[i][j] = 1;
	  }
	}
      }
    }
  }

  //Printing the adjacency matrix into a output file
  fadm = fopen("Output/Adjacency.pdt", "w");
  if (fadm == NULL)
  {
    printf("\nError in opening a new file ");
  }
  fprintf(fadm, "\n\n----------------------------------------------\nAdjacency Matrix is as follows");
  fprintf(fadm,"\n----------------------------------------------\n");
  for (i = 1; i <= NBus; i++)
  {
    for (j = 1; j <= NBus; j++)
    {
	fprintf(fadm, "%d\t", adj[i][j]);
    }
    fprintf(fadm, "\n");
  }
  fclose(fadm);
  fadm == NULL;
}


