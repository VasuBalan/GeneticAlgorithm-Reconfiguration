/************************************************************************************************
*				particle swarm optimization														*
*written by Vasudevan.B (vasubdevan@yahoo.com)													*
*Last Modified on 27.05.15																		* 
/************************************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#define NBus 15
#define NPART 10
#define Hi 100
#define Low 0
#define ntie 1
#define MITER 1

//Structure pointer
typedef struct {
	int fb;
	int tb;
	int status;
}Feederdata;

typedef struct {
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
void fitness_function(int **graph1);
void particle_swarm();
int random_bus();
double rand_number();


//Pointer declaration
Feederdata *fd; Busdata *bd;
double *load_weight,*generation,*gen_node_value;
int **adj,*gen_node;


//File Pointer
FILE *fpdt,*fgen,*fadm,*fwei,*test;


//Variable declaration
int lweight_choice,No_gen,i,j,k,l,ref_mat[3],iter;
double net_obj_value,Avg_ls;

clock_t t; //clock time 


int main(int argc,char *argv[]) {
   int i,j,**graph1;
   printf("\n********************************************************************************\n");
   printf("\t\tProgram to determine optimal Tie switch using PSO\n");
   printf("********************************************************************************\n");
   Adjacency_matrix_formation();  
   read_generation_data();
   
   
   graph1 = (int **)malloc((NBus+1) * sizeof(int*));
   if(adj == NULL)
   {
      printf("out of memory\n");
      exit(1);
   }
   for(i = 1; i <= (NBus+1); i++)
   {
      graph1[i] = (int*)malloc((NBus+1) * sizeof(int));
      if(graph1[i] == NULL)
      {
	  printf("out of memory\n");
	  exit(1);
      }
   }

   //Initialize to zero the copy adjacency matrix 
   for (i = 1; i <= NBus; i++) {
     for (j = 1; j <= NBus; j++) {
	 graph1[i][j] = 0;
         graph1[i][j] = adj[i][j];
     }
   }

   test = fopen("Output/Testing.pdt", "w");
   if (test == NULL) {
      printf("\nError in creating output file\n");
   }

   fitness_function(graph1);

   for(i = 1; i <= (NBus+1); i++){
       free(graph1[i]);
   }
   free(graph1); 
   
   t = clock();		//Time calculation starts
   particle_swarm();
   
   t = clock() - t;
   //printf ("It took me %d clicks (%lf seconds).\n",t,((double)t)/CLOCKS_PER_SEC);
   
   fclose(test);

   
   //Freeing memory
   free(fd); free(bd);   
   for(i = 1; i <= (NBus+1); i++)
   {
       free(adj[i]);
   }
   free(adj); free(load_weight); free(generation);
   free(gen_node); free(gen_node_value);
}

void particle_swarm()
{
   abc particle[NPART + 1];
   FILE *fpso;
   int p,**graph1;
   int oneswitch[3];
  
  
  
   for (i = 1; i <= NPART; i++)
   {
      particle[i].p = (int*)calloc(NBus+1,sizeof(int));
      particle[i].updated_p = (int*)calloc(NBus+1, sizeof(int));
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
		particle[i].p[j] = random_bus();
		fprintf(fpso, "%d\t", particle[i].p[j]);
	}
	//Initial Velocity
	particle[i].v_int = rand_number();
	fprintf(fpso, "\t%lf", particle[i].v_int);
	fprintf(fpso, "\n");
    }// step-1 ends
    
   /***************************************************************************
    *  		Step 2 : Run objective fucntion and find pbest		     *
    ***************************************************************************/
   iter = 1;
   do
   {
     if(ntie == 1)
     {
       fprintf(fpso, "\n\n**************************\nIteration = %d\n**************************\n\n", iter);
       fprintf(fpso, "------------------------------------------\nParticle\t\t\tAverage Load Served\n------------------------------------------\n");
       for (p = 1; p <= NPART; p++)	//NPART
       {
	 graph1 = (int **)malloc((NBus+1) * sizeof(int*));
         if(adj == NULL)
         {
             printf("out of memory\n");
             exit(1);
         }
         for(i = 1; i <= (NBus+1); i++)
         {
             graph1[i] = (int*)malloc((NBus+1) * sizeof(int));
             if(graph1[i] == NULL)
             {
		printf("out of memory\n");
		exit(1);
	     }
	 }
        
        //Initialize to zero the copy adjacency matrix 
	for (i = 1; i <= NBus; i++) {
	for (j = 1; j <= NBus; j++) {
	    graph1[i][j] = 0;
	    graph1[i][j] = adj[i][j];
	} }
	
	oneswitch[0] = particle[p].p[1]; oneswitch[1] = particle[p].p[2];
	printf("\n\tparticle = %d",p);
	
	graph1[oneswitch[0]][oneswitch[1]] = 1;
	graph1[oneswitch[1]][oneswitch[0]] = 1;    
	
	fitness_function(graph1);
	fprintf(fpso,"\nAverage load served = %lf",Avg_ls);
	//clearing inside for loop
	for(i = 1; i <= (NBus+1); i++){
           free(graph1[i]);
        }
        free(graph1);  
	 
	 
       }//part loop ends
       
     }//ntie = 1 ends
     
     
     
     
     
     printf("\nIter = %d",iter);
     iter = iter +1;
   }while(iter <=MITER);
   
  
    
    
    
   
  
  
  
  fclose(fpso);
  //Final Freeing memory
  for (i = 1; i <= NPART; i++)
  {
    free(particle[i].p);
    free(particle[i].updated_p);
  }
      
  printf("\nWelcome to swarm intelligence\n");
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

void fitness_function(int **graph1) {
   int nc,add,t1,t2,t3,No_Island,node_is2,node_is1,temp2,overall_ls;
   int **graph,*Is1_Node, *Is2_Node,*x, *y, *temp, *temp1,*LAF;
   double I_gen, I_tot_gen,I_load,I_tot_load,num_term;
   double *prod_mat,*num_obj,sum_nt,*ind_bus_ls,sum_ind_bus_ls,*tot_ls;
   
   
   num_obj = (double*)calloc((NBus+1), sizeof(double));
   tot_ls = (double*)calloc((NBus+1), sizeof(double));
   for (nc = 2; nc <= (NBus-1); nc++)	//considering the grid connected faults alone
   {
     ref_mat[0] = fd[nc].fb;		//From bus and to bus to be opened
     ref_mat[1] = fd[nc].tb;
     
     graph = (int **)malloc((NBus+1) * sizeof(int*));
     if(graph == NULL) {
	printf("out of memory\n");
	exit(1);
     }
     for(i = 1; i <= (NBus+1); i++) 
     {
       	graph[i] = (int*)malloc((NBus+1) * sizeof(int));
	if(graph[i] == NULL) 
	{
	  printf("out of memory\n");
	  exit(1); 
	 }
     }
     
     //Initialize to zero - malloc 
     for (i = 1; i <= NBus; i++) {
     for (j = 1; j <= NBus; j++) {
	graph[i][j] = 0; } }
     
     for (i = 1; i <= NBus; i++)		
     {
       for (j = 1; j <= NBus; j++)
       {
	 graph[i][j] = graph1[i][j];
       }
     }
  
     graph[ref_mat[0]][ref_mat[1]] = 0;	//graph
     graph[ref_mat[1]][ref_mat[0]] = 0;
     //printf("graph[%d][%d]\n",ref_mat[0],ref_mat[1]);

     /************************************************************
     *   Determination of nodes in each every island            * 
     ************************************************************/
     x = (int *)calloc((NBus+1),sizeof(int));
     y = (int *)calloc((NBus+1),sizeof(int));
     temp = (int *)calloc((NBus+1),sizeof(int));
     temp1 = (int *)calloc((NBus+1),sizeof(int));
     Is1_Node = (int *)calloc((NBus+1),sizeof(int));
     Is2_Node = (int *)calloc((NBus+1), sizeof(int));
     LAF = (int *)calloc((NBus+1), sizeof(int));
     prod_mat = (double*)calloc((NBus+1), sizeof(double));
     ind_bus_ls = (double*)calloc((NBus+1), sizeof(double));
     
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
      *	 Matrix Multiplication	        *
      **********************************/
     add = 0;
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
     *	Matrix Addition 	        *
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
	/*fprintf(test, "\n\nNode in Island-1 are as follows\n\n");
	for (i = 1; i <= node_is1; i++)
	   fprintf(test, "%d\t", Is1_Node[i]);
	fprintf(test, "\n\nNode in Island-2 are as follows\n\n");
	for (i = 1; i <= node_is2; i++)
	   fprintf(test, "%d\t", Is2_Node[i]);*/
     }
     
     /************************************************************
      *         Formation of Load Affordability Matrix           * 
      ************************************************************/
      
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
            I_gen = 0.0;
            I_load = 0.0;
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
      
      //Determination of risk factor
      //printf("\n\n");
      for (i = 1; i <= NBus; i++)	
      {
          prod_mat[i] = load_weight[i] * LAF[i];
          ind_bus_ls[i] = LAF[i] * bd[i].Pl; 
          //printf("\nprod_mat[%d]= %lf",i,prod_mat[i]);
          //printf("\nind_bus_ls[%d]= %lf",i,ind_bus_ls[i]);  
      }
      
      //Individual line opening risk factor
      num_term = 0; sum_ind_bus_ls = 0;
      for (i = 1; i <= NBus; i++)	
      {
	  num_term = num_term + prod_mat[i];
	  sum_ind_bus_ls = sum_ind_bus_ls + ind_bus_ls[i];
      }
      //fprintf(test,"\nRisk factor corresponding to line opening between bus %d  and bus %d  is %lf",ref_mat[0],ref_mat[1],num_term);
      //fprintf(test,"\nTotal load served for a single contingency is %lf",sum_ind_bus_ls);
      //printf("\nIndividual Risk factor = %lf",num_term);
      //printf("\nTotal load served  = %lf",sum_ind_bus_ls);
      
      num_obj[nc] = num_term;
      tot_ls[nc] = sum_ind_bus_ls;
      
      //Memory clearing indide fitness function 
     for(i = 1; i <= (NBus+1); i++)
     {
        free(graph[i]);
     }
     free(graph); 
     free(x); free(y); free(temp); free(temp1);
     free(Is1_Node); free(Is2_Node); free(LAF);
     free(prod_mat); free(ind_bus_ls);
   }
   
   //Calculating objective function value -Risk factor alone
   sum_nt = 0; net_obj_value = 0; overall_ls = 0; Avg_ls = 0;
   
   for (nc = 2; nc <= (NBus-1); nc++)  // 1 - indi risk factor 
   { 
       num_obj[nc] = 1 - num_obj[nc];
   }
   
   for (nc = 2; nc <= (NBus-1); nc++)	
   { 
     //printf("\nnum_obj[%d] = %lf",nc,num_obj[nc]);
     sum_nt = sum_nt + num_obj[nc];
     overall_ls = overall_ls + tot_ls[nc];
   }
   
   //printf("\n\nSum of objective function value is %lf",sum_nt);
   net_obj_value = (sum_nt/(NBus-2));
   Avg_ls = (overall_ls/(NBus-2));
   
   fprintf(test,"\nNet objective function value = %lf\n",net_obj_value);
   fprintf(test,"\nAverage load served = %lf\n\n",Avg_ls);
   
   //Freeing memory outside loop (line opening)
   free(num_obj);free(tot_ls);
}

void read_generation_data(){
  int aa;
  fgen = fopen("Input/gen_data15.pdt", "r");
  if (fgen == NULL)
  {
    printf("\nError in opening generator data file\nPlease check");
    exit(1);
  }
  fscanf(fgen, "%d", &No_gen);

  gen_node = (int*)calloc((10),sizeof(int));
  generation = (double*)calloc((NBus+1),sizeof(double));
  gen_node_value = (double*)calloc((NBus+1),sizeof(double));
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
}

void Adjacency_matrix_formation(){
  int  j, l,i;
  double cnt,inv_cnt;			//inverse of cnt variable;
  double Tot_load_given;

  //Memory allocation
  fd = (Feederdata*)calloc((NBus+1), sizeof(Feederdata)); 
  bd = (Busdata*)calloc((NBus+1), sizeof(Busdata));	
  if(fd == NULL)
     printf("\nError in memory Allocation\n\n");

  
  //Reading pdt file 
  fpdt = fopen("Input/PdTemp15.pdt", "r");
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
 //user preferences
 printf("\n----------------------------\nPrefernce for Load Weightage\n----------------------------\n");
 printf("\nPress 0 for Equal load weightage\nPress 1 for unequal load weightage:\t");
 scanf("%d", &lweight_choice);

 load_weight = (double *)malloc((NBus+1) * sizeof(double)); 
 if (lweight_choice == 0)	//Equal load weightage
 {
   cnt = 0; //cnt - count variable
   inv_cnt = 0; Tot_load_given = 0;
   for (i = 1; i <= NBus; i++)
   {
    if (bd[i].Pl > 0)
      {
	cnt++;
	Tot_load_given = Tot_load_given + bd[i].Pl;
      }	
   }
   printf("\nThe number of buses contain load are %2.0lf and total load = %lf\n\n",cnt,Tot_load_given); 
   inv_cnt = (1/cnt);	
   //printf("\nEqual load weighhtage = %lf",inv_cnt);
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
   cnt = 0;
   fwei = fopen("Input/Weight15.pdt", "r");
   if (fwei == NULL)
   {
     printf("The file weightage.pdt was not opened\n");
     exit(1);
   }
   else
   {
     for (i = 1; i <= NBus; i++)
     {
      if (bd[i].Pl > 0)
      {
        cnt++;
        Tot_load_given = Tot_load_given + bd[i].Pl;
      }	
     }
     printf("\nThe number of buses contain load are %2.0lf and total load = %lf\n\n",cnt,Tot_load_given); 
     //printf("\nLoad Weightage are as follows\n");
     for (i = 1; i <= NBus; i++)
     {
	fscanf(fwei, "%lf", &load_weight[i]);
	//printf("\nload[%d] = %lf",i,load_weight[i]);
     }
     fclose(fwei);
   }
  }

  adj = (int **)malloc((NBus+1) * sizeof(int*));
  if(adj == NULL)
  {
    printf("out of memory\n");
    exit(1);
  }
  for(i = 1; i <= (NBus+1); i++)
  {
    adj[i] = (int*)malloc((NBus+1) * sizeof(int));
    if(adj[i] == NULL)
    {
	printf("out of memory\n");
	exit(1);
    }
  }

  //Initialize to zero - malloc 
  for (i = 1; i <= NBus; i++)
  {
    for (j = 1; j <= NBus; j++)
    {
	adj[i][j] = 0;
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
}



