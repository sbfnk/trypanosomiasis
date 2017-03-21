//---PREPROCESSOR COMMANDS----------------------------------
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <ctime>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
using namespace std;

int no_vars=22;
int var1=0;
int scores=0;
int score=0;
double lambda=0;
int rowmatrix=0;
int colmatrix=0;
gsl_rng *r; //sets a pointer for the RNG


//---FUNCTIONS-----------------------------------------------
//-------------------------------------------------------------

//---TRIANGULAR DISTRIBUTION---------------------------------
double rand_triang_dist(double min,double mode,double max)
{
  //Generate quasi random value from a triangular distribution of min mode max
  //but will need the following include commands; #include <ctime> ; #include <cstdlib>
  double prob_of_mode;
  double x;
  double prob;
  //srand(time(0));  // Initialize random number generator implemented outside this function.
  prob = rand();
  prob = prob/RAND_MAX;
  prob_of_mode = pow((mode-min),2)/((max-min)*(mode-min));
  if(prob<=prob_of_mode)x =sqrt(prob*(max-min)*(mode-min))+min;
  if(prob>prob_of_mode)x = max-sqrt((1-prob)*(max-min)*(max-mode));
  return x;
}

//---ROUNDING OFF-----------------------------------------
double round(double x)
{
  if((x-floor(x))>=0.5)
    {return ceil(x);}
  else
    {return floor(x);}
}

//---RANDOM R1 OR R2 VALUES-------------------------------
double fr(double x, double y)
{
	double frReturnValue = 0;
  frReturnValue = gsl_ran_lognormal(r,x,y);
  return frReturnValue;
}


//---READING VARIABLES FILE--------------------------------
void read_in_VARdata(gsl_vector *vars)
{
  double var4 =0;
  FILE *VARdata;
  VARdata = fopen("variables_lambda.txt","r");
  if (VARdata == NULL)
    {
      perror("Error");
    }
	else
    {
      for (int n=0;n<no_vars;n++)
        {
          fscanf(VARdata,"%lf\n",&var4);
          gsl_vector_set(vars,n,var4) ;
        }
    }
  fclose(VARdata);
}


//---READING DATA FILE--------------------------------------
void read_in_HATdata(gsl_matrix *inputdata,int screening,int data_cols)
{
  double var2 =0;
  FILE *HATdata1;
  HATdata1 = fopen("hatdata_lambda.csv","r");
  if (HATdata1 == NULL)
    {
      perror("Error");
    }
  else
    {
      for (int r=0;r<screening;r++)
        {
          for (int c=0;c<data_cols;c++)
            {
              fscanf(HATdata1,"%lf",&var2);
              gsl_matrix_set(inputdata,r,c,var2) ;
            }
        }
    }
  fclose(HATdata1);
}


//---SCALING A LIKELIHOOD HITS DISTRIBUTION SUCH THAT ITS TOTAL IS ONE
void dist_scale(gsl_vector *somedist, int somelength) {
  double x=0;
  for (int y=0; y<somelength; y++)
    {
      x = x + gsl_vector_get(somedist, y);
    }

  if (x > 0) {gsl_vector_scale (somedist, 1/x);}
}


//---RANDOMLY ADJUSTING PASSIVE CASE DETECTIONS FOR STAGING INACCURACY

void staging_adjust(int i, gsl_matrix *finputdata, gsl_matrix *finputdata_temp) {

  int ftmax = gsl_matrix_get(finputdata, i, 1);

  double okclasss1min = gsl_matrix_get(finputdata, i, 2);
  double okclasss1mode = gsl_matrix_get(finputdata, i, 3);
  double okclasss1max = gsl_matrix_get(finputdata, i, 4);
  double okclasss2min = gsl_matrix_get(finputdata, i, 5);
  double okclasss2mode = gsl_matrix_get(finputdata, i, 6);
  double okclasss2max = gsl_matrix_get(finputdata, i, 7);


  //True number of passive detections per month, adjusting for staging misclassification
  int fD1adj = 0;
  int fD2adj = 0;
  int fD1obs = 0;
  int fD2obs = 0;
  int fD1poss = 0;
  int fD2poss = 0;
  int xmin = 0;
  int xmax = 0;
  int ymin = 0;
  int ymax = 0;
  int delta = 0;
  int vectorpos1 = 0;
  int solution1 = 0;
  int Drange = 0;
  double fokclasss1 = rand_triang_dist(okclasss1min, okclasss1mode, okclasss1max);
  double fokclasss2 = rand_triang_dist(okclasss2min, okclasss2mode, okclasss2max);
  double cumprob = 0;
  double rand9 = 0;


  for (int u=0;u<ftmax;u++)
    {
      fD1adj = 0;
      fD2adj = 0;
      fD1obs = gsl_matrix_get(finputdata, i, u+14);
      fD2obs = gsl_matrix_get(finputdata, i, u+14+34);

      if ((fD1obs + fD2obs)>0)
        {
          cumprob = 0;

          Drange = fD1obs + fD2obs + 1;
          gsl_matrix *misclass = gsl_matrix_calloc(Drange, 3);

          for (int v=0;v<Drange;v++)
            {
              fD1poss=fD1obs + fD2obs - v;
              fD2poss=0+v;
              gsl_matrix_set(misclass, v, 0, fD1poss);
              gsl_matrix_set(misclass, v, 1, fD2poss);

              delta = fD2poss - fD2obs;
              xmin = max(0, -delta);
              xmax = min(fD1poss, fD2poss - delta);
              ymin = max(0, delta);
              ymax = min(fD2poss, fD1poss + delta);

              cumprob=0;
              for (int w=xmin; w<(xmax+1); w++)
                {
                  cumprob=cumprob+( pow((fokclasss1), fD1poss-w) * pow((1-fokclasss1), w) * pow((fokclasss2), fD2poss-w+delta) * pow((1-fokclasss2), w+delta) );
                }


              if (v==0) {gsl_matrix_set(misclass,v,2,cumprob);}
              else {gsl_matrix_set(misclass,v,2, gsl_matrix_get(misclass,v-1,2)+cumprob);}

            }

          rand9=gsl_ran_flat(r,0,gsl_matrix_get(misclass, Drange-1, 2));
          vectorpos1=0;
          solution1=0;

          for (int z=0;z<Drange;z++) {
            if (z==0) {if (rand9<=gsl_matrix_get(misclass,z,2)) {vectorpos1=z; solution1=1;}}
            else {if ((rand9>gsl_matrix_get(misclass,z-1,2)) && (rand9<=gsl_matrix_get(misclass,z,2)))
                {vectorpos1=z; solution1=1;} }
            if (solution1==1) break;
          }


          fD1adj = gsl_matrix_get(misclass,vectorpos1,0);
          fD2adj = gsl_matrix_get(misclass,vectorpos1,1);

          gsl_matrix_set(finputdata_temp, i, u+14, fD1adj);
          gsl_matrix_set(finputdata_temp, i, u+14+34, fD2adj);

          gsl_matrix_free(misclass);
        }
    }

}





//---FP STOCHASTIC MODEL----------------------------------
int fp(int i, double flambda, double fr1mean, double fr1SD, double fr2mean, double fr2SD, gsl_matrix *finputdata)
{
  int total=0;

  //Reading variables from the data, including S1 and S2 prevalent and detected, list of D1, D2 and N per month, up to tmax
  int ftmax=gsl_matrix_get(finputdata,i,1);

  int fS1previous=gsl_matrix_get(finputdata,i,8);
  int fS2previous=gsl_matrix_get(finputdata,i,9);
  int fS1start=fS1previous-gsl_matrix_get(finputdata,i,10);
  int fS2start=fS2previous-gsl_matrix_get(finputdata,i,11);
  int fS1end=gsl_matrix_get(finputdata,i,12);
  int fS2end=gsl_matrix_get(finputdata,i,13);

  gsl_vector *fD1 = gsl_vector_alloc(ftmax);
  for (int u=0;u<ftmax;u++)
    {
      gsl_vector_set(fD1,u,gsl_matrix_get(finputdata,i,u+14));
    }

  gsl_vector *fD2 = gsl_vector_alloc(ftmax);
  for (int u=0;u<ftmax;u++)
    {
      gsl_vector_set(fD2,u,gsl_matrix_get(finputdata,i,u+14+34));
    }

  gsl_vector *fN = gsl_vector_alloc(ftmax);
  for (int u=0;u<ftmax;u++)
    {
      gsl_vector_set(fN,u,gsl_matrix_get(finputdata,i,u+14+68));
    }


  //Other needed parameters
  int fS1pred=fS1start;
  int fS2pred=fS2start;
  int ftcum = 0;
  int var33=0;

  //Random lognormal value of r1
  double fr1daily = fr(fr1mean,fr1SD);
  double fr1 = 1-exp(-30.41*fr1daily);
  //Random lognormal value of r2
  double fr2daily = fr(fr2mean,fr2SD);
  double fr2 = 1-exp(-30.41*fr2daily);



  //Implementing model up to tmax
  while(ftcum<ftmax){

    //If statement that decides whether treatment or natural progression occurs first, at random
    if (gsl_ran_flat(r,0,1)<0.5) {

      //S1 and S2 cases detected during month
      fS1pred=fS1pred-gsl_vector_get(fD1,ftcum);
      fS2pred=fS2pred-gsl_vector_get(fD2,ftcum);

      //Out of remaining S2 cases, cases progressing to death
      for (int l=0;l<fS2pred;l++)
        {
          if (gsl_ran_flat(r,0,1)<fr2)
            {
              fS2pred=fS2pred-1;
            }
        }

      //Out of remaining S1 cases, cases progressing to S2
      for (int m=0;m<fS1pred;m++)
        {
          if (gsl_ran_flat(r,0,1)<fr1)
            {
              fS1pred=fS1pred-1;
              fS2pred=fS2pred+1;
            }
        }
    }
    else {

      //Out of remaining S2 cases, cases progressing to death
      for (int l=0;l<fS2pred;l++)
        {
          if (gsl_ran_flat(r,0,1)<fr2)
            {
              fS2pred=fS2pred-1;
            }
        }

      //Out of remaining S1 cases, cases progressing to S2
      for (int m=0;m<fS1pred;m++)
        {
          if (gsl_ran_flat(r,0,1)<fr1)
            {
              fS1pred=fS1pred-1;
              fS2pred=fS2pred+1;
            }
        }

      //S1 and S2 cases detected during month
      fS1pred=fS1pred-gsl_vector_get(fD1,ftcum);
      fS2pred=fS2pred-gsl_vector_get(fD2,ftcum);

    }


    //Addition of new incident S1 cases during the present month
    fS1pred=fS1pred+gsl_ran_binomial(r, flambda, (gsl_vector_get(fN,ftcum)-fS1pred-fS2pred));

    //Breaking early if prevalence reaches 100%; done mainly to speed up model, esp. since the result in this case  is always going to be a no hit
    if ((fS1pred+fS2pred)>=(gsl_vector_get(fN,ftcum)))
      {ftcum=ftmax;
        var33=1;
      }

    //Breaking early if prevalence becomes negative:
    if ((fS1pred<0) || (fS2pred<0))
      {ftcum=ftmax;
        var33=1;
      }

    //Advancing time
    ftcum=ftcum+1;
    //Breaking when time reaches the entire analysis period
    if (ftcum>=ftmax) break;
  }

  //Checking whether model predictions fit observed data
  if ((fS1pred==fS1end) && (fS2pred==fS2end) && (var33==0)) {total = 1;} else {total=0;}


  gsl_vector_free(fD1);
  gsl_vector_free(fD2);
  gsl_vector_free(fN);

  return total;
}



//---WRITING OUTPUT-----------------------------------

//Function to write matrix output
void Fun_outputmatrix(gsl_matrix *scoresmatrix, int lambda_n, int screening, int index)
{
  FILE *outputmatrix;
  char* filename(0);
  sprintf(filename, "output_lambda%d.csv", index);
  outputmatrix = fopen(filename,"w");
  if (outputmatrix == NULL) perror ("Error");
  for(int r=0;r<(1+lambda_n);r++)
    {
      for(int c=0;c<(1+screening);c++)
        {
          fprintf(outputmatrix,"%6.4f ", gsl_matrix_get(scoresmatrix,r,c));
        }
      fprintf(outputmatrix,"\n");
    }
  fclose(outputmatrix);
}

//Function to write vector output
void outputvector(gsl_matrix *inputdata, gsl_vector *resultvector, int veclength, const char *filename)
{
  FILE *myfile;
  myfile = fopen(filename,"w");
  if (myfile == NULL) perror ("Error");
  fprintf(myfile, " result");
  for(int r=0;r<veclength;r++)
    {
      fprintf(myfile," %6.0f",gsl_matrix_get(inputdata,r,0));
    }
  fprintf(myfile,"\n");
  for(int r=0;r<veclength;r++)
    {
      fprintf(myfile,"%12.4f", gsl_vector_get(resultvector,r));
    }
  fprintf(myfile,"\n");
  fclose(myfile);
  return;
}



//---IMPLEMENTATION OF FUNCTIONS----------------------
//-----------------------------------------------------

int main(int argc, char** argv)
{
  int index = atoi(argv[1]);

  //---READING IN PARAMETERS FROM VARIABLES FILE-----------
  gsl_vector *vars = gsl_vector_calloc(no_vars);
  read_in_VARdata(vars);
  int iterations = gsl_vector_get(vars,0);

  double lambda_min = gsl_vector_get(vars,1);
  double lambda_int1 = gsl_vector_get(vars,2);
  double lambda_node1 = gsl_vector_get(vars,3);
  double lambda_int2 = gsl_vector_get(vars,4);
  double lambda_node2 = gsl_vector_get(vars,5);
  double lambda_int3 = gsl_vector_get(vars,6);
  double lambda_node3 = gsl_vector_get(vars,7);
  double lambda_int4 = gsl_vector_get(vars,8);
  double lambda_node4 = gsl_vector_get(vars,9);
  double lambda_int5 = gsl_vector_get(vars,10);
  double lambda_max = gsl_vector_get(vars,11);

  double r1mean = gsl_vector_get(vars,12);
  double r1SD = gsl_vector_get(vars,13);
  double r2mean = gsl_vector_get(vars,14);
  double r2SD = gsl_vector_get(vars,15);
  int screening = gsl_vector_get(vars,16);
  int data_cols = gsl_vector_get(vars,17);

  int runs = gsl_vector_get(vars, 18);
  int runs_max = gsl_vector_get(vars, 19);
  int hits_target = gsl_vector_get(vars, 20);
  int hits_min = gsl_vector_get(vars, 21);

  double var3=0;


  //---CREATING UNKNOWN PARAMETER VECTORS-----------

  //Creating vector of the range of lambda
  int lambda_n1=ceil(((lambda_node1-lambda_min)/lambda_int1));
  int lambda_n2=ceil(((lambda_node2-lambda_node1)/lambda_int2));
  int lambda_n3=ceil(((lambda_node3-lambda_node2)/lambda_int3));
  int lambda_n4=ceil(((lambda_node4-lambda_node3)/lambda_int4));
  int lambda_n5=ceil(((lambda_max-lambda_node4)/lambda_int5)+1);
  int lambda_n=lambda_n1+lambda_n2+lambda_n3+lambda_n4+lambda_n5;
  gsl_vector *lambdarange = gsl_vector_calloc(lambda_n);

  var3=lambda_min;
  for(int i=0;i<lambda_n1;i++)
    {
      gsl_vector_set(lambdarange,i,var3);
      var3 = var3+lambda_int1;
    }
  var3=lambda_node1;
  for(int i=lambda_n1;i<(lambda_n1+lambda_n2);i++)
    {
      gsl_vector_set(lambdarange,i,var3);
      var3 = var3+lambda_int2;
    }
  var3=lambda_node2;
  for(int i=(lambda_n1+lambda_n2);i<(lambda_n1+lambda_n2+lambda_n3);i++)
    {
      gsl_vector_set(lambdarange,i,var3);
      var3 = var3+lambda_int3;
    }
  var3=lambda_node3;
  for(int i=(lambda_n1+lambda_n2+lambda_n3);i<(lambda_n1+lambda_n2+lambda_n3+lambda_n4);i++)
    {
      gsl_vector_set(lambdarange,i,var3);
      var3 = var3+lambda_int4;
    }
  var3=lambda_node4;
  for(int i=(lambda_n1+lambda_n2+lambda_n3+lambda_n4);i<lambda_n;i++)
    {
      gsl_vector_set(lambdarange,i,var3);
      var3 = var3+lambda_int5;
    }



  //---SEED FOR RANDOM NUMBER GENERATOR
  long seed;//Declare a variable to hold the seed
  seed = time(NULL);
  const gsl_rng_type * T;//sets a pointer to the RNG type
  T = gsl_rng_default;//Sets the rng to the default type
  r = gsl_rng_alloc (T);//allocate the rng to the pointer
  gsl_rng_set(r,seed);//set the seed for the rng


  //---CREATING COMPONENTS OF INPUT FILES-----------
  gsl_matrix *inputdata = gsl_matrix_calloc(screening,data_cols);
  gsl_matrix *inputdata_temp = gsl_matrix_calloc(screening,data_cols);


  //---CREATING COMPONENTS OF OUTPUT FILE-----------

  //Temporary output vector
  gsl_vector *scoresvector = gsl_vector_calloc(lambda_n);

  //As a matrix with column 1=lambda value, other columns=screenings
  gsl_matrix *scoresmatrix = gsl_matrix_calloc(lambda_n+1, screening+1);

  //Attributing lambdarange as the first column of scoresmatrix
  for (int v=0;v<lambda_n;v++) {gsl_matrix_set(scoresmatrix,v+1,0,gsl_vector_get(lambdarange,v));}

  //Attributing screening ID numbers as the first row of scoresmatrix
  for (int z=0;z<screening;z++) {gsl_matrix_set(scoresmatrix,0,z+1,gsl_matrix_get(inputdata,z,0));}

  //Output vector to track for how many iterations time series fail to yield a solution
  gsl_vector *solutions = gsl_vector_calloc(screening);


  //---READING IN DATA-------------------------------------
  read_in_HATdata(inputdata,screening,data_cols);
  read_in_HATdata(inputdata_temp, screening, data_cols);



  //---IMPLEMENTING FP-----------------------------------------

  int hits=0;
  int runs_cum=0;

  //for each iteration (each iteration = one set of random D1 and D2 values for each village time series)
  for (int k1=0; k1<iterations; k1++)
    {
      printf("Now running iteration: %4i  \n",k1);

      //Resetting some parameters to zero
      hits = 0;

      //for each village time series, randomly adjusting D1 and D2 for staging misclassification...
      for (int i=0; i<screening; i++)
        {
          staging_adjust(i, inputdata, inputdata_temp);
        }

      //for each village time series...
      for (int i=0; i<screening; i++)
        {
          printf("Now modelling time series: %4i  id: %8.0f\n", i, gsl_matrix_get(inputdata, i, 0));

          gsl_vector_set_zero(scoresvector);

          runs_cum = 0;
          hits = 0;

          //for each model run (each run = running model once in a given village time series for all possible lambda values)
          for(int k2=0; k2<runs; k2++)
            {
              runs_cum = runs_cum + 1;

              //for each candidate lambda value
              for(int j=0; j<lambda_n; j++)
                {

                  lambda=gsl_vector_get(lambdarange,j);

                  score=fp(i, lambda, r1mean, r1SD, r2mean, r2SD, inputdata_temp);

                  if (score==1)
                    {
                      hits = hits + 1;
                      gsl_vector_set(scoresvector, j, gsl_vector_get(scoresvector, j)+1);
                    }
                }

              //Breaking early if at least the target number of likelihood hits is reached for the time series
              if (hits > (hits_target - 1)) break;

            }


          //If the runs were not enough to reach the minimum number of likelihood hits at least, then go on until hits_min are reached
          if (hits < hits_min) {

            while (runs_cum < runs_max) {

              runs_cum = runs_cum + 1;

              //for each candidate lambda value
              for(int j=0; j<lambda_n; j++)
                {

                  lambda=gsl_vector_get(lambdarange,j);

                  score=fp(i, lambda, r1mean, r1SD, r2mean, r2SD, inputdata_temp);

                  if (score==1)
                    {
                      hits = hits + 1;
                      gsl_vector_set(scoresvector, j, gsl_vector_get(scoresvector, j)+1);
                    }
                }

              //Breaking early if at least the minimum number of likelihood hits is reached for the time series
              if (hits > (hits_min - 1)) break;

            }
          }

          //Scaling the likelihood distribution for the time series so that its maximum is one, but only if there is any hit
          dist_scale(scoresvector, lambda_n);

          cout<<" runs_cum= "<<runs_cum<<" \n";

          //Adding the vector results for the time series to the overall output matrices
          for (int j1=0; j1<lambda_n; j1++)
            {
              gsl_matrix_set(scoresmatrix, j1+1, i+1, gsl_matrix_get(scoresmatrix, j1+1, i+1) + gsl_vector_get(scoresvector, j1) );
            }

          //Updating solutions vector
          if (hits > 0) {gsl_vector_set(solutions, i, gsl_vector_get(solutions, i)+1);}
        }
    }





  //---WRITING OUTPUT FILES-------------------------------------
  Fun_outputmatrix(scoresmatrix, lambda_n, screening, index);
  char* filename(0);
  sprintf(filename, "solutions%d.csv", index);
  outputvector(inputdata, solutions, screening, filename);

  //---CLEANING UP----------------------------------------------
  gsl_matrix_free(inputdata);
  gsl_matrix_free(inputdata_temp);
  gsl_matrix_free(scoresmatrix);
  gsl_vector_free(scoresvector);
  gsl_vector_free(lambdarange);
  gsl_vector_free(vars);

  return 0;
}
