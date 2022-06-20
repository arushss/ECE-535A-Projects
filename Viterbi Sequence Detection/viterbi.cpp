#include<stdio.h>
#include<stdlib.h>
#include<random> // use later
#include<string.h>
#include<chrono>
#include<time.h>
#include<float.h>
#include<fstream>
#include<algorithm>
#include<math.h>
/**Some global variables**/
FILE *fp1;
int diffcounter;
int nstates = 1;
int *finalstatetransition, *viterbinp, *test;
int i, j, l, k, num, _index, nsim, n_i;
double snr, pmetric, p_e, soutput, stddev, noise;
double *nrm;
char ch;
/*************************/
/*Function prototypes*/
void description();
void reserve_size(int **ptr); // reserve 2d memory in the array passed as an argument
int numstates(); // this function calculates how many states I have in my software
void fill_array(FILE *fp, int** ptr); //fills in the array
void startviterbi(int **smatrix, int **imatrix, int **omatrix);
void calculateprobability(int **imatrix, double **pathmatrix, int **statematrix, int *input, int* viterbinp, double snr, int n_i);
//deletes all the dynamically allocated memory at the end of the simulations
void deleteall(int *input, int **smatrix, int **imatrix, int **omatrix, double **pathmatrix, int **statematrix);
/*********************/
int main ()
{
    FILE *fp2, *fp3;
    fp2 = fopen("input_transitions.txt", "r");
    if (fp2 == NULL)
    {
        printf("Check input_transitions.txt file\n");
        return -1;
    }
    fp3 = fopen("output_transitions.txt", "r");
    if (fp3 == NULL)
    {
        printf("Check output_transitions.txt file\n");
        return -1;
    }
    description();
    nstates=numstates();
     // 2d arrays each for storing state transition, input transitions and output transitions respectively from three files
    int *smatrix[nstates], *imatrix[nstates], *omatrix[nstates];
    reserve_size(smatrix);
    reserve_size(imatrix);
    reserve_size(omatrix);
    //fill in the state transition in smatrix
    fill_array(fp1, smatrix);
    //fill in the input transition in imatrix
    fill_array(fp2, imatrix);
    //fill in the output transition in omatrix
    fill_array(fp3, omatrix);
    startviterbi(smatrix, imatrix, omatrix);
    fclose(fp1);
    fclose(fp2);
    fclose(fp3);
    return 0;
}
void description()
{
    printf("This program consists of three input files:\n");
    printf("state_connections text file consists of 1s and 0s where 1 indicates connectivity\n");
    printf("and 0 indicates disconnectivity\n");
    printf("input_transitions text file has 0/1 input transition from one state to another.\n");
    printf("output transitions indicates output transitions from one state to another\n");
    printf("999 in both the above files indicate disconnectivity between those two states\n");
}
int numstates()
{
    fp1 = fopen("state_connections.txt", "r");
    if (fp1 == NULL)
    {
        printf("Check state_connections.txt file\n");
        return -1;
    }
    while(!feof(fp1))
    {
        if(ch=fgetc(fp1) == '\n')
            nstates++;
    }
    fseek(fp1, 0, SEEK_SET);
    return nstates;
}
void reserve_size(int **ptr)
{
    for (i=0; i<nstates; i++)
        ptr[i]=(int *)malloc(nstates*sizeof(int));
}
void fill_array(FILE *fp, int** ptr)
{
    while(!feof(fp))
    {
        for(i=0;i<nstates;i++)
        {
            for(j=0;j<nstates;j++)
            {
                fscanf(fp, "%d", &num);
                ptr[i][j] = num;
            }
        }
    }
}
void startviterbi(int **smatrix, int **imatrix, int **omatrix)
{
    k=rand()%nstates; // randomly picks any state as a starting state from finite state machine
    printf("How many bits of information do you want to send from source?\n");
    scanf("%d", &num);
    // input array will be used in calculating probability of error
    int *input = (int*)malloc(num*sizeof(int));
    // this will show us maximum likelihood states at the end of viterbi
    finalstatetransition = (int*)malloc((num+1)*sizeof(int));
    // based on max likelihood states, this will output the decoded input
    viterbinp = (int*)malloc(num*sizeof(int));
    // now reserve the space for path metrics and previous state matrix
    double *pathmatrix[nstates]; // 2d array to store path metric of every state in all iterations
    int *statematrix[nstates]; // 2d array to store previous state in all iterations
    for(i=0;i<nstates;i++)
    {
        pathmatrix[i] = (double *)malloc((num+1)*sizeof(double));
        statematrix[i] = (int *)malloc(num*sizeof(int));
    }
    test=(int*)malloc(nstates*sizeof(int)); // used to store intermediate results
    nrm=(double*)malloc(nstates*sizeof(double)); // used to store intermediate results
    printf("Enter how many times you want to run the monte carlo simulations?\n");
    scanf("%d", &nsim);
    for(n_i=0;n_i<nsim;n_i++)
    {
        //noise should be added here
        snr = -500; //consider this having units as dB
        while (snr<60)
        {
            diffcounter=0;
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator(seed);
            stddev = sqrt(1.0/pow(10, snr/10.0));
            std::normal_distribution<double> distribution(0.0, stddev);
            noise = distribution(generator);
            //everything else needs to be added here
            //fill the path matrix with INFINITY except the beginning state
            for(i=0;i<nstates;i++)
            {
                for(j=0;j<num+1;j++)
                {
                    pathmatrix[i][j]=INFINITY;
                }
            }       
            pathmatrix[k][0]=0;
            //fill the statematrix with some redundant value
            for(i=0;i<nstates;i++)
            {
                for(j=0;j<num;j++)
                {
                    statematrix[i][j]=99;
                }
            }
            for(i=0;i<num;i++)
            {
                for(j=0;j<nstates;j++)
                {
                    test[j] = smatrix[k][j];
                }
                // randomly pick the element from test matrix
                _index = rand()%nstates;
                while(test[_index]==0)
                    _index = rand()%nstates;
                input[i] = imatrix[k][_index]; // store the input sequence bit by bit because we need to calculate probability of error
                soutput = omatrix[k][_index] + noise; // soutput gives output from the finite state machine. we need to add noise part
                k = _index;
                //iteration in context of path matrix starts from here (very main loop)
                for(j=0;j<nstates;j++)
                {
                    if(pathmatrix[j][i] == INFINITY)
                        continue;
                    else
                    {
                        //scan the connections for every state
                        for(l=0;l<nstates;l++)
                        {
                            if(smatrix[j][l]==0)
                                continue;
                            else
                            {
                                pmetric=pathmatrix[j][i]+pow(omatrix[j][l] - soutput,2);
                                if(pathmatrix[l][i+1]>pmetric)
                                {
                                    pathmatrix[l][i+1]=pmetric;
                                    statematrix[l][i]=j;
                                }
                                // normalize the pathmatrix in order to avoid overflow conditions
                                if(i!=0 && i%100==0)
                                {
                                    int m;
                                    double min;
                                    for(m=0;m<nstates;m++)
                                    {
                                        nrm[m] = pathmatrix[m][i];
                                    }
                                    for(m=0; m<nstates; m++)
                                    {
                                        min = nrm[0];
                                        if(nrm[m] <= min)
                                            min = nrm[m];
                                    }
                                    for(m=0;m<nstates;m++)
                                    {
                                        pathmatrix[m][i] -= min;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            // do the calculation of probability by traceback
            calculateprobability(imatrix, pathmatrix, statematrix, input, viterbinp, snr, n_i);
            snr = snr + 5;
        }
    }
    //should delete all the dynamically allocated memory
    deleteall(input, smatrix, imatrix, omatrix, pathmatrix, statematrix);
}
void calculateprobability(int **imatrix, double **pathmatrix, int **statematrix, int *input, int* viterbinp, double snr, int n_i)
{
    //get the minimum element from path matrix last column
    int index_;
    double min = pathmatrix[0][num];
    for(i=0;i<nstates;i++)
    {
        if(pathmatrix[i][num]<=min)
        {
            min=pathmatrix[i][num];
            index_=i;
        }
    }
    //extract all the information from state matrix
    int counter = num;
    int statevalue;
    finalstatetransition[num]=index_;
    while(counter!=0)
    {
        statevalue=statematrix[index_][counter-1];
        finalstatetransition[counter-1]=statevalue;
        index_=statevalue;
        counter--;
    }
    for(i=0;i<num;i++)
    {
        viterbinp[i]=imatrix[finalstatetransition[i]][finalstatetransition[i+1]];
    }
    //calculating probability of error
    for(i=0;i<num;i++)
    {
        if (viterbinp[i] != input[i])
            diffcounter++;
    }
    p_e = (double)diffcounter/num;
    //putting everything in a file
    std::ofstream myfile;
    std::string name="output" + std::to_string(n_i) + ".txt";
    myfile.open (name, std::fstream::app);
    myfile << snr << "\t" << p_e << std::endl;
    myfile.close();
}
void deleteall(int *input, int **smatrix, int **imatrix, int **omatrix, double **pathmatrix, int **statematrix)
{
    int i;
    free(finalstatetransition);
    free(viterbinp);
    free(test);
    free(input);
    free(nrm);
    for(i=0;i<nstates;i++)
    {
        free(smatrix[i]);
        free(imatrix[i]);
        free(omatrix[i]);
        free(pathmatrix[i]);
        free(statematrix[i]);
    }
}