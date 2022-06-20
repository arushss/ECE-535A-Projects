#include<stdio.h>
#include<math.h>
#include<random>
#include<chrono>
#include<stdlib.h>
#include<time.h>
#include<fstream>

void generate_inputs(int*, int);
void add_noise(int*, double*, int);
void probability_error(int*, double*, int, int, int);
int main ()
{
    int Nbits, *arr;
    double *outarr;
    printf("Enter the total number of bits you would want source to send\n");
    scanf("%d", &Nbits);
    arr = (int *)malloc(Nbits*sizeof(int));
    outarr = (double *)malloc(Nbits*sizeof(double));
    generate_inputs(arr, Nbits);
    add_noise(arr,outarr,Nbits);
    free(arr);
    free(outarr);
    return 0;
}
void generate_inputs(int *arr, int Nbits)
{
    int i, bipolar;
    for(i=0; i<Nbits; i++)
    {
        bipolar=rand()%2;
        if (bipolar==0)
            bipolar=-1;
        arr[i]=bipolar;
    }
}
void add_noise(int *arr, double *outarr, int Nbits)
{
    int i, j, snr, counter;
    double noise, stddev;
    printf("Enter how many times you want to run the monte carlo simulations?\n");
    scanf("%d", &counter);
    for(i=0; i<counter; i++)
    {
        snr = -500; //considering the units as dB
        while (snr<55)
        {
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator(seed);
            stddev = sqrt(1.0 / pow(10.0, snr/10.0));
            std::normal_distribution<double> distribution(0.0, stddev);
            noise = distribution(generator);
            for(j=0; j<Nbits; j++)
            {
                outarr[j] = arr[j] + noise;
            }
            probability_error(arr, outarr, Nbits, snr, i);
            snr = snr + 5;
        }
    }
}
void probability_error(int *arr, double* outarr, int Nbits, int snr, int i)
{
    int counter = 0;
    int j;
    double p_e, q_val;
    FILE *fp;
    for(j=0; j<Nbits; j++)
    {
        if(outarr[j]<0)
            outarr[j] = -1;
        else
            outarr[j] = 1;
    }
    for(j=0; j<Nbits; j++)
    {
        if (arr[j]*outarr[j]<0)
            counter++;
    }
    q_val = 0.5*erfc(sqrt(pow(10.0, snr/10.0)));
    p_e = (double)counter/Nbits;
    counter = 0;
    std::ofstream myfile;
    std::string name="output" + std::to_string(i) + ".txt";
    myfile.open (name, std::fstream::app);
    myfile << snr<<"\t"<<p_e<<"\t"<<q_val<<std::endl;
    myfile.close();
}