#include<stdio.h>
#include<iostream>
#include<math.h>
#include<string>
#include<random>
#include<map>
#include<chrono>
#include<vector>
#include<fstream>
#include<stdlib.h>
#include<algorithm>
#define PI 3.14159265
/*****Some Global Variables**********/
int arrctr;
int *ptr1, *ptr2;
std::vector<std::string> arr; //arr stores all the grey codes
std::map<std::string, double> constellation; // maps all the grey codes to phases
/*****Function Prototypes**********/
int check_validity(int );
void adjust_length(int *, int); // adjusts the buffer such that it is integer multiple of log_2(M)
void greycoding(int); // generates log_2_(M) bit grey codes with total # of M
void modulation(); // assigns every codeword an angle
void generate_input(int *, int); // generates input 0 or 1
void add_noise(int *, int *, int, int); // adds AWGN noise
std::string map_estimator(double, double, int);
void put_in_array(int *output, std::string ans);
void probability_error(int, int *, int *, int, int, int);
void freeall(int *, int *); // deletes all dynamic alloted memory
void print();
void printarr(int *, int);
/***************/
int main ()
{
    int M, num, *input, *output;
    printf("Enter the modulation scheme (M-ary) which you would like to simulate\n");
    scanf("%d", &M);
    while (check_validity(M) != 1)
    {
        printf("Please enter the valid modulation scheme which is a power of 2\n");
        scanf("%d", &M);
    }
    printf("Enter the number of bits you would like the source to send\n");
    scanf("%d", &num);
    adjust_length(&num, M);
    printf("the new value of number of bits is %d\n", num);
    input = (int *)malloc(num*sizeof(int));
    output = (int *)malloc(num*sizeof(int));
    greycoding(log(M)/log(2));
    modulation();
    print();
    //generate input and add noise
    generate_input(input, num);
    //added to debug the code
    // printarr(input, num);
    add_noise(input, output, num, M);
    freeall(input ,output);
    return 0;
}
int check_validity(int M)
{
    if (M == 0)
        return 0;
    while (M != 1)
    {
        if (M%2 != 0)
            return 0;
        M = M/2;
    }
    return 1;
}
void adjust_length(int *ptr, int M)
{
    int test = log(M)/log(2);
    while((*ptr)%test !=0)
        *ptr = *ptr+1;
}
void greycoding(int n)
{
    if(n<=0)
        return;
    // start with one bit grey code
    arr.push_back("0");
    arr.push_back("1");
    int i, j;
    //Every iteration of this loop generates 2*i codes
    for(i=2; i< (1<<n); i=i<<1)
    {
        // previously generated codes are again in arr[] in reverse order
        for(j=i-1; j>=0; j--)
            arr.push_back(arr[j]);
        // append 0 to the first half
        for (j = 0 ; j < i ; j++)
            arr[j] = "0" + arr[j];
        // append 1 to the second half
        for (j = i ; j < 2*i ; j++)
            arr[j] = "1" + arr[j];
    }
}
void modulation()
{
    int i;
    for(i=0; i<arr.size(); i++)
    {
        constellation[arr[i]] = 360*i/(arr.size());
    }
}
void print()
{
    for(std::map<std::string, double>::iterator it = constellation.begin(); it !=constellation.end(); ++it)
    {
        std::cout << it->first <<"\t" << it->second;
        printf("\n");
    }
}
void generate_input(int *input, int num)
{
    int i, bits;
    for(i=0;i<num;i++)
    {
        bits = rand()%2;
        input[i] = bits;
    }
}
void add_noise(int *input, int *output, int num, int M)
{
    int i, j, n_sim, snr;
    int *ptr1, *ptr2;
    double noise, stddev, angle, c1, c2;
    int symbol_length = log(M)/log(2);
    std::string temp;
    std::string ans = "";
    printf("Enter how many times you want to run monte carlo simulations?\n");
    scanf("%d", &n_sim);
    for(i=0;i<n_sim;i++)
    {
        snr = -50; //units in dB
        while(snr < 50)
        {
            temp = "";
            arrctr = 0;
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator(seed);
            stddev = sqrt(1.0 / pow(10.0, snr/10.0));
            std::normal_distribution<double> distribution(0.0, stddev);
            noise = distribution(generator);
            ptr1 = input;
            ptr2 = input+(symbol_length-1);
            while(ptr2 <= input+(num-1))
            {
                while(ptr1 <= ptr2)
                {
                    temp += std::to_string(*ptr1);
                    ptr1++;
                }
                angle = constellation[temp];
                c1 = cos(angle*PI/180.0);
                c2 = sin(angle*PI/180.0);
                c1 += noise;
                c2 += noise;
                ans = map_estimator(c1, c2, M);
                put_in_array(output, ans);
                temp = "";
                ptr2 += symbol_length;
            }
            probability_error(snr, input, output, num, log(M)/log(2), i);
            snr += 5;
        }
    }
}
std::string map_estimator(double c1, double c2, int M)
{
    std::string ans = "";
    double min_dist, temp, c1_i, c2_i;
    int i;
    c1_i = cos(PI*constellation[arr[0]]/180.0);
    c2_i = sin(PI*constellation[arr[0]]/180.0);
    min_dist = sqrt((c1-c1_i)*(c1-c1_i)+(c2-c2_i)*(c2-c2_i));
    ans=arr[0];
    for(i=1; i<M; i++)
    {
        c1_i = cos(PI*constellation[arr[i]]/180.0);
        c2_i = sin(PI*constellation[arr[i]]/180.0);
        temp = sqrt((c1-c1_i)*(c1-c1_i)+(c2-c2_i)*(c2-c2_i));
        if (temp < min_dist)
        {
            min_dist = temp;
            ans = arr[i];
        }
    }
    return ans;
}
void put_in_array(int *output, std::string ans)
{
    int i;
    for(i=0; i<ans.length(); i++)
    {
        output[arrctr] = ans[i] - '0';
        arrctr++;
    }
}
void probability_error(int snr, int *input, int *output, int num, int n, int i)
{
    int j;
    int k = 0;
    int bitctr = 0;
    int symctr = 0;
    int nsymbols = num/n;
    double p_be, p_se;
    ptr1 = (int *)malloc(n*sizeof(int));
    ptr2 = (int *)malloc(n*sizeof(int));
    for(j=0;j<num;j++)
    {
        if (input[j] != output[j])
            bitctr++;
    }
    for(j=0;j<num;j++)
    {
        if(k<=n)
        {
            ptr1[k] = input[j];
            ptr2[k] = output[j];
            k++;
            if(k == n)
            {
                if (std::equal(ptr1, ptr1+n, ptr2) == false)
                    symctr++;
                k = 0;
            }
        }
    }
    p_be = bitctr/(double)num;
    p_se = symctr/(double)nsymbols;
    std::ofstream myfile;
    std::string name="output" + std::to_string(i) + ".txt";
    myfile.open (name, std::fstream::app);
    myfile << snr<<"\t"<<p_be<<"\t"<<p_se<<std::endl;
    myfile.close();
}
void freeall(int *input, int *output)
{
    free(input);
    free(output);
    free(ptr1);
    free(ptr2);
}
void printarr(int *input, int num)
{
    int i;
    for(i=0;i<num;i++)
    {
        printf("%d", input[i]);
    }
    printf("\n");
}