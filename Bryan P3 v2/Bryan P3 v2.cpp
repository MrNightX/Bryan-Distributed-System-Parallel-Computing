// Bryan P3 v2.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <stdio.h>
#include <omp.h>
#include <iomanip>

using namespace std;

static long num_steps = 100000000;

double Ori_Pi(int numsteps)
{
    int i;
    double x, pi, sum = 0.0;
    double start_time, run_time;
    double step;
    step = 1.0 / (double)num_steps;
    start_time = omp_get_wtime();

    for (i = 1; i <= num_steps; i++) {
        x = (i - 0.5) * step;
        sum = sum + 4.0 / (1.0 + x * x);
    }

    pi = step * sum;
    run_time = omp_get_wtime() - start_time;
    cout << "Ori Code" << endl;
    cout << "Thread ID: " << omp_get_thread_num() << endl;
    cout << "Pi with " << num_steps << setprecision(16) << " steps is "
        << pi << " in " << run_time << " seconds\n";
    return run_time;
}

void find_pi_multi(int num_steps)
{
    const int TOTAL_THREADS = 16;
    double step;

    double pi = 0.0;
    double start_time, run_time;

    step = 1.0 / (double)num_steps;
    double sum[TOTAL_THREADS] = { 0.0 };

    start_time = omp_get_wtime();

#pragma omp parallel num_threads(TOTAL_THREADS)
    {
        int id = omp_get_thread_num();
        double partial_sum = 0.0;
        double x = 0.0;

        for (int i = id; i < num_steps; i += TOTAL_THREADS) 
        {
            x = (i + 0.5) * step;
            partial_sum += 4.0 / (1.0 + x * x);
        }

        sum[id] = partial_sum;
    }

    for (int i = 0; i < TOTAL_THREADS; i++) 
    {
        pi += step * sum[i];
    }

    run_time = omp_get_wtime() - start_time;
    cout << "\n pi with " << num_steps << setprecision(16)
        << " step is " << pi
        << " in " << run_time << " seconds\n";
    const int Padding_Size = 64;
    

}

float find_pi_multi_V1(int numstep)
{
    //Insert Process Here
    const int Total_Threads = 16;
    const int Padding_Size = 64;
    double PartialSum[Total_Threads][Padding_Size];

    double pi = 0.0;
    double step;
    step = 1.0 / (double)numstep;

    double StartTime, EndTime = 0.0;
    StartTime = omp_get_wtime(); //start time

#pragma omp parallel num_threads(Total_Threads)
    {
        
        int id = omp_get_thread_num();
        double partial_sum = 0.0;
        double x = 0.0;

        for (int i = id; i < num_steps; i += Total_Threads) 
        {
            x = (i + 0.5) * step;
            partial_sum += 4.0 / (1.0 + x * x);
        }
        
        /*
        for(i stuff)
        {
            PartialSum[i][0] = 0;
        }
        */
        //double PartialSum[Total_Threads][64]
        //sum[i][64];

        //sum[id] = partial_sum;
    }

    for (int i = 0; i < Total_Threads; i++)
    {
        //pi += step * sum[i];
    }

    EndTime = omp_get_wtime() - StartTime; //end time
    return EndTime;
}

double find_pi_multi_V2(int numstep)
{
    double StartTime, EndTime = 0.0;
    StartTime = omp_get_wtime();

    //Insert Process here
#pragma omp parallel for reduction(+ : sum)
    {

    }

    EndTime = omp_get_wtime() - StartTime;
    return EndTime;
}

// I DONT FUCKING GET IT AT ALL
// WHERE IS THE SAMPLE? HOW THE FUCK CAN WE SOLVE THIS?
// ALL THEORY BUT NO PRACTICAL LESSON

//pragma omp parallel for reduction(+:sum)
int main()
{
    std::cout << "Hello World!\n";
    Ori_Pi(num_steps);
    //find_pi_multi_V1(num_steps);
    find_pi_multi(num_steps);
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
