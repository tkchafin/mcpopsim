/* 
 * File:   mcpopsim.cpp
 * Author: Tyler K. Chafin
 *
 * Created on July 15, 2015, 9:10 AM
 */

#include <cstdlib>
#include <set>
#include <iostream>
#include <fstream>
#include "locus.h"
#include <ctime>
#include <vector>
#include <climits>
#include "event.h"
#include "mpi.h"

using namespace std;

void MCgeneration(locus &loc); 
void MPI_parseArgs(int rank, int argc, char* argv[], string* in, string* out, int* num, double* mean, int* samp, int* len, int* popsize, unsigned int* rng); 
void show_usage(std::string name); 


int main(int argc, char* argv[]) {

   /*notes 

    * Need to add some error handling
    * 
    * Need to go back and add comments
    * 
    * */
    
    std::string infile, out; 
    int* seeds; 
    int numloci=10; 
    double mumean = 0.0000001;
    int num_samples=10; 
    int seqlength=100; 
    int popsize = 1000; 
    unsigned int seed = time(0);
    int my_rank, p, start, end; 
    
    MPI::Init(argc, argv);
    my_rank = MPI::COMM_WORLD.Get_rank(); //Get rank
    MPI_Comm_size(MPI::COMM_WORLD, &p); //Get number of processes
    
    /*Parse command-line arguments*/
    MPI_parseArgs(my_rank, argc, argv, &infile, &out, &numloci, &mumean, &num_samples, &seqlength, &popsize, &seed); 
    
    event event_table(infile); 
    srand(seed); 
    int local_seed; 
    
    if (my_rank ==0 ){
        event_table.testPrintTable(); 
        /*Create a seed for each locus*/ 
        seeds = new int[p];
        for (int i=0; i<p; i++)
            seeds[i] = rand(); 
        start = 1; 
        end = ((numloci/p)+(20%p));
    }else{ 
        /*Not sure these do what I want, check later*/
        start = ((((numloci/p)+((numloci/p)*(my_rank-1)))+(numloci%p))+1);
        end = ((((numloci/p)+((numloci/p)*my_rank))+(numloci%p)));
    }
    
    /*All ranks wait to assign local seed
     prevents read from seeds before rank 0 populates*/
    //MPI_SCATTER THE SEEDS!!!
    MPI_Bcast(seeds, p, MPI_INT,0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    local_seed = seeds[my_rank];
    //cout << my_rank << " starts at " << start << " and ends at " << end << endl; 
    //cout << my_rank <<" has seed "<<local_seed<<endl; 
    for (int i=start; i<=end; i++){
        srand(local_seed); 
        local_seed = rand(); 
        locus testlocus(locus::logPoissonMu(mumean,rand()), seqlength, popsize,0.5 ,i,local_seed);
        //locus testlocus(0.00001, 80, 1000, 0.5, i);

        cout << "Locus "<< i << " - mu is " << testlocus.getMu() << endl;


        for (int k=0; k<=event::END;k++)
            MCgeneration(testlocus);
   
        testlocus.finalize(); 
    
        cout << flush; 
        //cout <<"\nEnding simulation at time " << testlocus.getGen()-1 << "\n\n";
    }
    
    if (my_rank==0){ 
        delete[] seeds; 
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI::Finalize();
    return 0;
}

void MCgeneration(locus &loc){ 
    loc.eventParser(); 
    loc.MCpermuteAllPops();  
    loc.incrementGen(); 
    loc.parseMigMatrix();
}

void MPI_parseArgs(int rank, int argc, char* argv[], string* in, string* out, int* num, double* mean, int* samp, int* len, int* popsize, unsigned int* rng){
    bool input = false; 
    bool is_final = false; 
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        //cout << i << " is "<<argv[i] << "\n"; 
        if ((arg == "-h") || (arg == "--help")) {
            if (rank == 0){
                show_usage(argv[0]);
            }
            is_final = true; 
        } else if ((arg == "-i") || (arg == "--infile")) {
            string temp = argv[i+1];
            input = true; 
            (*in) = temp;  
        } else if ((arg == "-n") || (arg == "--popsize")) {
                (*popsize) = stoi(argv[i+1]); 
        }else if ((arg == "-o") || (arg == "--out")) {
            string temp = argv[i+1];
            input = true; 
            (*out) = temp; 
        }else if ((arg == "-m") || (arg == "--mu")) {
                (*mean) = stod(argv[i+1]); 
        }else if ((arg == "-s") || (arg == "--seed")) { 
                (*rng) = stoi(argv[i+1]); 
        }else if ((arg == "-l") || (arg == "--loci")) {
                (*num) = stoi(argv[i+1]); 
        }else if ((arg == "-b") || (arg == "--length")) {
                (*len) = stoi(argv[i+1]); 
        }else if ((arg == "-f") || (arg == "--samp")) {
                (*samp) = stoi(argv[i+1]); 
        }
    }
    if (input == false){ 
        if (rank==0){
            cerr << "\nKilled: Input must be specified!\n\n"; 
            show_usage(argv[0]);
        }
        is_final = true; 
    }
    if (is_final == true){
        MPI_Finalize(); 
        exit(0); 
    }
}

void show_usage(std::string name){
    std::cerr << "Program: MCpopsim v 0.1\n\nAuthor: Tyler K. Chafin, University of Arkansas\n\nLast modified: July 27, 2015"
              << "\n\nUsage: " << name << " -i </path/to/input> <optional arguments>\n"
              << "\nMandatory arguments:\n"
              << "\t-i,--infile\t- Input file path (tab-delimited event table)"
              << "\n\nOptional arguments:\n"
              << "\t-h,--help\t- Show this help message\n" 
              << "\t-n,--popsize\t- Static population size [default=1000]\n"
              << "\t-m,--mu\t\t- Mean mutation rate; lambda for Poisson distribution [default=0.0000001]\n"
              << "\t-s,--seed\t- Integer seed for random number generator [default=time(0)]\n"  
              << "\t-l,--loci\t- Number of loci to simulate [default=10]\n"
              << "\t-b,--length\t- Static sequence length [default=100]\n"
              << "\t-o,--out\t- Output file path and prefix [default=./]\n"
              << "\t-f,--samp\t- Number of individuals to randomly from end state [default=10]\n"
              << std::endl;
}