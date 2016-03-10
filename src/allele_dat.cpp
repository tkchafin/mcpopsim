/* 
 * File:   allele_dat.cpp
 * Author: Tyler K. Chafin
 * 
 * Created on July 16, 2015, 5:41 AM
 */

#include "allele_dat.h"
#include "locus.h"
#include <iostream>

const char allele_dat::NUCS[4] = {'A', 'G', 'T', 'C'}; 

//Default constructor
allele_dat::allele_dat() {
    key = 0; 
    freq = 0; 
    seqlength = locus::DEFLEN; 
    seq = randomSeq(seqlength);
}

//Constructor, building random sequence
allele_dat::allele_dat(int k, int f, int len) {
    key = k; 
    freq = f; 
    seqlength = len; 
    seq = randomSeq(seqlength);
}

//Constructor where sequence is provided
allele_dat::allele_dat(int k, int f, int len, int src, char* s) {
    key = k; 
    freq = f; 
    seqlength = len;
    source = src;
    seq = s; 
}


allele_dat::allele_dat(const allele_dat& orig) {
}

allele_dat::~allele_dat() {
    delete[] seq; //Look into using std::memory instead for shr ptr
}

char* allele_dat::randomSeq(int len){
    char* temp = new char [len]; 
    for (int i=0; i<len;i++){
        temp[i] = randomNuc(); 
    }
    return temp; 
}

char* allele_dat::fixedMutateSeq(const int num_mut, const int len, char* &ancestor){
    char* temp = new char[len]; //Declare new char*
    //Copy ancestor sequence
    for (int i=0; i< len;i++){
        temp[i] = ancestor[i]; 
    }
    //Randomly permute it by fixed number of substitutions
    int pos; 
    for (int i=0; i<num_mut; i++){
        pos = rand()% len;
        char newNuc; 
        //Capture a random nuc (jukes cantor)
        do{
            newNuc = randomNuc();
        }while (temp[pos] == newNuc);
        //Until we get a new one (don't allow mutation from e.g. A->A)
        
        temp[pos] = newNuc; 
    }
    return temp; 
}

char* allele_dat::mallocSeq(const int len, char* &ancestor){
    char* temp = new char[len]; //Declare new char*
    //Copy ancestor sequence
    for (int i=0; i< len;i++){
        temp[i] = ancestor[i]; 
    }
    return temp; 
}

//assumes jukes-cantor (1969),maybe add other models later
char allele_dat::randomNuc(){
    int randindex = rand()% 4; 
    //return randindex; 
    return NUCS[randindex];
}

int allele_dat::getKey()const{
    return this->key; 
}

double allele_dat::getFreq() const{
    return this->freq;
}

void allele_dat::printSeq() const{
    for (int i=0; i< this->seqlength; ++i){
        std::cout << this->seq[i]; 
    }
    std::cout << std::endl;
}


bool allele_dat::operator< (const allele_dat &other)const{
    return key < other.key;
}

//Construct an allele_dat with a random sequence
void allele_dat::construct(int k, double f, int len) {
    key = k; 
    freq = f; 
    seqlength = len; 
    seq = randomSeq(seqlength);
}

void allele_dat::printSeqDirect(char* print, int len){
    for (int i=0; i< len; ++i){
        std::cout << print[i]; 
    }
    std::cout << std::endl;
}
