/* 
 * File:   locus.h
 * Author: Tyler K. Chafin
 *
 * Created on July 16, 2015, 5:38 AM
 */

#ifndef LOCUS_H
#define	LOCUS_H

#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>
#include "deme.h" 
#include "allele_dat.h"
#include <memory>
#include <unordered_map>

class locus {
public:
    /*Constructors etc*/
    locus();
    locus(double,int,int,double,int,unsigned int s); 
    locus(const locus& orig);
    virtual ~locus();
    
    /*Utilities*/
    void setLength(const int len);
    void setMu(const double rate);
    void initPopAlleles();
    void testPrintPop()const; 
    void cladogenesis(const int seedpop,int id,double bot); 
    void unidirGeneFlow(int source,int sink,double rate);
    void intOverflowCheck()const;
    void finalize();
    void hybridogenesis(int par1, double pr1, int par2, double pr2, int new_idx);
    void incrementGen();
    int getGen()const;
    int setGen(int x);
    void eventParser();
    void MCpermuteAllPops();
    void parseMigMatrix();
    void writeRandomFasta(); 
    void writeCompleteFasta();
    
    /*Static functions*/
    static bool doubleCompare(double x, double y);
    static double logPoissonMu(const double lambda,unsigned int seed); 
    
    /*Static vars*/
    static const double MAXMU; 
    static const double DEFMU;
    static const int MAXLEN; 
    static const int MINLEN;
    static const int DEFLEN;
    static const double MINBOT;
    static const int MINPOPSIZE; 
    static const int DEFPOPSIZE;
    static const double DEFBOT; 
    static bool MUT_OVERFLOW; 
    static int NUM_SAMPLES; 
    static std::string FILEBASE; 
    static const double MINMU; 
    /*Public objects*/
    int getLength() const; 
    double getMu() const; 
    unsigned int RNGseed; 
    
    std::vector<std::shared_ptr<deme>> pops; 
    std::unordered_map<int,int> demeids; //key=demeid; value=pops index

    
private:
    double mu; 
    int seqlength; 
    int popsize; 
    int locus_id; 
    int numpops=0; 
    double bottleneck; 
    int locusid; 
    int gen; 
    int allele_index=0; 
    


 
     
};

#endif	/* LOCUS_H */

