/* 
 * File:   allele.h
 * Author: Tyler K. Chafin
 *
 * Created on July 16, 2015, 5:43 AM
 */

#ifndef DEME_H
#define	DEME_H

#include <set>
#include <iostream>
#include <algorithm>
#include <memory>
#include "allele_dat.h"
#include <unordered_set>
#include <map>
//#include "locus.h"


class deme {
public:

    deme();
    deme(int len,int n,double locusmu,int id);
    deme(const deme& orig);
    deme(const std::shared_ptr<deme>* orig,double bottleneck,int id,int src);
    deme(const std::shared_ptr<deme>* par1, double pr1, const std::shared_ptr<deme>* par2, double pr2, int new_idx);
    virtual ~deme();
    
    std::map<unsigned int, allele_dat*> allele_set;
    std::unordered_set<int> allele_ids; //For quick lookup of which alleles are present
    int demeid; 
    unsigned int pop2N; 
    unsigned int current2N = 0;
    int bookkeep; 
    
    void initAlleles(int mut, int len, std::map<unsigned int,allele_dat*> &set_addr);
    //void initAlleles(int);
    void testPrintAlleles() const;
    void printAlleleIDs() const;
    int MCpermutePop(int idx, double pool_size=0.0);
    void correctPopSize(double current, double target);
    int mutate(const double locusmu,int idx);
    void weightedRandomIt(std::map<unsigned int,allele_dat*>::iterator *i);
    void mergeAlleles(const std::shared_ptr<deme>* orig, double ratio, int &total);
    std::string printWeightedRand();
    std::string writeAlleles(int loc) const;
    
    static double getUniform(); 
    
private:
    int seqlength; 
    double mu; 
    static const int INITMUT; 
    unsigned int allele_index=0; 
    unsigned int set_contents=0; //How many nodes in set?
    unsigned int popsize; 

    static const int INITMAX;
    static const int INITMIN; 
    static const double POPSIZEVAR; 
    
};

#endif	/* ALLELE_H */

