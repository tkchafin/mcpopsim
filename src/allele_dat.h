/* 
 * File:   allele_dat.h
 * Author: Tyler K. Chafin
 *
 * Created on July 16, 2015, 5:41 AM
 */

#ifndef ALLELE_DAT_H
#define	ALLELE_DAT_H


class allele_dat {
public:
    allele_dat();
    allele_dat(int,int,int);
    allele_dat(int,int,int,int,char*);
    allele_dat(const allele_dat& orig);
    virtual ~allele_dat();
    unsigned int key; 
    int source; 
    double freq; 
    char* seq; 
    static char* randomSeq(int len);
    static char randomNuc(); 
    bool operator<(const allele_dat &other)const;
    void construct(int,double,int);
    int getKey()const; 
    double getFreq() const; 
    void printSeq() const; 
    static char* mallocSeq(const int len, char* &ancestor);
    static char* fixedMutateSeq(const int num_mut, const int len, char* &ancestor);
    static void printSeqDirect(char* print, int len);
    
private:
    int seqlength; 
    static const char NUCS[4];
};

#endif	/* ALLELE_DAT_H */

