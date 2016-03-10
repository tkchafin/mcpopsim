/* 
 * File:   allele.cpp
 * Author: Tyler K. Chafin
 * 
 * Created on July 16, 2015, 5:43 AM
 */

#include "deme.h"
#include "locus.h"
#include <memory>
#include "allele_dat.h"
#include <cstdlib> 
#include <cmath>
#include <algorithm> //advance() set iterator
#include <unordered_set>
#include <map>
#include <climits>


const int deme::INITMUT = 2; //Initial mutations for seed alleles 
const int deme::INITMAX = 10; //Max initial alleles (not including seed allele)
const int deme::INITMIN = 1; //Min initial alleles (not including seed allele)
const double deme::POPSIZEVAR = 0.10; //How much pop size can randomly vary; 0.10 represents 10%

//bool operator< (allele_dat &left, allele_dat &right){
//    return left.key < right.key;
//}

deme::deme() { 
    seqlength = locus::DEFLEN; 
    popsize = locus::DEFPOPSIZE; 
    pop2N = popsize *2; 
    initAlleles(INITMUT, seqlength, allele_set);    
}

//Initializes population 0
deme::deme(int len,int n,double locusmu,int id){
    if (len < locus::MINLEN)
        len = locus::MINLEN; 
    if (len > locus::MAXLEN)
        len = locus::MAXLEN; 
    seqlength = len; 
    popsize = n; 
    pop2N = n*2; 
    mu = locusmu; 
    demeid=id; 
    initAlleles(INITMUT, seqlength, allele_set); 
}

deme::deme(const deme& orig){
    //Fill in later
}


//Copy constructor
deme::deme(const std::shared_ptr<deme>* orig, double bottleneck, int id,int src) {
    seqlength = (*orig)->seqlength;
    popsize = (*orig)->popsize; 
    pop2N = popsize*2; 
    mu = (*orig)->mu; 
    demeid=id; 
    allele_index = (*orig)->allele_index; 
    set_contents = (*orig)->set_contents; 
    
    //Copy alleles from source
    std::map<unsigned int,allele_dat*>::iterator it;//Create iterator variable 
    int k, f;
    char* temp; 

    for (it=(*orig)->allele_set.begin(); it != (*orig)->allele_set.end();++it){
        k = it->first;
        f = it->second->freq; 
        char* temp = it->second->seq; 
        temp = allele_dat::mallocSeq(seqlength,temp);
        allele_set[k] = (new allele_dat(k,f,seqlength,src,temp));
    }
    
}

//"Merging" copy constructor
deme::deme(const std::shared_ptr<deme>* par1, double pr1, const std::shared_ptr<deme>* par2, double pr2, int new_idx){
    /*Uses parameters from 1st parent*/
    seqlength = (*par1)->seqlength;
    popsize = (*par1)->popsize; 
    pop2N = popsize*2; 
    mu = (*par1)->mu; 
    demeid=new_idx; 
    int total = 0; 
    
    /*Sample each parental population, randomly choosing 
     * to keep or pass on alleles of frequency 1/2N*/
    mergeAlleles(par1, pr1, total);
    mergeAlleles(par2, pr2, total); 
    //correctPopSize(total,pop2N); 
    //this->testPrintAlleles(); 
}

void deme::mergeAlleles(const std::shared_ptr<deme>* orig, double ratio, int &total){
    std::map<unsigned int,allele_dat*>::iterator it;//Create iterator variable 
    int k, f,s;
    char* temp; 
    for (it=(*orig)->allele_set.begin(); it != (*orig)->allele_set.end();++it){
        k = it->first;
        s= it->second->source; 
        //std::cout << s << "\n"; 
        f = (double(it->second->freq*ratio)+0.5); 
        char* temp = it->second->seq;
        if (f >= 1){
            if (this->allele_set.count(k) != 0){
                allele_set[k]->freq+= f; 
            }else{
                temp = allele_dat::mallocSeq(seqlength,temp);
                allele_set[k] = (new allele_dat(k,f,seqlength,s,temp));
                this->set_contents++; 
            }
            total += f; 
        }
    }  
}

deme::~deme() {
    std::map<unsigned int,allele_dat*>::iterator it;//Create iterator variable
    for (it=allele_set.begin(); it != allele_set.end(); ++it){ 
        //std::cout <<"Address " << (it->second);
        //std::cout << " Erasing allele "<< it->second->key<<std::endl;
        delete (it->second); 
    }
}


void deme::initAlleles(int mut, int len, std::map<unsigned int,allele_dat*> &set_addr){
    
    //Store a random seed sequence
    char* temp = allele_dat::randomSeq(len); 
    
    //Frequency of seed allele
    int f = rand()%pop2N; 
    int all = pop2N - f;  
    int idx=0; 
    
    //Insert seed allele
    set_addr[idx]= new allele_dat(idx,f,len,0,temp);
    allele_ids.insert(idx);
    this->set_contents++; 
    
    //Initialize with random number of additional alleles 
    int limit = rand()%INITMAX + INITMIN;
    for (int i=1; i<=limit; i++){
        idx++; 
        bookkeep++; 
        this->set_contents++; 
        if (i==limit){
            f=all; 
            all = 0; 
        }else{
            f = rand()%all;
            all -= f; 
        }
        temp = allele_dat::fixedMutateSeq(mut,len,temp);
        //std::cout << "Making allele " << &temp <<"  "<< temp->getKey() << std::endl;
        set_addr[idx] = new allele_dat(idx,f,len,0,temp); 
        allele_ids.insert(idx);
    }
}

void deme::testPrintAlleles() const{
    //std::map<int,allele_dat*>::iterator it;//Create iterator variable
    for (auto it=allele_set.begin(); it != allele_set.end(); ++it){ 
        std::cout <<"Allele " << it->first << "("
        <<it->second->source<< ") " 
        <<" ; Freq " << it->second->freq/pop2N << " - " 
        << it->second->freq << "; ";
        it->second->allele_dat::printSeq();             
    }
}

void deme::printAlleleIDs() const{
    for (auto it=allele_ids.begin(); it != allele_ids.end(); ++it){ 
        std::cout <<(*it)<<" ";                    
    }
    std::cout << std::endl;
}

int deme::MCpermutePop(int idx, double bot){
    std::map<unsigned int,allele_dat*>::iterator it;//Create iterator variable 
    double p,q,x,newp; 
    int newfreq; 
    int alt; 
    int pool_size; 
    this->current2N=0; 
    if (locus::doubleCompare(bot, 0.0))
        pool_size = this->pop2N;
    else if (bot <= locus::MINBOT)
        pool_size = pop2N*locus::MINBOT; 
    else
        pool_size = bot*pop2N; 
    for (it=allele_set.begin(); it != allele_set.end();){ 
        p = it->second->freq /pop2N ; 
        q = 1-p; 
        x = getUniform(); 
        if (p<=0)
            newp = 0; 
        else
            newp = (p+(((2*x)-1)*(sqrt((3*p*q)/pool_size))));
        newfreq = newp * pop2N; 
        this->current2N += newfreq; 
     
        if (newfreq <= 0){ 
            allele_ids.erase(it->first);
            delete (it->second);
            it = allele_set.erase(it);
            
            this->set_contents--;  
            
        }else if (newfreq >= pop2N){
            it->second->freq = pop2N; 
            ++it;
        }else{
            it->second->freq = newfreq; 
            ++it;
        }
        //std::cout << "p= " << p*pop2N << " newp= " << newfreq << " Total = " << current2N<< "\n"; 
    }
    
    alt = deme::mutate(mu,idx); 
    
    //If popsize has drifted further than allowed, correct
    if ((current2N >= double(pop2N + double(pop2N*POPSIZEVAR))) 
        || (current2N <= double(pop2N - double(pop2N*POPSIZEVAR)))){
            //std::cout << "Correcting 2N of "<< current2N << "\n";
            correctPopSize(current2N,pop2N);
            current2N = pop2N; 
    }
   // std::cout << "\nNext generation\n\n";
    return alt; 
}

void deme::correctPopSize(double current, double target){ 
    double scale = ((double(target))/double(current));
    //std::cout << "Current: "<<current<<" Target: "<<target <<"\n";
    //std::cout << "\nScale: " << current/target << std::endl; 
    
    //std::cout << "Scale is " << scale << std::endl; 
    std::map<unsigned int,allele_dat*>::iterator it;//Create iterator variable 
    int newfreq; 
    for (it=allele_set.begin(); it != allele_set.end();){ 
        newfreq = (double(it->second->freq * scale)+0.5); //+0.5 because c++ always truncates to 0
 
        if (newfreq <= 0){
            allele_ids.erase(it->first);
            delete (it->second); 
            it = allele_set.erase(it);
            this->set_contents--; 
        }else{ 
            it->second->freq = newfreq; 
            ++it;
        }
    }
}


int deme::mutate(const double locusmu, int idx){
    int alt=0; 
    double x = getUniform(); 
    double muprob = (locusmu*pop2N*seqlength);
    int pre;  
    //std::cout << "From allele::mutate " << locusmu << "\n";
    /*Create iterator variable*/
    std::map<unsigned int,allele_dat*>::iterator it = allele_set.begin();
    
    /*While mutating (decrement local mu for each event
     because each mutation is independent thus their combined 
     probability is the sum*/
    while(x <= muprob){

        it = allele_set.begin();

        /*Capture random iterator*/
        do{
            this->weightedRandomIt(&it);
        }while((it->second->freq) < 1 );
        
        //std::cout<<"Before: "<<(*it)->freq <<"\n";
        //std::cout << "MUTATE! " << x << " allele " << (*it)->key <<std::endl; 

        /*Grab random source*/
        char* temp = it->second->seq; // </editor-fold>
 
        it->second->freq--; //Decrement source allele frequency
        //std::cout << "Allele " << (*it)->key << " before: ";allele_dat::printSeqDirect(temp,seqlength);
        
        this->set_contents++;
        idx++;
        alt++; 
        
        /*Mutate randomly chosen sequence and insert as new allele_dat*/
        temp = allele_dat::fixedMutateSeq(1, seqlength,temp);
 
        //std::cout << "After: "; allele_dat::printSeqDirect(temp,seqlength);
        allele_set[idx] = new allele_dat(idx,1,seqlength,demeid,temp);
        allele_ids.insert(idx);
        
        /*Check if we have overflowed UINT_MAX*/
        if ((idx >= UINT_MAX) && (locus::MUT_OVERFLOW == false))
            locus::MUT_OVERFLOW = true; 
        /*Capture some new values*/
        muprob -= x; 
        x = getUniform();
        
    }
    return alt; 
}

void deme::weightedRandomIt(std::map<unsigned int,allele_dat*>::iterator *i){
    std::map<unsigned int,allele_dat*>::iterator it;
    unsigned int total = current2N; 
    unsigned int picked = rand()%total; 
    for (it=allele_set.begin();it!=allele_set.end();++it){
        total -= it->second->freq; 
        if (total <= picked){
            (*i)=it;
            break;
        }
    }
}
std::string deme::printWeightedRand(){
    std::string ret = ""; 
    auto it = allele_set.begin();
    do{
        this->weightedRandomIt(&it);
    }while((it->second->freq) < 1 );
    ret += "_Allele." + std::to_string(it->first);
    ret += "(" + std::to_string(it->second->source) + ")";
    ret += "_Freq." + std::to_string((it->second->freq)/pop2N) + "\n";
    for (int i=0; i<seqlength; i++){
        ret += it->second->seq[i];
    }
    return ret;   
}


std::string deme::writeAlleles(int loc) const{
    std::string ret = ""; 
    for (auto it=allele_set.begin(); it != allele_set.end(); ++it){ 
        ret+= ">Locus." + std::to_string(loc)+"_Population."+std::to_string(demeid); 
        ret += "_Allele." + std::to_string(it->first);
        ret += "(" + std::to_string(it->second->source) + ")";
        ret += "_Freq." + std::to_string((it->second->freq)/pop2N) + "\n";
        for (int i=0; i<seqlength; i++){
            ret += it->second->seq[i];
        }
        ret += "\n";                       
    }
    return ret; 
}

/*Samples uniform(0,1) dist*/
double deme::getUniform(){
    return rand()/(double)RAND_MAX;
}

