/* 
 * File:   locus.cpp
 * Author: Tyler K. Chafin
 * 
 * Created on July 16, 2015, 5:38 AM
 */

#include "locus.h"
#include <random> //poisson
#include <cmath> //log10 + abs(double)
#include <iostream>
#include <iomanip> 
#include <ctime> //time
#include "deme.h"
#include "allele_dat.h"
#include <memory> //shared_ptr
#include <unordered_map>
#include <cassert>
#include <map>
#include <climits>
#include "event.h"
#include <string>
#include <fstream>
#include <stdlib.h>


//Define constants
const int locus::MAXLEN = 10000;
const int locus::MINLEN = 10; 
const double locus::MAXMU = 0.005; 
const double locus::MINMU = 0.0;
const int locus::MINPOPSIZE = 10; 
const int locus::DEFPOPSIZE = 10000; 
const double locus::DEFMU = 0.00000001;
const int locus::DEFLEN = 100; 
const double locus::DEFBOT = 0.5; 
bool locus::MUT_OVERFLOW = false;
const double locus::MINBOT = 0.1;
int locus::NUM_SAMPLES = 10; 
std::string locus::FILEBASE = "";
 

//Default constructor
locus::locus() {
    //Set defaults
    mu = locus::DEFMU;
    seqlength = locus::DEFLEN; 
    popsize = locus::DEFPOPSIZE;
    bottleneck = locus::DEFBOT; 
    gen=0; 
}

//Custom constructor
locus::locus(double rate,int len,int N,double bot, int id, unsigned int s){
    if (rate > 0.0){ 
        mu = rate; 
    }else{ 
        mu = 0.0; 
    } 
    if (rate > locus::MAXMU){
        mu = locus::MAXMU; 
    } 
    if (len <= locus::MAXLEN){
        seqlength = len; 
    }else{
        seqlength = locus::MAXLEN;
    }
    if (len < locus::MINLEN){ 
        seqlength = locus::MINLEN; 
    } 
    if (N >= MINPOPSIZE)
        popsize = N;
    else 
        popsize = MINPOPSIZE;  
    if (bot >=1.0)
        bot = 1.0;
    else if (bot <= 0.0)
        bot = 0.1; 
    bottleneck = bot; 
    locusid=id; 
    //allele alleles = new allele(10); 
    //allele alleles(10);
    //This isn't entirely working, but the version in main has no memory leaks
   
    pops.push_back(std::make_shared<deme>(seqlength,popsize,mu,0));
    demeids[0]=0;
    numpops++;
    gen=0; 
    allele_index += pops[0]->bookkeep; 
    //allele::allele
    RNGseed = s; 
}


locus::locus(const locus& orig) {
}

locus::~locus() {

    //std::vector<allele*>::iterator it;//Create iterator variable
    //for (it=pops.begin(); it != pops.end(); ++it) ;
    //pops[0]->~allele();
    //delete it; 
      //(*it)->allele::~allele();
       //it->~locus();    
}

double locus::getMu() const{ 
    return mu; 
}

int locus::getLength() const{
    return seqlength; 
}

void locus::setLength(int len){
    if (len <= locus::MAXLEN){
        seqlength = len; 
    }else{
        seqlength = locus::MAXLEN;
    }
    if (len < locus::MINLEN){ 
        seqlength = locus::MINLEN; 
    } 
}

void locus::setMu(const double rate){
    if (rate > 0.0){ 
        mu = rate; 
    }else{ 
        mu = 0.0; 
    } 
    if (rate > locus::MAXMU){
        mu = locus::MAXMU; 
    } 

}


//Send mean mu adjusted for seqlength and pop size
double locus::logPoissonMu(const double lambda, unsigned int seed){

    typedef std::mt19937 G;
    typedef std::poisson_distribution<> D;

    G generator(seed);  // seed if you want with integral argument
    D mudist(-log10(lambda));
    double tempmu = pow(10,-(mudist(generator))); 
    if (tempmu > locus::MAXMU)
        return logPoissonMu(lambda,rand());
    else if (tempmu < locus::MINMU)
        return logPoissonMu(lambda,rand()); 
    else
        return tempmu; 
}
//double locus::gammaMu(const double k, const double theta){

void locus::testPrintPop()const{ 
    for (auto i=demeids.begin(); i != demeids.end(); ++i){ 
        std:: cout << "Population " << i->first << std::endl;  
        pops[i->second]->testPrintAlleles(); 
    }
}

void locus::cladogenesis(const int seedpop,int id,double bot){
    
    int source_idx = demeids[seedpop];
    pops.push_back(std::make_shared<deme>(&pops[source_idx],bottleneck,id,seedpop));
    numpops++; 
    demeids[id]=numpops-1;
    pops[numpops-1]->MCpermutePop((pops[source_idx]->pop2N)*bot);
    
}

void locus::unidirGeneFlow(int source,int sink, double rate){
    
    double x = rand()/(double)RAND_MAX; 
    /*Capture which populations are intended for gene flow*/
    //assert(demeids.count(source));
    //assert(demeids.count(sink)); 
    int source_idx = demeids[source];
    int sink_idx = demeids[sink];
    unsigned int k;
    char* temp;
    int orig;
    std::map<unsigned int,allele_dat*>::iterator it;

    while(x <= rate){
        /*Exchange two alleles b/c diploid*/
        for (int j=0; j<2;j++){
            it = pops[source]->allele_set.begin();
    
            do{
                pops[source_idx]->weightedRandomIt(&it);
            } while((it->second->freq) < 1 );
    
            k = it->first; 
            temp = it->second->seq; 
            orig = it->second->source; 
            
            it ->second->freq--; //decrement source frequency

            if ((pops[sink_idx]->allele_set.count(k)) != 0){ //If allele exists in sink
                it = (pops[sink_idx]->allele_set.find(k));  
                it->second->freq; //increment sink frequency
            }else{
                //std::cout << " - Adding new allele!\n"; 
                temp = allele_dat::mallocSeq(seqlength,temp);
                pops[sink_idx]->allele_set[k] = new allele_dat(k,1,seqlength,orig,temp);
                pops[sink]->allele_ids.insert(k);
            }
        }
        /*Capture some new values*/
        rate -= x; 
        x = deme::getUniform();
    }
}

void locus::finalize(){ 
    for (int i=0; i < numpops; i++){ 
        pops[i]->correctPopSize(pops[i]->current2N,pops[i]->pop2N);
    }
    intOverflowCheck();
    writeRandomFasta();
    writeCompleteFasta(); 
}

void locus::intOverflowCheck()const{
    if (locus::MUT_OVERFLOW == true)
        std::cerr << "WARNING: Allele numbers overflowed UINT_MAX. Common allele indices may not represent identical sequences.\n";
}
 
void locus::hybridogenesis(int par1, double pr1, int par2, double pr2, int new_idx){
    int par1_idx = demeids[par1];
    int par2_idx = demeids[par2]; 
    double sum = pr1+pr2; 
    
    if (!doubleCompare(sum,1.0)){
        /*If they are NOT equal, correct values*/
        double scale = sum/1.0; 
        pr1 = pr1/scale; 
        pr2 = pr2/scale;
    }
    pops.push_back(std::make_shared<deme>(&pops[par1_idx],pr1, &pops[par2_idx],pr2, new_idx));
    numpops++; 
    demeids[new_idx]=numpops-1;
    //pops[numpops-1]->MCpermutePop((pops[source_idx]->pop2N)*bottleneck);
}

bool locus::doubleCompare(double x, double y){
    /*Note that this function is asymmetric;
    e.g. doubleCompare(x,y) MAY BE != doubleCompare(y,x)*/
    const double epsilon = 0.0000001; 
    return std::abs(x - y) <= epsilon * std::abs(x);
    /* see Knuth section 4.2.2 pages 217-218
     * Code modified from: 
     * https://isocpp.org/wiki/faq/newbie#floating-point-arith
     * */
}

void locus::incrementGen(){
    gen++; 
}

int locus::getGen()const{
    return gen; 
}

int locus::setGen(int x){
    gen = x; 
}

void locus::eventParser(){
    int g = this->getGen(); 
    int pop1, pop2, idx; 
    double rate1, rate2; 
    char type; 
    if ((event::important_dates.count(g))!=0){
        //std::cout << "Event at time " << g << std::endl; 
        for (int j=0; j<event::event_table.size(); j++){
            if (std::stoi(event::event_table[j][1]) == g){
                type = event::event_table[j][0][0]; 
               // std::cout << "type " << type << " ";
                
                if (type == 't'){
                    pop1 = std::stoi(event::event_table[j][2]); //Need to throw exception if invalid source pop given!!
                    pop2 = std::stoi(event::event_table[j][3]);
                    rate1 = std::stod(event::event_table[j][4]); //Even better we could check all these things when reading in file!??
                    cladogenesis(pop1,pop2,rate1);
                }else if (type== 'h'){
                    pop1 = std::stoi(event::event_table[j][2]);
                    rate1 = std::stod(event::event_table[j][3]);
                    pop2 = std::stoi(event::event_table[j][4]);
                    rate2 = std::stod(event::event_table[j][5]);
                    idx = std::stoi(event::event_table[j][6]);
                    //std::cout << pop1<< " " << rate1 << " " << pop2 <<" "<<rate2<<" "<<idx<<"\n";
                    hybridogenesis(pop1,rate1,pop2,rate2,idx);  //Need to throw exception if attempted population overwrite!
                }else if (type == 'm'){
                    pop1 = std::stoi(event::event_table[j][2]);
                    pop2 = std::stoi(event::event_table[j][3]);
                    rate1 = std::stod(event::event_table[j][4]);
                    event::mig_matrix[pop1][pop2] = rate1; 
                    //event::printMigMatrix(); 
                }
            }
        }
        //std::cout << std::endl;
    }
}

void locus::MCpermuteAllPops(){
    int x; 
    for (int i = 0; i< pops.size(); i++){
        x = pops[i]->MCpermutePop(allele_index,1.0);    
        allele_index += x; 
    }
}

void locus::parseMigMatrix(){ 
    double rate; 
    for (int source=0; source < event::matrix_dim; ++source){
        for (int sink=0; sink< event::matrix_dim; ++sink){
            rate = event::mig_matrix[source][sink]; 
            if (!doubleCompare(rate,0.0)){                
               this->unidirGeneFlow(source,sink,rate);
            }
        }
    }
}

void locus::writeRandomFasta(){ 
    std::ofstream outfile; 
    std::string loc = std::to_string(locusid); 
    std::string outname = locus::FILEBASE + "locus_" + loc;
    std::map<unsigned int,allele_dat*>::iterator map_it; 
    outname += ".fasta";
    outfile.open(outname);
    if (outfile.is_open()){
        for (auto it=demeids.begin(); it != demeids.end(); ++it){ 
            for (int i=0; i<NUM_SAMPLES; i++){
                //do{
                //int id = demeids[it->first];   
                //}while(((*map_it).second->freq) < 1 );
                outfile << ">Locus." <<locusid<< "_Population." << it->first << "_Individual."<< i; 
                outfile << pops[it->second]->printWeightedRand(); 
                outfile << "\n"; 
            }
        }
    }else{
        std::cerr << "Error! Could not open file " << outname << "!\n";
    }
}

void locus::writeCompleteFasta(){
    std::ofstream outfile; 
    std::string loc = std::to_string(locusid); 
    std::string outname = locus::FILEBASE + "locus_" + loc;
    std::map<unsigned int,allele_dat*>::iterator map_it; 
    outname += ".alleles";
    outfile.open(outname);
    if (outfile.is_open()){
        for (auto it=demeids.begin(); it != demeids.end(); ++it){   
            outfile << pops[it->second]->writeAlleles(locusid);  
        }
    }else{
        std::cerr << "Error! Could not open file " << outname << "!\n";
    }
}
