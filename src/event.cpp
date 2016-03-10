/* 
 * File:   event.cpp
 * Author: tyler
 * 
 * Created on July 26, 2015, 2:09 AM
 */

#include "event.h"
#include <fstream>
#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <unordered_set>

int event::END; 
std::vector<std::vector<std::string>> event::event_table;
std::unordered_set<int> event::important_dates; 
double** event::mig_matrix; 
int event::matrix_dim = 0; 

event::event() {
    filename = "event_table.tsv";
    END = 1; 
}

event::event(std::string infile) {
    filename = infile; 
    readEventTable(); 
}

event::event(const event& orig) {
}

event::~event() {
    for (int i=0; i<matrix_dim; i++){
        delete[] mig_matrix[i]; 
    }
    delete[] mig_matrix;
}

void event::readEventTable(){ 
    std::ifstream filehandle; 
    filehandle.open(filename);
    std::string line; 
    std::string partial; 
    std::vector<std::vector<std::string>*> lines; 
    std::vector<std::string> tokens; 
    
    if(filehandle.is_open()){
        while(std::getline(filehandle, line)){
            
            std::vector<std::string> tokens;
            int count = 0; 
            //tokens = new std::vector<std::string>; 
            
            std::istringstream iss(line); 
            std::string token; 
            
            if ((line[0] != ' ') && (line[0] != '#') && (line[0] != '/') && (line[0] != '\n')){
                while (std::getline(iss,  token, '\t')){
                        tokens.push_back(token); 
                }
                if (tokens.size() > 0){
                    //std::cout << stoi(tokens[1])<< "\n";
                    if (tokens[2] == "end"){
                        event::END = std::stoi(tokens[1]);                      
                    }else if (tokens[0] == "m"){
                        std::vector<std::string> temp; 
                        temp = tokens;
                        temp.erase(temp.begin()+2);
                        event::important_dates.insert(std::stoi(temp[1]));
                        event::event_table.push_back(temp);
                        temp = tokens;
                        temp.erase(temp.begin()+1);
                        temp[4] = "0"; 
                        event::important_dates.insert(std::stoi(temp[1]));
                        event::event_table.push_back(temp);
                    }else{
                        event::important_dates.insert(std::stoi(tokens[1]));
                        event::event_table.push_back(tokens);
                    }
                    if((tokens[0] == "t") || (tokens[0] == "h"))
                        event::matrix_dim++; 
                }
            }
        } 
        allocMatrix(event::matrix_dim);
    }
}

void event::testPrintTable()const{
    std::cout << "Loading demographic events:\n"; 
    for (int i=0; i < event::event_table.size(); ++i){
        for (int j=0; j< event::event_table[i].size(); ++j)
            std::cout << std::setw(10) << event::event_table[i][j];
        std::cout << std::endl; 
    }
    std::cout << "Simulation will end at time " << END << "\n\n";
}

void event::printMigMatrix(){
    std::cout << std::endl;
    for (int i=0; i < event::matrix_dim; ++i){
        for (int j=0; j< event::matrix_dim; ++j)
            std::cout << std::setw(10) << event::mig_matrix[i][j];
        std::cout << std::endl; 
    }
}

void event::allocMatrix(int dim){
    mig_matrix = new double*[dim];
    for (int i=0; i<dim; ++i){
        mig_matrix[i]= new double[dim]; 
        for (int j=0; j<dim; ++j){
            mig_matrix[i][j] = 0.0; 
        }
    }    
}
