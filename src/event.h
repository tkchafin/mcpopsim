/* 
 * File:   event.h
 * Author: Tyler K. Chafin
 *
 * Created on July 26, 2015, 2:09 AM
 */

#ifndef EVENT_H
#define	EVENT_H

#include <string>
#include <vector>
#include <map>
#include <unordered_set>

class event{
public:
    event();
    event(std::string);
    event(const event& orig);
    virtual ~event();
    void readEventTable();
    void testPrintTable()const;
    void allocMatrix(int);
    
    static void printMigMatrix(); 
    
    static std::vector<std::vector<std::string>> event_table;
    static std::unordered_set<int> important_dates; 
    static double** mig_matrix;  
    static int matrix_dim; 
    static int END; 
    
private:
    std::string filename;
    

};

#endif	/* EVENT_H */

