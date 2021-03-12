#include <iostream>
#include <iomanip>
#include <cmath>
#include "Crank_Nicolson.h"
#include <chrono>

int main()
{
    
    std::cout << std::setprecision(10);
    
    /*
    double X=100.,T=1.,r=0.06,sigma=2, delta=0;
    int iMax = 30, jMax = 1000;
    // Crank_Nicolson_Call(100., 95, 3., 0.0043, 17.1, 0.0043, 200, 50000);
    Crank_Nicolson_Put  (15., 12, 3., 0.06, 17, 0, 200, 50000);
    */
    
    // Crank_Nicolson_Bond(100, 2000);
    
    // analytic_estimate();
    
    // get_abs_err();
    
    // get_convergange_rate();
    
    // get_convergange_rate_extrap();
    
    // testing_Smax();
    
    // testing_Smax_analytic();
    
    
    for (int n{1}; n<100; n++)
    {
        double dS = 5./n;
        double j1 = 115/dS, j2 = 95/dS;
        std::cout << "n: " << n << std::endl;
        std::cout << "j1: " << j1 << " j2: " << j2 << std::endl;
    }
    
    return 0;
}


