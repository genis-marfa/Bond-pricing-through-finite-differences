#include <iostream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <chrono>
#include <vector>
#include <cassert>
#include <fstream>

// Define parameters globally:
double T = 3., F = 190., R = 2., r = 0.0043, K = 0.08333333333, mu = 0.0076, X = 95.24, C = 0.408, alpha = 0.03, beta = 0.142, sigma = 17.1, t_0 = 1.2654, Cp = 230;

template <class T>
T max(T n1, T n2)
{
    return (n1 > n2) ? n1 : n2;
}

std::vector<double> bunch_time_points(double t_0, double iMax)
{
    std::vector<double> time_grid_points;
    double dt = T/iMax;
    for (int i{}; i<=iMax; i++)
    {
        if ( t_0 >= (i-0.5)*dt  && t_0 <=  (i+0.5)*dt)
        {
            for (int j=0; j<=100; j++)
            {
                time_grid_points.push_back((i-0.5+j/100.)*dt);
            }
            
        }
        else time_grid_points.push_back(i*dt);
    }
    return time_grid_points;
}

std::vector<double> LU_decompose (std::vector<double> a, std::vector<double> b,
                                  std::vector<double> c, std::vector<double> d_i)
{
    assert(d_i.size() == c.size() && c.size() == b.size() && b.size() == a.size() );
    size_t size = d_i.size();
    
    double b_0 = b[0];
    
    std::vector<double> beta(size);
    beta[0] = b_0;
    for (int j{1}; j< size; j++) beta[j] = b[j] - (a[j]*c[j-1])/beta[j-1];
    double beta_max = beta[size-1];
    
    std::vector<double> D(size);
    D[0] = d_i[0];
    for (int j{1}; j< size; j++) D[j] = d_i[j] - a[j]/beta[j-1] * D[j-1];
    double D_max = D[size-1];
    
    std::vector<double> V_i(size);
    V_i[size-1] = D_max/beta_max;
    for (int j=size-2; j>=0; j--) V_i[j] = 1/beta[j]*(D[j]-c[j]*V_i[j+1]);
    
    return V_i;
}

// A generic lagrange interpolation function
double lagrangeInterpolation(const std::vector<double>& y,const std::vector<double>& x,double x0,unsigned int n)
{
    if(x.size()<n) return lagrangeInterpolation(y,x,x0,x.size());
    if(n==0) throw;
    int nHalf = n/2;
    int jStar;
    double dx=x[1]-x[0];
    if(n%2==0)
        jStar = int((x0 - x[0])/dx) -(nHalf-1);
    else
        jStar = int((x0 - x[0])/dx+0.5)-(nHalf);
    jStar=std::max(0,jStar);
    jStar=std::min(int(x.size()-n),jStar);
    if(n==1)return y[jStar];
    double temp = 0.;
    for(unsigned int i=jStar;i<jStar+n;i++){
        double  int_temp;
        int_temp = y[i];
        for(unsigned int j=jStar;j<jStar+n;j++){
            if(j==i){continue;}
            int_temp *= ( x0 - x[j] )/( x[i] - x[j] );
        }
        temp += int_temp;
    }
    // end of interpolate
    return temp;
}

// CONVERTIBLE BOND:

// Functions for convertible bond:
double theta(double t) {return (1+mu)*X*exp(mu*t);}

// Functions for boundary at S -> inf:
double A(double t)
{
    return R*exp((K+r)*(t-T));
}
double B(double t)
{
    return C/(alpha+r)*(exp(-alpha*t)-exp(-alpha*T+r*(t-T))) + X*R*(exp(r*(t-T)) - exp((K+r)*(t-T)));
}

// Numerical scheme coefficients:
double get_a(int j, int i, double dS, double dt)
{
        return 0.25*(sigma*sigma*pow(j, 2*beta)*pow(dS, 2*(beta-1))-K*(theta(i*dt)/dS-j));
}

double get_b(int j, int i, double dS, double dt)
{
    return -0.5*sigma*sigma*pow(j, 2*beta)*pow(dS, 2*(beta-1)) - 0.5*r - 1./dt;
}

double get_c(int j, int i, double dS, double dt)
{
    return 0.25*(sigma*sigma*pow(j, 2*beta)*pow(dS, 2*(beta-1))+K*(theta(i*dt)/dS-j));
}

double get_d(int j, int i, double dS, double dt, std::vector<double> vOld)
{
    return -get_a(j, i, dS, dt)*vOld[j-1] - (-0.5*sigma*sigma*pow(j, 2*beta)*pow(dS, 2*(beta-1)) - 0.5*r + 1./dt)*vOld[j] - get_c(j,i,dS, dt)*vOld[j+1] - C*exp(-alpha*i*dt);
}

// Implement Crank-Nicolson method for the convertible bond:
double Crank_Nicolson_American_Bond(double S_want, int iMax, int jMax, double S_max)
{
    // Declare and Initialise local variables: dS,dt, Smax:
    double dS = S_max/jMax;
    double dt = T/iMax ;
    // std::cout << "dS: " << dS << " dt: " << dt << std::endl;
    
    // First appply terminal condition V(S,t = T):
    std::vector<double> vNew(jMax+1), S(jMax+1);
    for (int j{}; j<=jMax; j++) S[j] = j*dS;
    for (int j{}; j<=jMax; j++) vNew[j] = max(F, R*S[j]);
    std::vector<double> vOld = vNew;
    
    for (int i = iMax-1; i >= 0; i--)
    {
        std::vector<double> a(jMax+1), b(jMax+1), c(jMax+1), d(jMax+1);
        
        // Set up boundary condition at S = 0:
        a[0] = 0;
        b[0] = -(1./dt+K*theta(i*dt)/dS+r/2);
        c[0] = K*theta(i*dt)/dS;
        d[0] = (-1./dt+r/2)*vOld[0]-C*exp(-alpha*i*dt);
        
        for (int j{1}; j<jMax; j++)
        {
            a[j] = get_a(j, i, dS, dt);
            b[j] = get_b(j, i, dS, dt);
            c[j] = get_c(j, i, dS, dt);
            d[j] = get_d(j, i, dS, dt, vOld);
        }
        
        // Set up boundary condition at S -> inf.
        a[jMax] = 0;
        b[jMax] = 1;
        c[jMax] = 0;
        d[jMax] = R*jMax*dS;
        
        // Apply penalty method:
        double penalty=1.e8, tolerance = 1.e-8; int iterMax = 100, penaltyIt = 0;
        for(int penaltyIt=0; penaltyIt< iterMax; penaltyIt++)
        {
            std::vector<double> bHat(b),dHat(d);
            for(int j=1;j<jMax;j++)
            {
                // turn on penalty if V < RS for conversion:
                if(vNew[j] < S[j]*R)
                {
                    bHat[j] = b[j] - penalty;
                    dHat[j] = d[j] - penalty*S[j]*R;
                }
                
                // If t < t_0, turn on penalty if V > max(Cp, RS) for call Option:
                if (i*dt < t_0 && vNew[j] > max(Cp, R*S[j]) )
                {
                  bHat[j] = b[j] - penalty;
                  dHat[j] = d[j] - penalty*max(Cp, R*S[j]);
                }
                
            }
            // solve matrix equations with LU decompomposition:
            std::vector<double> y = LU_decompose(a,bHat,c,dHat);
            
            // calculate difference from last time
            double error=0.;
            for(int j=0;j<=jMax;j++) error += fabs(vNew[j] - y[j]);
            vNew = y;
            if(error < tolerance) break;
        }
        if (penaltyIt >= iterMax)
            {std::cout << "Failed to converge after " << iterMax << " iterations. " << '\n';}
        vOld = vNew;
    }
    double value = lagrangeInterpolation(vNew,S,S_want,4);
    // std::cout << "V(S=" << S_want << ", t=0) = " << value << std::endl;
    return value;
}


// Implement Crank-Nicolson method for the convertible bond:
double Crank_Nicolson_American_Bond_bunched(double S_want, int iMax, int jMax, double S_max)
{
    // Declare and Initialise local variables: dS,dt, Smax:
    double dS = S_max/jMax;
    std::vector<double> time_grid_points = bunch_time_points(t_0, iMax);
    double iMax_star = time_grid_points.size()-1;
    // std::cout << "dS: " << dS << " dt: " << dt << std::endl;
    
    // First appply terminal condition V(S,t = T):
    std::vector<double> vNew(jMax+1), S(jMax+1);
    for (int j{}; j<=jMax; j++) S[j] = j*dS;
    for (int j{}; j<=jMax; j++) vNew[j] = max(F, R*S[j]);
    std::vector<double> vOld = vNew;
    

    for (int k = iMax_star-1; k >= 0; k--)
    {
        double t = time_grid_points[k];
        double dt = time_grid_points[k+1] - time_grid_points[k];
        // std::cout << "t: " << t << " dt: " << dt << std::endl;
        
        std::vector<double> a(jMax+1), b(jMax+1), c(jMax+1), d(jMax+1);
        
        // Set up boundary condition at S = 0:
        a[0] = 0;
        b[0] = -(1./dt+K*theta(t)/dS+r/2);
        c[0] = K*theta(t)/dS;
        d[0] = (-1./dt+r/2)*vOld[0]-C*exp(-alpha*t);
        
        for (int j{1}; j<jMax; j++)
        {
            a[j] = 0.25*(sigma*sigma*pow(j, 2*beta)*pow(dS, 2*(beta-1))-K*(theta(t)/dS-j));
            b[j] = -0.5*sigma*sigma*pow(j, 2*beta)*pow(dS, 2*(beta-1)) - 0.5*r - 1./dt;
            c[j] = 0.25*(sigma*sigma*pow(j, 2*beta)*pow(dS, 2*(beta-1))+K*(theta(t)/dS-j));
            d[j] = -a[j]*vOld[j-1] - (-0.5*sigma*sigma*pow(j, 2*beta)*pow(dS, 2*(beta-1)) - 0.5*r +                         1./dt)*vOld[j] - c[j]*vOld[j+1] - C*exp(-alpha*t);
        }
        
        // Set up boundary condition at S -> inf.
        a[jMax] = 0;
        b[jMax] = 1;
        c[jMax] = 0;
        d[jMax] = R*jMax*dS;
        
        // Apply penalty method:
        double penalty=1.e8, tolerance = 1.e-8; int iterMax = 100, penaltyIt = 0;
        for(int penaltyIt=0; penaltyIt< iterMax; penaltyIt++)
        {
            std::vector<double> bHat(b),dHat(d);
            for(int j=1;j<jMax;j++)
            {
                // turn on penalty if V < RS for conversion:
                if(vNew[j] < S[j]*R)
                {
                    bHat[j] = b[j] - penalty;
                    dHat[j] = d[j] - penalty*S[j]*R;
                }
                
                // If t < t_0, turn on penalty if V > max(Cp, RS) for call Option:
                if (t < t_0 && vNew[j] > max(Cp, R*S[j]) )
                {
                  bHat[j] = b[j] - penalty;
                  dHat[j] = d[j] - penalty*max(Cp, R*S[j]);
                }
                
            }
            // solve matrix equations with LU decompomposition:
            std::vector<double> y = LU_decompose(a,bHat,c,dHat);
            
            // calculate difference from last time
            double error=0.;
            for(int j=0;j<=jMax;j++) error += fabs(vNew[j] - y[j]);
            vNew = y;
            if(error < tolerance) break;
        }
        if (penaltyIt >= iterMax)
            {std::cout << "Failed to converge after " << iterMax << " iterations. " << '\n';}
        vOld = vNew;
    }
    double value = lagrangeInterpolation(vNew,S,S_want,4);
    // std::cout << "V(S=" << S_want << ", t=0) = " << value << std::endl;
    return value;
}

void results_table()
{
    std::cout << "iMax jMax        V          Diff         Time  " << std::endl;
    
    double vOld{1};
    for (int n{1}; n<=30; n+=3)
    {
        int iMax = n*76, jMax = n*76;
        auto start = std::chrono::high_resolution_clock::now();
        double vNew = Crank_Nicolson_American_Bond_bunched(95.24, iMax, jMax, 4*F/R);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        
        double diff = vNew - vOld;
        
        std::cout << iMax << "   " << jMax << "    " << vNew << "   " << fabs(diff) << "    " << duration.count()/1000. << std::endl;
        vOld = vNew;
    }
}

int main()
{
    
    std::ofstream DataFile;
    DataFile.open("t_bunched.csv");
    
    std::vector<double> bunched_t = bunch_time_points(t_0, 60);
    for (int i{}; i<bunched_t.size(); i++)
    {
        DataFile.width(10);  DataFile <<  bunched_t[i] << std::endl;
    }
    
    std::cout << std::setprecision(10);
   
    // Crank_Nicolson_American_Bond_bunched(95.24, 1292, 1292, 4*F/R);
    
    results_table();
    
    // DataFile.width(20);  DataFile << S_want;
    // DataFile.width(20); DataFile << value << std::endl;
    
    return 0;
}

