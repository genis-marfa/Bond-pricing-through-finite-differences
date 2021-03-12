#include <vector>
#include <cassert>
#include <fstream>

template <class T>
T max(T n1, T n2)
{
    return (n1 > n2) ? n1 : n2;
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
    if(x.size()<n)return lagrangeInterpolation(y,x,x0,x.size());
    if(n==0)throw;
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
// Define parameters globally:
double T = 3., F = 190., R = 2., r = 0.0043, K = 0.08333333, mu = 0.0076, X = 95.24, C = 0.408, alpha = 0.03, beta = 0.142, sigma = 17.1;

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
double Crank_Nicolson_Bond(double S_want, int iMax, int jMax, double S_max)
{
    // Declare and Initialise local variables: dS,dt, Smax:
    // double S_max = 5*F/R;
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
        d[jMax] = (jMax*dS)*A(i*dt)+B(i*dt);
        
        // Call solver:
        vNew = LU_decompose(a, b, c, d);
        vOld = vNew;
    }
    return lagrangeInterpolation(vNew,S,S_want,4);;
}

double N(double x){return 0.5*erfc(-x/sqrt(2));}

double CallOption(double St, double X, double T, double r, double D, double sigma, double t)
{
    double d1 = (log(St/X)+(r-D+0.5*sigma*sigma)*(T-t) )/ (sigma*sqrt(T-t));
    double d2 = d1-sigma*sqrt(T-t);
    return St*exp(-D*(T-t)) * N(d1) - X * exp(-r*(T-t))*N(d2);
}

void analytic_estimate()
{
    double est{};
    std::ofstream DataFile;
    DataFile.open("Values_an.csv");
    DataFile.width(20);  DataFile << "S";
    DataFile.width(20); DataFile << "Value" << std::endl;
    
    for (int S = 0; S<200; S+=20)
    {
        double call =  CallOption(S, F/R, T, r, r, sigma, 0);
        est = R*call+F*exp(-r*T)+C/(alpha+r)*(1-exp(-(alpha+r)*T));
        DataFile.width(20);  DataFile << S;
        DataFile.width(20);  DataFile << est << std::endl;
    }
    std::cout << "Succesfully saved to file." << std::endl;
}

void get_abs_err()
{
    double call =  CallOption(25, F/R, T, r, r, sigma, 0);
    double est = R*call+F*exp(-r*T)+C/(alpha+r)*(1-exp(-(alpha+r)*T));
    
    std::ofstream DataFile;
    DataFile.open("Abs_Err_table.csv");
    DataFile.width(20);  DataFile << "iMax";
    DataFile.width(20);  DataFile << "jMax";
    DataFile.width(20);  DataFile << "dt";
    DataFile.width(20);  DataFile << "dS";
    DataFile.width(20);  DataFile << "Value_50";
    DataFile.width(20); DataFile << "Abs_Err" << std::endl;
    
    for (int n{1}; n <= 30; n++)
    {
        double Smax = 5*F/R;
        int jmax = Smax/25*n;
        double vNumerical = Crank_Nicolson_Bond(25, 200, jmax, Smax);
        double dS =  Smax/jmax;
        int j_want = 25/dS;
        
        if (j_want*dS == 25)
        {
            DataFile.width(20);  DataFile << jmax;
            DataFile.width(20);  DataFile << jmax;
            DataFile.width(20);  DataFile << 200;
            DataFile.width(20);  DataFile << dS;
            
            double num = vNumerical;
            DataFile.width(20);  DataFile << num;
            double abs_err = fabs(est-num);
            DataFile.width(20);  DataFile << abs_err << std::endl;
        }
    }
    std::cout << "Succefully saved to file." << std::endl;
}

void get_convergange_rate()
{
    int k = 2;
    
    double vOld{1}, diffOld{1};
    std::cout << "iMax=jMax" << "    " << "vNew" << "     " << "     Extrap   " << "   diffNew   "  << "    " << "    R" << "        " << "      c" << "   " << "    Duration(s)     " << std::endl;
    for (int jmax{5}; jmax<3000; jmax*= k)
    {
        auto start = std::chrono::high_resolution_clock::now();
        double vNew = Crank_Nicolson_Bond(95.24, jmax, jmax, 4*F/R);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        
        double diffNew = vNew - vOld;
        double R = diffOld/diffNew;
        double c = log(R)/log(k);
        double extrapValue = (4*vNew - vOld)/3.;
        std::cout << jmax << "        " << vNew << "   " << extrapValue << "  " << diffNew  << "   " << R << "   " << c << "   " << duration.count()/1000. << std::endl;
        
        vOld = vNew;
        diffOld = diffNew;
    }
}

void get_convergange_rate_extrap()
{
    int k = 2;
    
    double vOld{1}, diffOld{1}, extrap_old{1};
    
    std::cout << "iMax=jMax    vNew     Diff    Extrap    Diff_Extrap    Rat   " << std::endl;
    
    for (int jmax{20}; jmax<3000; jmax*= k)
    {
        auto start = std::chrono::high_resolution_clock::now();
        double vNew = Crank_Nicolson_Bond(45 , jmax, jmax, 5*F/R);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        
        double diffNew = vNew - vOld;
        double R = diffOld/diffNew;
        double extrapValue = (4*vNew - vOld)/3.;
        double diffExtrap = extrapValue-extrap_old;
        
        std::cout << jmax << "        " << vNew << "    " << diffNew << "   " <<  extrapValue << "  " << diffExtrap  << "     " << duration.count()/1000 << std::endl;
        
        vOld = vNew;
        diffOld = diffNew;
        extrap_old = extrapValue;
    }
}


void testing_Smax()
{
    double diff{1}; double vOld{1};
    std::cout << "n    Smax   iMax   jMax   V(45)    Diff    Duration(s)" << std::endl;
    for (int n{2}; n<=10; n++)
    {
        double S_max = n*F/R;
        int jMax = S_max/0.5;
        
        auto start = std::chrono::high_resolution_clock::now();
        double vNew = Crank_Nicolson_Bond(45, 100, jMax, S_max);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        
        diff = vNew - vOld;
        
        std::cout << n << "    " << S_max << "    " << 100 << "    " << jMax << "    " << vNew << "  " << diff << "   " << duration.count()/1000. << std::endl;
        
        vOld = vNew;
    }
}
