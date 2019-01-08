#include <vector>     // std::vector
#include <math.h>     // pow, sqrt, fabs
#include <algorithm>  // std::nth_element, std::min_element, std::max_element
#include <cmath>      // std::copysign
#include <iostream>   // std::cout, std::cerr, std::endl

double sign(double x){
    if(x == 0){
        return 1.0;
    }
    else{
        return std::copysign(1.0, x);
    }
}

// Accurate Monotonicity Preserving Cubic Interpolation
// In this case x, y is a cdf, and always monotone increasing
// https://epubs-siam-org.proxy.lib.utk.edu/doi/pdf/10.1137/0904045
std::vector<double> spline(std::vector<double> x, std::vector<double> y, std::vector<double> t){
    // x - monotone increasing
    // evaluate at t
    // note: using i=0, ..., n-1 instead of 1, ..., n as in paper

    int n = x.size();
    std::vector<double> dx(n-1);
    std::vector<double> S(n);
    double dy;
    for(int i=0; i < n-1; i++){
        dx[i] = x[i+1] - x[i];
        dy = y[i+1] - y[i];
        S[i] = dy / dx[i];
    }

    S[n-1] = S[n-2];


    // dy/dx or f_dot approximation
    std::vector<double> f_dot(n);
    double f_dot_i;
    double S_i_min, S_i_max;
    for(int i=0; i < n; i++){
        //j = i+1/2

        // boundaries
        if(i==0){
            f_dot_i = ( (2*dx[i] + dx[i+1])*S[i] - dx[i]*S[i+1] )/( dx[i] + dx[i+1] );
        }
        else if(i==n-1){
            f_dot_i = ( (2*dx[i-1] + dx[i-2])*S[i-1] - dx[i-1]*S[i-2] )/( dx[i-1] + dx[i-2] );
        }
        else{
            // Use Fritsch-Butland
            S_i_min = std::min(S[i], S[i-1]);
            S_i_max = std::max(S[i], S[i-1]); 
        
            f_dot_i = 3*S_i_min * S_i_max / (S_i_max + 2*S_i_max); // ( dx[i-2]*S[i-1] + dx[i-1]*S[i-2] ) / ( x[i] - x[i-2]  );
        }
        
        // Hyman filter: eqn 2.6
        double sigma = sign(f_dot_i);
        if(sigma > 0.0){
            f_dot_i = std::min( std::max(0.0, f_dot_i),  3*std::min(fabs(S[i-1]), fabs(S[i])) );
        } 
        else{
            f_dot_i = std::max( std::min(0.0, f_dot_i), -3*std::min(fabs(S[i-1]), fabs(S[i])) );
        }
        f_dot[i] = f_dot_i;
    }
    
    // evalute
    int t_n = t.size();
    std::vector<double> P(t_n);
    int i, j;
    double t_j; 
    double c1, c2, c3, c4, p;

    for(j=0; j<t_n; j++){
        // x[i] <= t[j] < x[i+1]
        t_j = t[j];
        i = std::lower_bound(x.begin(), x.end(), t_j) - x.begin() -1;

        c1 = y[i];
        c2 = f_dot[i];
        c3 = (3*S[i] - f_dot[i+1] - 2*f_dot[i]) / (dx[i]);
        c4 = (2*S[i] - f_dot[i+1] - f_dot[i]) / pow(dx[i], 2.0); 
        p = c1 + (t_j - x[i])*c2 + pow(t_j - x[i], 2.0)*c3 + pow(t_j - x[i], 3.0)*c4;

        P[j] = p;
    }
    return P;
}

int main(){
    std::vector<double> x = {7.99, 8.09, 8.19, 8.7, 9.2, 10, 12.0, 15.0, 20.0};
    std::vector<double> y = {0, 2.76429E-05, 0.0437498, 0.169183, 0.469428, 0.94374, 0.998636, 0.999919, 0.999994};
    std::vector<double> t = {8.0, 8.20, 9.0, 10.2, 12.2, 16.0, 17.0, 18.0};
    
    std::vector<double> P = spline(x, y, t);

    int n = P.size();
    for(int i=0; i < n; i++){
        std::cout << t[i] << "\t" << P[i] << std::endl;
    }

    return 0;
}

