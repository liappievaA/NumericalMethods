#include <cmath>
#include <iostream>
#include "func.hpp"

double u(double x){
    // return 5*cos(M_PI/2. * (4. + 1.) * x) - 10 * cos(M_PI/2. * (10. + 1.) * x);
    return x*x*x*x - x*x*x;
    // return cos(M_PI/2. * (2. + 1.) * x);
    // return 1. / (double)(1. + 25.*x*x) - 1./26.;
    // return 1 ? x < 0.5 && x > 0.4: 0;
}
int makePoints(int NUMBER_OF_DOTS, double *MASSIVE_OF_DOTS){
    double temp = 0, h = 1. / ((double)NUMBER_OF_DOTS - 0.5);
    int i;
    if (NUMBER_OF_DOTS <= 2){
        std::cout << "NOT ENOUGH NUMBER OF DOTS" << std::endl;
        return -1;
    }
    MASSIVE_OF_DOTS[0] = 0.;
    MASSIVE_OF_DOTS[NUMBER_OF_DOTS] = 1. + h / 2.;    
    temp = h;
    for (i = 1; i < NUMBER_OF_DOTS; i++){
        MASSIVE_OF_DOTS[i] = temp;
        temp += h;
    }
    return 1;
}
double basis_function(int m, int k, int N){
    double h = 1 / ((double)N - 0.5);
    return cos(M_PI / 2. * (2. * m + 1.) * (k * h));

}
void write_bfm_for_all_points(int m, int N, double *All_k){
    for(int k = 0; k < N; k++){
        All_k[k] = basis_function(m, k, N);
    }
}
int makeKoef_C_nk(int NUMBER_OF_DOTS, double *MASSIVE_OF_U, double *MASSIVE_C_nk, double *_trashMassive){
    double result = 0, h =  1. / ((double)NUMBER_OF_DOTS - 0.5);
    int m, k;
    for(k = 0; k < NUMBER_OF_DOTS; k++){
        write_bfm_for_all_points(k, NUMBER_OF_DOTS, _trashMassive);
        result = 0;
        for(m = 1; m < NUMBER_OF_DOTS; ++m){
            result += _trashMassive[m] * MASSIVE_OF_U[m] * h;
        } 
        result += h / 2. * _trashMassive[0] * MASSIVE_OF_U[0];
        MASSIVE_C_nk[k] = 2. * result;
    }
    return 1;
}
double seriesOfFurierAtPoint(int N, double *koef, double var){
    double result = 0;
    int i;
    for(i = 0; i < N; i++){
        result += koef[i] * cos(M_PI / 2. * (2. * i + 1.) * var);
    }
    return result;
}
// double my_m_eps(){
//     double eps = 1.0;
//     while ((1.0 + eps) < 1.0) {
//         eps /= 2.0;
//     }
//     return eps;
// }
