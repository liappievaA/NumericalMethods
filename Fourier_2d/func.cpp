#include <cmath>
#include <iostream>
#include "func.hpp"

double u(double x, double y){
    // return (x*x*x - x*x)*(y*y*y - y*y);
      return cos(M_PI/2. * (2. + 1.) * x) * cos(M_PI/2. * (2. + 1.) * y);
    // return 1 ? x < 0.5 && x > 0.4 && y > 0.4 && y < 0.5: 0;
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
    // for (i = 0; i < NUMBER_OF_DOTS + 1; i++){
    //     std::cout<<MASSIVE_OF_DOTS[i]<<"\n";
    // }
    return 1;
}
int makeKoef_C_nk(int NUMBER_OF_DOTS, double *MASSIVE_OF_U, double *MASSIVE_OF_DOTS, double *MASSIVE_C_nk, double *_trashMassive){
    double result = 0, h =  1. / ((double)NUMBER_OF_DOTS - 0.5);
    int m, k;
    for(m = 0; m < NUMBER_OF_DOTS + 1; m++){
        result = 0;
        for(k = 0; k < NUMBER_OF_DOTS; k++){
            _trashMassive[k] = cos(M_PI / 2. * (2. * m + 1.) * MASSIVE_OF_DOTS[k]);
        }
        for(k = 1; k < NUMBER_OF_DOTS; k++){
            result += _trashMassive[k] * MASSIVE_OF_U[k] * h;
        } 
        result += h / 2. * _trashMassive[0] * MASSIVE_OF_U[0];
        MASSIVE_C_nk[m] = 2. * result;
    }
    return 1;
}
int make_C_mn(int NUMBER_OF_DOTS, double *MASSIVE_OF_DOTS, 
                        double *_trashMatrix, double *_trashMassive_1, double *_trashMassive_2, 
                        double *MATRIX_OF_U, double *MATRIX_C_nk){
    int i, j;
    for(i = 0; i < NUMBER_OF_DOTS + 1; i++){
        makeKoef_C_nk(NUMBER_OF_DOTS, MATRIX_OF_U + i*(NUMBER_OF_DOTS+1), 
                        MASSIVE_OF_DOTS, _trashMatrix + i*(NUMBER_OF_DOTS+1), _trashMassive_1);
    }

    for(i = 0; i < NUMBER_OF_DOTS + 1; i++){
        for(j = 0; j < NUMBER_OF_DOTS + 1; j++){
            _trashMassive_1[j] = _trashMatrix[j*(NUMBER_OF_DOTS + 1) + i];
        }

        makeKoef_C_nk(NUMBER_OF_DOTS, _trashMassive_1, MASSIVE_OF_DOTS, 
                        MATRIX_C_nk + i*(NUMBER_OF_DOTS + 1), _trashMassive_2);
    }

    return 1;
}
double seriesOfFurierAtPoint(int N, double *koef, double var_x, double var_y){
    double result = 0;
    int i, j;
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++){
            result += koef[i*(N+1) + j] * cos(M_PI / 2. * (2. * i + 1.) * var_x) * cos(M_PI / 2. * (2. * j + 1.) * var_y);
        }
    }
    return result;
}
