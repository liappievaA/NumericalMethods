#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include "stdio.h"
#include <string>
#include "math.h"
#include "func.hpp"

const int NUMBER_OF_DOTS = 20;
const int TEST_FOR_CONVERGE = 2;
const int MACHINE_EPS = 1.e-16;

int main(){
    int i, j, l, temp;
    double hForNorm, MaxForNorm, trash;
    std::string nameFile, GNUPLOT;
    std::ofstream outFile;
    double *MATRIX_C_nk;
    double *MATRIX_OF_U;
    double *_trashMatrix;

    double *_trashMassive_1;
    double *_trashMassive_2;
    double *MASSIVE_OF_U;
    double *MASSIVE_OF_DOTS;
    // -----
    double *MASSIVE_LOG_HEIGHT_ERROR;
    double *MASSIVE_NUMBER_OF_DOTS;

    try {
        MASSIVE_LOG_HEIGHT_ERROR = new double[TEST_FOR_CONVERGE];
        MASSIVE_NUMBER_OF_DOTS = new double[TEST_FOR_CONVERGE];
        MATRIX_C_nk = new double[(TEST_FOR_CONVERGE+2)*(TEST_FOR_CONVERGE+2)*(NUMBER_OF_DOTS + 1)*(NUMBER_OF_DOTS + 1)];
        MATRIX_OF_U = new double[(TEST_FOR_CONVERGE+2)*(TEST_FOR_CONVERGE+2)*(NUMBER_OF_DOTS+ 1)*(NUMBER_OF_DOTS + 1)];
        _trashMatrix = new double [(TEST_FOR_CONVERGE+2)*(TEST_FOR_CONVERGE+2)*(NUMBER_OF_DOTS + 1)*(NUMBER_OF_DOTS+ 1)];

        _trashMassive_1 = new double[(TEST_FOR_CONVERGE+2)*(NUMBER_OF_DOTS+1)];
        _trashMassive_2 = new double[(TEST_FOR_CONVERGE+2)*(NUMBER_OF_DOTS + 1)];
        MASSIVE_OF_U = new double[(TEST_FOR_CONVERGE+2)*(NUMBER_OF_DOTS+ 1)];
        MASSIVE_OF_DOTS = new double[(TEST_FOR_CONVERGE+2)*(NUMBER_OF_DOTS+1)];

        } catch(...){
		std::cout << "Some trouble with memory" << std::endl;
		return -2;
	}

    hForNorm = 1 / ((double)(NUMBER_OF_DOTS * 5));
    MaxForNorm = MACHINE_EPS;

    for(i = 1; i < TEST_FOR_CONVERGE + 1; i++){
        if(makePoints((NUMBER_OF_DOTS + 1)*i, MASSIVE_OF_DOTS) != 1){
            std::cout << "ERROR IN make_points.cpp" << std::endl;
            delete [] MASSIVE_NUMBER_OF_DOTS;
            delete [] MASSIVE_LOG_HEIGHT_ERROR;
            delete [] MATRIX_C_nk;
            delete [] MATRIX_OF_U;
            delete [] _trashMatrix;

            delete [] MASSIVE_OF_U;
            delete [] _trashMassive_1;
            delete [] _trashMassive_2;
            delete [] MASSIVE_OF_DOTS;

            return -3;
        }
        for(l = 0; l < i*(NUMBER_OF_DOTS + 1) + 1; l++){
            for(j = 0; j < i*(NUMBER_OF_DOTS + 1); j++){
                // добавить проверку на краевые условия
                temp = (l*(i*(NUMBER_OF_DOTS + 1)) + j);
                if(i == i*(NUMBER_OF_DOTS + 1)){
                    MATRIX_OF_U[temp] = -MATRIX_OF_U[(i*(NUMBER_OF_DOTS + 1) - 1)*(i*(NUMBER_OF_DOTS + 1)) + j];
                } else {
                    MATRIX_OF_U[temp] = u(MASSIVE_OF_DOTS[l], MASSIVE_OF_DOTS[j]);
                }
            }
            MATRIX_OF_U[l*(i*(NUMBER_OF_DOTS + 1)) + i*(NUMBER_OF_DOTS + 1)] = -MATRIX_OF_U[l*(i*(NUMBER_OF_DOTS + 1)) + (i*(NUMBER_OF_DOTS + 1)-1)];
        }

        make_C_mn((NUMBER_OF_DOTS+1)*i - 1, MASSIVE_OF_DOTS, _trashMatrix, _trashMassive_1, _trashMassive_2, MATRIX_OF_U, MATRIX_C_nk);
        
        for (double x = 0; x < 1; x += hForNorm){
            for(double y= 0; y < 1; y += hForNorm){
                trash = seriesOfFurierAtPoint(NUMBER_OF_DOTS, MATRIX_C_nk, x, y) - u(x, y);
                if (fabs(trash) > MaxForNorm)
                    MaxForNorm = fabs(trash);
            }

        }
        MASSIVE_LOG_HEIGHT_ERROR[i - 1] = log(MaxForNorm);
        MASSIVE_NUMBER_OF_DOTS[i - 1] = log(i*NUMBER_OF_DOTS);

    }
    nameFile = "out.txt";
    outFile.open(nameFile.c_str());
    if(outFile.is_open()){
        for(i = 0; i < TEST_FOR_CONVERGE; i++){
            outFile << std::setprecision(15)
                    << MASSIVE_NUMBER_OF_DOTS[i] << " " << MASSIVE_LOG_HEIGHT_ERROR[i] << std::endl;
        }
    } else {
        std::cout << "ERROR: u cant create out.txt" << std::endl;
        delete [] MASSIVE_NUMBER_OF_DOTS;
        delete [] MASSIVE_LOG_HEIGHT_ERROR;
        delete [] MATRIX_C_nk;
        delete [] MATRIX_OF_U;
        delete [] _trashMatrix;

        delete [] MASSIVE_OF_U;
        delete [] _trashMassive_1;
        delete [] _trashMassive_2;
        delete [] MASSIVE_OF_DOTS;
        return -4;
    }
    outFile.close();
    nameFile = "gnu.txt";
    
    outFile.open(nameFile.c_str());
    if(outFile.is_open()){
        // outFile << "set terminal png size 1000,1000 \n" << std::endl;
        // outFile << "set output \"test.png\" \n" << std::endl;
        outFile << "plot 'out.txt' u 1:2 w linesp title 'dependence' \\" << std::endl;
    } else {
        std::cout << "ERROR: u cant create gnu.txt" << std::endl;
        delete [] MASSIVE_NUMBER_OF_DOTS;
        delete [] MASSIVE_LOG_HEIGHT_ERROR;
        delete [] MATRIX_C_nk;
        delete [] MATRIX_OF_U;
        delete [] _trashMatrix;

        delete [] MASSIVE_OF_U;
        delete [] _trashMassive_1;
        delete [] _trashMassive_2;
        delete [] MASSIVE_OF_DOTS;
        return -5;
    }
    outFile.close();

    GNUPLOT = "gnuplot -persist \"" + nameFile + "\"";
    FILE *pipe = popen(GNUPLOT.c_str() , "w");

    if (pipe != NULL){
        fflush(pipe);
        pclose(pipe);
    } else {
        std::cout << "ERROR: u have fail in terminal" << std::endl;
        delete [] MASSIVE_NUMBER_OF_DOTS;
        delete [] MASSIVE_LOG_HEIGHT_ERROR;
        delete [] MATRIX_C_nk;
        delete [] MATRIX_OF_U;
        delete [] _trashMatrix;

        delete [] MASSIVE_OF_U;
        delete [] _trashMassive_1;
        delete [] _trashMassive_2;
        delete [] MASSIVE_OF_DOTS;
        return -6;
    }

    for(i = 1; i < TEST_FOR_CONVERGE; i++){
        if(abs(MASSIVE_NUMBER_OF_DOTS[i] - MASSIVE_NUMBER_OF_DOTS[i-1]) < MACHINE_EPS){
            std::cout << "ERROR: u have zero value at log" << std::endl;
            delete [] MASSIVE_NUMBER_OF_DOTS;
            delete [] MASSIVE_LOG_HEIGHT_ERROR;
            delete [] MATRIX_C_nk;
            delete [] MATRIX_OF_U;
            delete [] _trashMatrix;

            delete [] MASSIVE_OF_U;
            delete [] _trashMassive_1;
            delete [] _trashMassive_2;
            delete [] MASSIVE_OF_DOTS;
            return -7;
        } 
    }
    trash = 0;
    for(i = 1; i < TEST_FOR_CONVERGE; i+=2){
        trash += (2. / TEST_FOR_CONVERGE) * ((MASSIVE_LOG_HEIGHT_ERROR[i - 1] - MASSIVE_LOG_HEIGHT_ERROR[i]) / (MASSIVE_NUMBER_OF_DOTS[i] - MASSIVE_NUMBER_OF_DOTS[i - 1]));
    }

    std::cout << -trash << std::endl;
    delete [] MASSIVE_NUMBER_OF_DOTS;
    delete [] MASSIVE_LOG_HEIGHT_ERROR;
    delete [] MATRIX_C_nk;
    delete [] MATRIX_OF_U;
    delete [] _trashMatrix;

    delete [] MASSIVE_OF_U;
    delete [] _trashMassive_1;
    delete [] _trashMassive_2;
    delete [] MASSIVE_OF_DOTS;

    return 1;
}





