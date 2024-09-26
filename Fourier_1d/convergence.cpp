#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "func.hpp"

const int NUMBER_OF_DOTS = 50; 
const int TEST_FOR_CONVERGE = 10; 

int main(){
    int i, j;
    double hForNorm, MaxForNorm, temp;
    std::string nameFile, GNUPLOT;
    std::ofstream outFile;
    double *MASSIVE_C_nk;
    double *_trashMassive;
    double *MASSIVE_OF_U;
    double *MASSIVE_OF_DOTS;
    double *MASSIVE_NUMBER_OF_DOTS;
    double *MASSIVE_LOG_HEIGHT_ERROR;
    try {
        MASSIVE_NUMBER_OF_DOTS = new double[TEST_FOR_CONVERGE];
        MASSIVE_LOG_HEIGHT_ERROR = new double[TEST_FOR_CONVERGE];
        MASSIVE_C_nk = new double[TEST_FOR_CONVERGE * (NUMBER_OF_DOTS + 1)];
        _trashMassive = new double[TEST_FOR_CONVERGE * (NUMBER_OF_DOTS + 1)];
        MASSIVE_OF_U = new double[TEST_FOR_CONVERGE *(NUMBER_OF_DOTS + 1)];
        MASSIVE_OF_DOTS = new double[TEST_FOR_CONVERGE * (NUMBER_OF_DOTS + 1)];
	} catch(...){
		std::cout << "Some trouble with memory" << std::endl;
		return -2;
	}
    // заполняем массивы
    for(i = 1; i < TEST_FOR_CONVERGE + 1; i++){
        if(makePoints((NUMBER_OF_DOTS + 1) * i - 1, MASSIVE_OF_DOTS) != 1){
            std::cout << "ERROR IN func.cpp" << std::endl;
            delete [] MASSIVE_NUMBER_OF_DOTS;
            delete [] MASSIVE_LOG_HEIGHT_ERROR;
            delete [] MASSIVE_C_nk;
            delete [] MASSIVE_OF_U;
            delete [] _trashMassive;
            delete [] MASSIVE_OF_DOTS;
            return -3;
        }
        for(j = 0; j < (NUMBER_OF_DOTS + 1) * i; j++){
            MASSIVE_OF_U[j] = u(MASSIVE_OF_DOTS[j]);
        }
        makeKoef_C_nk((NUMBER_OF_DOTS + 1) * i - 1, MASSIVE_OF_U, MASSIVE_C_nk, _trashMassive);
        hForNorm = 1 / ((double)(NUMBER_OF_DOTS * i));
        MaxForNorm = 1.e-16;
        for (double x = 0; x < 1; x += hForNorm){
            temp = seriesOfFurierAtPoint((NUMBER_OF_DOTS + 1) * i - 1, MASSIVE_C_nk, x) - u(x);
            if (fabs(temp) > MaxForNorm)
                MaxForNorm = fabs(temp);
        }
        MASSIVE_LOG_HEIGHT_ERROR[i - 1] = log(MaxForNorm);
        MASSIVE_NUMBER_OF_DOTS[i - 1] = log(NUMBER_OF_DOTS * i);
        // MASSIVE_LOG_HEIGHT_ERROR[i - 1] = MaxForNorm;
        // MASSIVE_NUMBER_OF_DOTS[i - 1] = NUMBER_OF_DOTS * i;
    }
    nameFile = "out.txt";
    outFile.open(nameFile);
    if(outFile.is_open()){
        for(i = 0; i < TEST_FOR_CONVERGE; i++){
            outFile << std::setprecision(15) << MASSIVE_NUMBER_OF_DOTS[i] << " " << MASSIVE_LOG_HEIGHT_ERROR[i] << std::endl;
        }
    } else {
        std::cout << "ERROR: u cant create out.txt" << std::endl;
        delete [] MASSIVE_LOG_HEIGHT_ERROR;
        delete [] MASSIVE_NUMBER_OF_DOTS;
        delete [] MASSIVE_OF_U;
        delete [] _trashMassive;
        delete [] MASSIVE_C_nk;
        delete [] MASSIVE_OF_DOTS;
        return -4;
    }
    outFile.close();
    //рисуем
    nameFile = "gnu.txt";
    outFile.open(nameFile);
    if(outFile.is_open()){
        outFile << "set terminal png size 1000,1000 \n" << std::endl;
        outFile << "set output \"test.png\" \n" << std::endl;
        // outFile << "set xlabel 'number of dots' \n" << std::endl;
        // outFile << "set ylabel 'error' \n" << std::endl;
        outFile << "set xlabel 'log(number of dots)' \n" << std::endl;
        outFile << "set ylabel 'log(error)' \n" << std::endl;
        outFile << "plot 'out.txt' u 1:2 w linesp title 'Dependence'\\" << std::endl;
    } else {
        std::cout << "ERROR: u cant create gnu.txt" << std::endl;
        delete [] MASSIVE_LOG_HEIGHT_ERROR;
        delete [] MASSIVE_NUMBER_OF_DOTS;
        delete [] MASSIVE_OF_U;
        delete [] _trashMassive;
        delete [] MASSIVE_C_nk;
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
        delete [] MASSIVE_LOG_HEIGHT_ERROR;
        delete [] MASSIVE_NUMBER_OF_DOTS;
        delete [] MASSIVE_OF_U;
        delete [] _trashMassive;
        delete [] MASSIVE_C_nk;
        delete [] MASSIVE_OF_DOTS;
        return -6;
    }
    for(i = 1; i < TEST_FOR_CONVERGE; i++){
        if(fabs(MASSIVE_NUMBER_OF_DOTS[i] - MASSIVE_NUMBER_OF_DOTS[i - 1]) < 1.e-16){
            std::cout << "ERROR: u have zero value at log" << std::endl;
            delete [] MASSIVE_LOG_HEIGHT_ERROR;
            delete [] MASSIVE_NUMBER_OF_DOTS;
            delete [] MASSIVE_OF_U;
            delete [] _trashMassive;
            delete [] MASSIVE_C_nk;
            delete [] MASSIVE_OF_DOTS;
            return -7;
        }
    }
    temp = 0;
    for(i = 1; i < TEST_FOR_CONVERGE; i += 2){
        temp += (2. / TEST_FOR_CONVERGE) * ((MASSIVE_LOG_HEIGHT_ERROR[i - 1] - MASSIVE_LOG_HEIGHT_ERROR[i]) / (MASSIVE_NUMBER_OF_DOTS[i] - MASSIVE_NUMBER_OF_DOTS[i - 1]));
    }
    std::cout << temp << std::endl;
    delete [] MASSIVE_LOG_HEIGHT_ERROR;
    delete [] MASSIVE_NUMBER_OF_DOTS;
    delete [] MASSIVE_OF_U;
    delete [] _trashMassive;
    delete [] MASSIVE_C_nk;
    delete [] MASSIVE_OF_DOTS;
    return 1;
}

