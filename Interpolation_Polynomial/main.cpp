#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include "cmath"
#include "func.hpp"

const int TYPE_OF_DOTS = 3;
const int NUMBER_OF_DOTS = 10;
const double LEFT_SIDE = -1.0;
const double RIGHT_SIDE = 1.0;
const double M_EPS = 1.e-12;

int main(){
    int i, j, k;
    std::string nameTypeDots[3]{"Uniform", "Chebyshev", "Random"};
    std::string nameFile, GNUPLOT;
    std::ofstream outFile;
    double *x, *y;
    double *coeffOfPolinom, *matrix;
    double **MASSIVE_OF_DOTS;
    double *_trashMassive;
    if(NUMBER_OF_DOTS < 2 || LEFT_SIDE - RIGHT_SIDE > M_EPS){
        std::cout << "N, LEFT_SIDE or RIGHT_SIDE are incorrect!" << std::endl;
        return -1;
    }
    //---------Выделяем память----------
    try {
         x = new double[NUMBER_OF_DOTS];
         y = new double[NUMBER_OF_DOTS];
         matrix = new double[NUMBER_OF_DOTS * NUMBER_OF_DOTS];
         coeffOfPolinom = new double[NUMBER_OF_DOTS];
        _trashMassive = new double[NUMBER_OF_DOTS];
        MASSIVE_OF_DOTS = new double*[TYPE_OF_DOTS];
        for(i = 0; i < 3; i++){
            MASSIVE_OF_DOTS[i] = new double[(3)*(NUMBER_OF_DOTS - 1) + 1];
        }
	} catch(...){
		std::cout << "Some trouble with memory" << std::endl;
		return -2;
	}
    //----------------------------------------Заполняем массивы----------------------------------------------
    if(makePoints(NUMBER_OF_DOTS, RIGHT_SIDE, LEFT_SIDE, MASSIVE_OF_DOTS, TYPE_OF_DOTS, _trashMassive) != 1){
        std::cout << "ERROR IN makePoints.cpp" << std::endl;
         delete [] coeffOfPolinom;
         delete [] y;
         delete [] x;
         delete [] matrix;
        delete [] _trashMassive;
        for (int i = 0; i < TYPE_OF_DOTS; i++) {
            delete [] MASSIVE_OF_DOTS[i];
        }
        delete [] MASSIVE_OF_DOTS;
        return -3;
    }
    //--Записываем все данные в файл--
    for(i = 0; i < TYPE_OF_DOTS; i++){
        for(j = 0; j < 3*(NUMBER_OF_DOTS - 1) + 1; j+=3){
            x[j / 3] = MASSIVE_OF_DOTS[i][j];
            y[j / 3] = f(x[j / 3]);
        }
        for(j = 0; j < NUMBER_OF_DOTS; j++){
            for(k = 0; k < NUMBER_OF_DOTS; k++){
                matrix[j*NUMBER_OF_DOTS + k] = __pow(x[j], k);
            }
        }
        findSolutionWithJordanGauss(NUMBER_OF_DOTS, matrix, y, coeffOfPolinom);
        for(j = 0; j < 3*(NUMBER_OF_DOTS - 1) + 1; j+=3){
            x[j / 3] = MASSIVE_OF_DOTS[i][j];
            y[j / 3] = f(x[j / 3]);
        }
        nameFile = "out" + std::to_string(i) + ".txt";
        outFile.open(nameFile);
        if(outFile.is_open()){
            for(j = 0; j < 3*(NUMBER_OF_DOTS - 1) + 1; j++){
                outFile << std::setprecision(15)
                        << MASSIVE_OF_DOTS[i][j] << " "
                        << f(MASSIVE_OF_DOTS[i][j]) << " "
                        << lagrange(NUMBER_OF_DOTS, MASSIVE_OF_DOTS[i][j], y, x) << " "
                        << canonPolinom(NUMBER_OF_DOTS, MASSIVE_OF_DOTS[i][j], coeffOfPolinom) << " " 
                        << std::endl;
            }
        }
        outFile.close();
    }
    //------Рисуем-------
    nameFile = "gnu.txt";
    outFile.open(nameFile);
    if(outFile.is_open()){
        outFile << "set multiplot \\" << std::endl //несколько графиков на одном рисунке
                << "title 'IMAGE' \\" << std::endl //заголовок всего изображения
                << "layout 1, 3 \\" << std::endl << std::endl; //три графика в ряд
        for(i = 0; i < TYPE_OF_DOTS; i++){ //каждая итерация строит один график
            outFile << "set title '"<< nameTypeDots[i].c_str() << "'" <<std::endl; //заголовок каждого графика
            outFile << "plot 'out" << i << ".txt' u 1:3 w linesp title 'L' lt rgb 'red' pt 2 ps 2, 'out" //большие крестики 
                                    << i << ".txt' u 1:4 w linesp title 'P' pt 1 ps 2 lt rgb 'green', 'out" 
                                    << i << ".txt' u 1:2 w linesp title 'f(x)' lw 2 lt rgb 'black' pt -1" << std::endl << std::endl; // маленькие точки
        }
        outFile << "unset multiplot" << std::endl;
    }
    outFile.close();
    GNUPLOT = "gnuplot -persist \"" + nameFile + "\"";
    FILE *pipe = popen(GNUPLOT.c_str() , "w");
    if (pipe != NULL){
        fflush(pipe);
        pclose(pipe);
    }
    delete [] coeffOfPolinom;
    delete [] y;
    delete [] x;
    delete [] matrix;
    delete [] _trashMassive;
    for (int i = 0; i < TYPE_OF_DOTS; i++) {
        delete [] MASSIVE_OF_DOTS[i];
    }
    delete [] MASSIVE_OF_DOTS;
    return 1;
}
