#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include "func.hpp"

const int NUMBER_OF_DOTS = 50; 

int main(){
    int i;
    std::string nameFile, GNUPLOT;
    std::ofstream outFile;
    double *MASSIVE_C_nk;
    double *_trashMasssive;
    double *MASSIVE_OF_U;
    double *MASSIVE_OF_DOTS;
    try {
        MASSIVE_C_nk = new double[NUMBER_OF_DOTS + 1];
        _trashMasssive = new double[NUMBER_OF_DOTS + 1];
        MASSIVE_OF_U = new double[NUMBER_OF_DOTS + 1];
        MASSIVE_OF_DOTS = new double[NUMBER_OF_DOTS + 1];
	} catch(...){
		std::cout << "Some trouble with memory" << std::endl;
		return -2;
	}
    if(makePoints(NUMBER_OF_DOTS, MASSIVE_OF_DOTS) != 1){
        std::cout << "ERROR IN funk.cpp" << std::endl;
        delete [] MASSIVE_C_nk;
        delete [] MASSIVE_OF_U;
        delete [] _trashMasssive;
        delete [] MASSIVE_OF_DOTS;
        return -3;
    }
    for(i = 0; i < NUMBER_OF_DOTS + 1; i++){
        MASSIVE_OF_U[i] = u(MASSIVE_OF_DOTS[i]);
    }
    makeKoef_C_nk(NUMBER_OF_DOTS, MASSIVE_OF_U, MASSIVE_C_nk, _trashMasssive);
    nameFile = "out.txt";
    outFile.open(nameFile);
    //добавим в файл промежуточные точки для оценки погрешности
    double dot1;
//     double dot2;
    if(outFile.is_open()){
        for(i = 0; i < NUMBER_OF_DOTS; i++){ 
            dot1 = MASSIVE_OF_DOTS[i];
//             dot2 = MASSIVE_OF_DOTS[i + 1];
            outFile << std::setprecision(15) << dot1 << " " << seriesOfFurierAtPoint(NUMBER_OF_DOTS, MASSIVE_C_nk, dot1) << " " << u(dot1) << std::endl;
//             outFile << std::setprecision(15) << 2*dot1 / 3. + dot2 / 3. << " " << seriesOfFurierAtPoint(NUMBER_OF_DOTS, MASSIVE_C_nk, 2*dot1 / 3. + dot2 / 3.) << " "  << u(2*dot1 / 3. + dot2 / 3.) << std::endl;
//             outFile << std::setprecision(15) << dot1 / 3. + 2*dot2 / 3. << " " << seriesOfFurierAtPoint(NUMBER_OF_DOTS, MASSIVE_C_nk, dot1 / 3. + 2*dot2 / 3.) << " " << u(dot1 / 3. + 2*dot2 / 3.) <<  std::endl;
        }
    } else {
        std::cout << "ERROR: u can't create out.txt" << std::endl;
        delete [] MASSIVE_OF_U;
        delete [] _trashMasssive;
        delete [] MASSIVE_C_nk;
        delete [] MASSIVE_OF_DOTS;
        return -4;
    }
    dot1 = MASSIVE_OF_DOTS[NUMBER_OF_DOTS];
    outFile << std::setprecision(15) << dot1 << " " << seriesOfFurierAtPoint(NUMBER_OF_DOTS, MASSIVE_C_nk, dot1) << " " << u(dot1) << std::endl;
    outFile.close();
    //рисуем
    nameFile = "gnu.txt";
    outFile.open(nameFile);
    if(outFile.is_open()){
        outFile << "set terminal png size 1000,1000 \n" << std::endl;
        outFile << "set output \"test.png\" \n" << std::endl;
        outFile << "plot 'out.txt' u 1:2 w linesp title 'Fourier', 'out.txt' u 1:3 w linesp title 'u' pt -1 \\" << std::endl;
    } else {
        std::cout << "ERROR: u cant create gnu.txt" << std::endl;
        delete [] MASSIVE_OF_U;
        delete [] _trashMasssive;
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
        delete [] MASSIVE_OF_U;
        delete [] _trashMasssive;
        delete [] MASSIVE_C_nk;
        delete [] MASSIVE_OF_DOTS;
        return -6;
    }
    delete [] MASSIVE_OF_U;
    delete [] _trashMasssive;
    delete [] MASSIVE_C_nk;
    delete [] MASSIVE_OF_DOTS;
    return 1;
}
