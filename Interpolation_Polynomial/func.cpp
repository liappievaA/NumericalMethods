#include "func.hpp"
#include "cmath"
#include <algorithm>
#include <iostream>

const double M_EPS = 1.e-12;

double f(double x){
    return 1/(1 + 25*x*x);
    // return sin(x);
    // return 3*x*x*x; 
    // return x*x;
    // return fabs(x);
}
//------------------------------------------Создаем точки в зависимости от выбранного типа сетки-------------------------------------------
int makePoints(int NUMBER_OF_DOTS, double RIGHT_SIDE, double LEFT_SIDE, double **MASSIVE_OF_DOTS, int TYPE_OF_DOTS, double *_trashMassive){
    int i, j;
    double temp, delta;
    int (*fill_methods[3])(int NUMBER_OF_DOTS, double RIGHT_SIDE, double LEFT_SIDE, double *_trashMassive) = {uniformDots, chebyshevDots, randomDots};
    for(i = 0; i < TYPE_OF_DOTS; i++){
        fill_methods[i](NUMBER_OF_DOTS, RIGHT_SIDE, LEFT_SIDE, _trashMassive);
        for(j = 0; j < 3*(NUMBER_OF_DOTS - 1) + 1; j+=3){
           temp = _trashMassive[j / 3];
            if(j != 0){
                delta = (temp - MASSIVE_OF_DOTS[i][j - 3]) / 3;
                MASSIVE_OF_DOTS[i][j - 2] = MASSIVE_OF_DOTS[i][j - 3] + delta;
                MASSIVE_OF_DOTS[i][j - 1] = MASSIVE_OF_DOTS[i][j - 3] + 2*delta;
            }
            MASSIVE_OF_DOTS[i][j] = temp;
        }
    }
    return 1;
}
//--------------------Первый тип сетки - равномерная сетка-------------------------------------
int uniformDots(int NUMBER_OF_DOTS, double RIGHT_SIDE, double LEFT_SIDE, double *_trashMassive){
    int i;
    double h = (double)(RIGHT_SIDE - LEFT_SIDE) / (NUMBER_OF_DOTS - 1);
    _trashMassive[0] = LEFT_SIDE;
    for(i = 1; i < NUMBER_OF_DOTS - 1; i++){
        _trashMassive[i] = LEFT_SIDE + (double)i*h;
    }
    _trashMassive[NUMBER_OF_DOTS - 1] = RIGHT_SIDE;
    return 1;
}
//-----------------------Второй тип сетки - узлы Чебышева-----------------------------------------
int chebyshevDots(int NUMBER_OF_DOTS, double RIGHT_SIDE, double LEFT_SIDE, double *_trashMassive){
    int i;
    double h;
    _trashMassive[0] = LEFT_SIDE;
    for(i = 1; i < NUMBER_OF_DOTS - 1; i++){
        h = (double)(2.0*(i + 1) - 1)/(2*NUMBER_OF_DOTS); //(i+1) потому что _trashMassive[0] уже задан и нам нужен следующий
        _trashMassive[NUMBER_OF_DOTS - i - 1] = (double)(1./2.)*(RIGHT_SIDE + LEFT_SIDE) +
                                                (double)(1./2.)*(RIGHT_SIDE - LEFT_SIDE)*cos(h*M_PI); //в начае формула замены переменной
    }
    _trashMassive[NUMBER_OF_DOTS - 1] = RIGHT_SIDE;
    return 1;
}
//-------------------------Третий тип сетки - случайные узлы-----------------------------------
int randomDots(int NUMBER_OF_DOTS, double RIGHT_SIDE, double LEFT_SIDE, double *_trashMassive){
    int i, j;
    double temp;
    _trashMassive[0] = LEFT_SIDE;
    for(i = 1; i < NUMBER_OF_DOTS - 1; i++){
        _trashMassive[i] = LEFT_SIDE + (RIGHT_SIDE - LEFT_SIDE)*((double)std::rand() / RAND_MAX); 
    }
    _trashMassive[NUMBER_OF_DOTS - 1] = RIGHT_SIDE;

     if(__checkEqual(_trashMassive, NUMBER_OF_DOTS) == -1){
             std::cout << "The same nodes were obtained, recompile the program" << std::endl;
             return -2;
         }
     for (i = NUMBER_OF_DOTS - 1; i > 0; i--) { //сортируем узлы по возрастанию алгоритом пузырька
        for (j = 0; j < i; j++) {
          if (_trashMassive[j] - _trashMassive[j + 1] > M_EPS) {
            temp = _trashMassive[j];
            _trashMassive[j] = _trashMassive[j + 1]; 
            _trashMassive[j + 1] = temp; 
          }
        }
      }
    return 1;
}
//--------------Решение СЛУ методом Жордана-Гаусса-----------------
int findSolutionWithJordanGauss(int n, double* a, double* b, double* x){
	int i, j, k, indMax;
	double tmp, max;

	for (i = 0; i < n; i++){ 
		max = fabs(a[i * n + i]);
		indMax = i;

		for (j = i + 1; j < n; j++)
			if (max < fabs(a[j * n + i])){ //максимум в столбце ниже текущей позиции
				max = fabs(a[j * n + i]);
				indMax = j;
			}

		for (j = 0; j < n; j++){ //перестановка строк
			tmp = a[i * n + j];
			a[i * n + j] = a[indMax * n + j];
			a[indMax * n + j] = tmp;
		}

		tmp = b[i];
		b[i] = b[indMax];
		b[indMax] = tmp;

		if (fabs(a[i * n + i]) < 1e-99) //проверка на вырожденность
			return -1;

		tmp = 1.0 / a[i * n + i]; //делим текущую строку на главный элемент
		for (j = i; j < n; j++){
			a[i * n + j] *= tmp;
		}
		b[i] *= tmp;

		for (j = 0; j < i; j++){ //для строк выше текущей
			tmp = a[j * n + i];
			for (k = i; k < n; k++){
				a[j * n + k] -= a[i * n + k] * tmp;
			}
			b[j] -= b[i] * tmp;
		}

		for (j = i + 1; j < n; j++){ //для строк ниже текущей
			tmp = a[j * n + i];
			for (k = i; k < n; k++){
				a[j * n + k] -= a[i * n + k] * tmp;
			}
			b[j] -= b[i] * tmp;
		}
	}

	for (i = n - 1; i >= 0; i--){ //запись решения в массив х
		x[i] = b[i];
	}
	return 0;
}
double canonPolinom(int deg, double var, double *coef){
    double result = 0;
    int i;
    for(i = 0; i < deg; i++){
        result += coef[i]*__pow(var, i);
    }
    return result;
}
//---------------------Многочлен Лагранжа-----------------------
double lagrange(int deg, double var, double *y, double *x){
    int i, j;
    double result = 0.0, temp = 1.0;
    for(i = 0; i < deg; i++){
        temp = 1.0;
        for(j = 0; j < deg; j++){
            if(j != i){
                temp *= (var - x[j])/(x[i] - x[j]);
            }
        }
        result += y[i]*temp;
    }
    return result;
}
//-------Дополнительные функции------
int __checkEqual(double *mas, int n){
    int i, j;
    for(i = 0; i < n; i++){
        for(j = i + 1; j < n; j++){
            if(fabs(mas[i] - mas[j]) < M_EPS)
                return -1;
        }
    }
    return 1;
}
double __pow(double x, int deg){
    int i;
    double temp = 1;
    for(i = 0; i < deg; i++){
        temp *= x;
    }
    return temp;
}