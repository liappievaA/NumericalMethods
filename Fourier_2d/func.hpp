double u(double x, double y);
int makePoints(int NUMBER_OF_DOTS, double *MASSIVE_OF_DOTS);
double seriesOfFurierAtPoint(int N, double *koef, double var_x, double var_y);
int make_C_mn(int NUMBER_OF_DOTS, double *MASSIVE_OF_DOTS, 
                        double *_trashMatrix, double *_trashMassive_1, double *_trashMassive_2, 
                        double *MATRIX_OF_U, double *MATRIX_C_nk);
int makeKoef_C_nk(int NUMBER_OF_DOTS, double *MASSIVE_OF_U, double *MASSIVE_OF_DOTS, double *MASSIVE_C_nk, double *_trashMassive);
