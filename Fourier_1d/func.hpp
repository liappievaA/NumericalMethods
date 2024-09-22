double u(double x);
int makePoints(int NUMBER_OF_DOTS, double *MASSIVE_OF_DOTS);
double basis_function(int m, int k, int N);
void write_bfm_for_all_points(int m, int N, double *ph);
int makeKoef_C_nk(int NUMBER_OF_DOTS, double *MASSIVE_OF_U, double *MASSIVE_C_nk, double *_trashMassive);
double seriesOfFurierAtPoint(int N, double *koef, double var);
double my_m_eps();
