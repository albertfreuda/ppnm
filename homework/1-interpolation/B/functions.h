typedef struct {int n; double *x,*y,*b,*c;} qspline;
qspline* qspline_alloc(int n, double *x, double *y);
void qspline_free(qspline *s);
double qspline_eval(qspline *s,double z);
double qspline_integ(qspline *s,double z);
double qspline_deriv(qspline *s,double z);
