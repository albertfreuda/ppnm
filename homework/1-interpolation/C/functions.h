typedef struct {int n; double *x,*y,*b,*c,*d;} cspline;
cspline* cspline_alloc(int n, double *x, double *y);
void cspline_free(cspline *s);
double cspline_eval(cspline *s,double z);
double cspline_integ(cspline *s,double z);
double cspline_deriv(cspline *s,double z);
