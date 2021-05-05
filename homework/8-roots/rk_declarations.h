
void rkstep12(
	void f(double t,gsl_vector* y,gsl_vector* dydt,double param),//dydt=f(t,y)
	double t,                                       //current value	t
	gsl_vector* yt,                                 //current value y(t)
	double h,                                       //step size h
	gsl_vector* yh,					//out: y(t+h)
	gsl_vector* err,
	double param);					//out: error estimate dely

void driver(void f(double t, gsl_vector* y,gsl_vector* dydt,double param),
		double a,//t0
		gsl_vector* ya, //y(a)
		double b,//end point of integration
		gsl_vector* yb, //y(b) to be calc
		gsl_vector* err, //err memory location
		double h,//initial step size
		double acc,//abs acc. goal
		double eps,
		double param);//relative acc. goal

void driver2(void f(double t, gsl_vector* y,gsl_vector* dydt,double param),
		double a,//t0
		gsl_vector* ya, //y(a)
		double b,//end point of integration
		gsl_vector* yb, //y(b) to be calc
		gsl_vector* err, //err memory location
		double h,//initial step size
		double acc,//abs acc. goal
		double eps,
		double param);//relative acc. goal

void print_vector(gsl_vector* v);

void print_matrix(gsl_matrix* A);
