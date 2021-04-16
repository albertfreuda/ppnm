double integral(double f(double),
		double a,
		double b,
		double abstol,
		double reltol,
		int* n_evaluations);

double quad(double f(double),
		double a,
		double b,
		double f2,
		double f4,
		double abstol,
		double reltol,
		int nrec,
		int* n_evaluations);


double clenshaw_curtis(double f(double),
		double a,
		double b,
		double abstol,
		double reltol,
		int* n_evaluations);


double quadC(double f(double),
		double a,
		double b,
		double f2,
		double f4,
		double abstol,
		double reltol,
		int nrec,
		double* err);

double integralC(double f(double),
		double a,
		double b,
		double abstol,
		double reltol,
		double* err);
