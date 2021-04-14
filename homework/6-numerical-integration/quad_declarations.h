double integral(double f(double),
		double a,
		double b,
		double abstol,
		double reltol);

double quad(double f(double),
		double a,
		double b,
		double f2,
		double f4,
		double abstol,
		double reltol,
		int nrec);


double clenshaw_curtis(double f(double),
		double a,
		double b,
		double abstol,
		double reltol);
