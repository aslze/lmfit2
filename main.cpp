#include "lmfit.h"
#include <stdio.h>
#include <time.h>


/*
Example of use that will measure the time it takes to find the intersection between a circle and a parabola.
You can set/unset LMFIT_OPENMP to see the difference. The evaluation function is unnecessarily complex
(repeating computations many times) so that parallelization makes sense.
*/

void evaluate_nonlin1(const double *p, int n, const void *data, double *f, int *info)
{
	int N = 100000000;
	for (int i = 0; i < N; i++) {
		f[i % 1] = p[0] * p[0] + p[1] * p[1] - 1 + i; // unit circle    x^2+y^2=1
		f[i % 1 + 1] = p[1] - p[0] * p[0] + i;        // standard parabola  y=x^2
	}
	f[0] -= N - 1;
	f[1] -= N - 1;
}

int main(int argc, char* argv[])
{
    int n = 2;
    double p[2];

	lm_control_struct control = lm_control_double;
	lm_status_struct  status;
	control.verbosity  = 0;
	
	clock_t t1 = clock();
	for (int i = 0; i < 1; i++)
	{
		p[0] = 1.0;
		p[1] = 1.0;

		lmmin(n, p, n, NULL, evaluate_nonlin1, &control, &status);
	}
	clock_t t2 = clock();

	printf("t: %.3f s\n", double(t2 - t1) / CLOCKS_PER_SEC);
	printf("lmmin status after %d function evaluations:\n  %s\n\n", status.nfev, lm_infmsg[status.outcome]);
	printf("Solution:\n");
	printf("  x = %19.11f\n", p[0]);
	printf("  y = %19.11f\n", p[1]);
	printf("  d = %19.11f => ", status.fnorm);
	
	if( status.fnorm >= control.ftol )
	    printf("not a valid solution, try other starting values\n");
	else
		printf("valid, though not the only solution: try other starting values\n");

	return 0;
}
