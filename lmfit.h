/*
Copywight: aslze (2019)
Copyright: Joachim Wuttke, Forschungszentrum Juelich GmbH (2004-2013)
*/

#ifndef LMMIN2_H
#define LMMIN2_H

#ifdef __cplusplus
extern "C" {
#endif

/**
Collection of input parameters for fit control.
*/
typedef struct
{
	double ftol;      /* Relative error desired in the sum of squares.
                         Termination occurs when both the actual and
                         predicted relative reductions in the sum of squares
                         are at most ftol. */
	double xtol;      /* Relative error between last two approximations.
                         Termination occurs when the relative error between
                         two consecutive iterates is at most xtol. */
	double gtol;      /* Orthogonality desired between fvec and its derivs.
                         Termination occurs when the cosine of the angle
                         between fvec and any column of the Jacobian is at
                         most gtol in absolute value. */
	double epsilon;   /* Step used to calculate the Jacobian, should be
                         slightly larger than the relative error in the
                         user-supplied functions. */
	double stepbound; /* Used in determining the initial step bound. This
                         bound is set to the product of stepbound and the
                         Euclidean norm of diag*x if nonzero, or else to
                         stepbound itself. In most cases stepbound should lie
                         in the interval (0.1,100.0). Generally, the value
                         100.0 is recommended. */
	int patience;     /* Used to set the maximum number of function evaluations
                         to patience*(number_of_parameters+1). */
	int scale_diag;   /* If 1, the variables will be rescaled internally.
                         Recommended value is 1. */
	void* msgfile;    /* Actually FILE*. Progress messages will be written to this file. */
	int verbosity;    /* OR'ed: 1: print some messages; 2: print Jacobian. */
	int n_maxpri;     /* -1, or max number of parameters to print. */
	int m_maxpri;     /* -1, or max number of residuals to print. */
} lm_control_struct;

/**
Collection of output parameters for status info.
*/
typedef struct
{
	double fnorm;  /* norm of the residue vector fvec. */
	int nfev;      /* actual number of iterations. */
	int outcome;   /* Status indicator. Nonnegative values are used as index
                      for the message text lm_infmsg, set in lmmin.c. */
	int userbreak; /* Set when function evaluation requests termination. */
} lm_status_struct;

/* Preset (and recommended) control parameter settings. */
extern const lm_control_struct lm_control_double;
extern const lm_control_struct lm_control_float;

/* Preset message texts. */
extern const char* lm_infmsg[];
extern const char* lm_shortmsg[];

/**
Levenberg-Marquardt minimization.

It minimizes the sum of the squares of m nonlinear functions
in n variables by a modified Levenberg-Marquardt algorithm.
The function evaluation is done by the user-provided routine 'evaluate'.
The Jacobian is then calculated by a forward-difference approximation.

\param n_par is the number of variables (INPUT, positive integer).

\param par is the solution vector (INPUT/OUTPUT, array of length n).
    On input it must be set to an estimated solution.
    On output it yields the final estimate of the solution.

\param m_dat is the number of functions to be minimized (INPUT, positive integer).
    It must fulfill m>=n.

\param data is a pointer that is ignored by lmmin; it is however forwarded
    to the user-supplied functions evaluate and printout.
    In a typical application, it contains experimental data to be fitted.

\param evaluate is a user-supplied function that calculates the m functions. Parameters:
      - n, x, m, data as above.
      - fvec is an array of length m; on OUTPUT, it must contain the
        m function values for the parameter vector x.
      - userbreak is an integer pointer. When userbreak is set to a
        nonzero value, lmmin will terminate.

\param control contains INPUT variables that control the fit algorithm,
    as declared and explained in lmstruct.h

\param status contains OUTPUT variables that inform about the fit result,
    as declared and explained in lmstruct.h
*/
void lmmin(const int n_par, double* par, const int m_dat, const void* data,
           void (*evaluate)(const double* par, const int m_dat, const void* data, double* fvec, int* userbreak),
           const lm_control_struct* control, lm_status_struct* status);

/**
Refined calculation of Eucledian norm.
*/
double lm_enorm( const int, const double* );

/**
A simplified API for curve fitting using the generic Levenberg-Marquardt routine lmmin.
*/
void lmcurve(const int n_par, double *par, const int m_dat, const double *t, const double *y,
             double (*f)(const double t, const double *par),
             const lm_control_struct *control, lm_status_struct *status);

/**
Implements lmcurve_tyd(), a variant of lmcurve() that weighs
data points y(t) with the inverse of the standard deviations dy.
*/
void lmcurve_tyd(const int n_par, double* par, const int m_dat,
                 const double* t, const double* y, const double* dy,
                 double (*f)(double t, const double* par),
                 const lm_control_struct* control, lm_status_struct* status);

#ifdef __cplusplus
}
#endif

#endif
