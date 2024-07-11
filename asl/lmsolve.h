// Copyright(c) 2023-2024 aslze
// Licensed under the MIT License (http://opensource.org/licenses/MIT)

#pragma once
#include <asl/Matrix.h>
#include "lmfit.h"

namespace asl
{
template<class F>
struct LMData
{
	F   f;
	int n;
};

template<class T, class F>
Matrix_<T> solveZeroLM(F f, const Matrix_<T>& x0, const SolveParams& p = SolveParams(100))
{
	auto evaluate = [](const double* p, int n, const void* data, double* f, int* info) {
		Matrixd x(((LMData<F>*)data)->n, 1, p);
		auto    r = ((LMData<F>*)data)->f(x);
		for (int i = 0; i < n; i++)
			f[i] = r[i];
	};

	lm_control_struct control = lm_control_double;
	lm_status_struct  status;
	control.verbosity = p.maxiter < 0 ? 1 : 0;
	control.patience = p.maxiter;
	Matrixd   x = x0.clone();
	LMData<F> data{ f, x.rows() };
	lmmin(x0.rows(), &x(0, 0), f(x0).rows(), &data, evaluate, &control, &status);
	return x;
}

#ifdef ASL_HAVE_INITLIST
template<class T, class F>
Matrix_<T> solveZeroLM(F f, const std::initializer_list<T>& x0, const SolveParams& p = SolveParams(100))
{
	return solveZeroLM(f, Matrix_<T>(x0), p);
}
#endif

}
