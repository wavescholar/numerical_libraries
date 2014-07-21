// TestArpackCPPLib.cpp : Defines the entry point for the console application.

#include "stdafx.h"

void solve_a_real_symmetric_standard_eigenvalue_problem_in_regular_mode_using_the_ARluSymStdEig_class();
void real_symmetric_standard_eigenvalue_problem_in_shift_and_invert_mode_using_ARluSymStdEig();
void real_symmetric_dense_standard_eigenvalue_problem_in_regular_mode();
int _tmain(int argc, _TCHAR* argv[])
{
	real_symmetric_dense_standard_eigenvalue_problem_in_regular_mode();
	solve_a_real_symmetric_standard_eigenvalue_problem_in_regular_mode_using_the_ARluSymStdEig_class();
	real_symmetric_standard_eigenvalue_problem_in_shift_and_invert_mode_using_ARluSymStdEig();

	return 0;
}

