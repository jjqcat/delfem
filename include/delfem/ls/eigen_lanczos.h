
#if !defined(EIGEN_LANCZOS_H)
#define EIGEN_LANCZOS_H

#include <vector>

namespace LsSol{
class CLinearSystem;
class CPreconditioner;

//bool ArnoldiQR(Fem::Ls::CLinearSystem_Eigen& ls );
//bool EigenValue_Lanczos( unsigned int nlambda, std::vector<double>& aLambda, unsigned int num_iter, Fem::Ls::CLinearSystem_Eigen& ls);
double MinimumEigenValueVector_InvPower(
	LsSol::CLinearSystem& ls , 
	LsSol::CPreconditioner& pls, 
	const std::vector<unsigned int>& aIdVec, 
	unsigned int& iter_num);
}

#endif
