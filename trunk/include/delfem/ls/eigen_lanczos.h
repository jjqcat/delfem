
#if !defined(EIGEN_LANCZOS_H)
#define EIGEN_LANCZOS_H

#include <vector>

namespace Fem{
namespace Ls{
	class CLinearSystem_Eigen;
	class CPreconditioner;
}

namespace Sol
{

	bool ArnoldiQR(Fem::Ls::CLinearSystem_Eigen& ls );
	bool EigenValue_Lanczos( unsigned int nlambda, std::vector<double>& aLambda, unsigned int num_iter, Fem::Ls::CLinearSystem_Eigen& ls);
	double MinimumEigenValueVector_InvPower(Fem::Ls::CLinearSystem_Eigen& ls , 
		Fem::Ls::CPreconditioner& pls, 
		const std::vector<unsigned int>& aIdVec, 
		unsigned int& iter_num);
}
}

#endif
