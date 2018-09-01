#include <Eigen/Core>
#include <singleslater.hpp>
namespace ChronusQ
{
	template<>
	void SingleSlater<double,double>::fliporbits()
	{
		if(this->nC==1)
		{
			size_t NB = aoints->basisSet().nBasis;
			Eigen::MatrixXd AMO=Eigen::Map<Eigen::Matrix<double,-1,-1,Eigen::RowMajor>>(this->mo1,NB,NB);
			for(int i=0;i<3;++i)
				std::cout << AMO(0,i) << std::endl;
			Eigen::MatrixXd BMO=Eigen::Map<Eigen::Matrix<double,-1,-1,Eigen::RowMajor>(this->mo2,NB,NB);
			Eigen::MatrixXd contracted_AMO(this->nOA,NB);
			contracted_AMO=AMO.block(0,0,this->nOA,NB);
			contracted_AMO.row(nOA)=AMO.row(nOA+1);
			Eigen::MatrixXd contracted_BMO(this->nOB,NB);
			contracted_BMO=BMO.block(0,0,this->nOB,NB);
			Eigen::MatrixXd ADM=contracted_AMO.transpose()*contracted_AMO;
			Eigen::MatrixXd BDM=contracted_BMO.transpose()*contracted_BMO;
			Eigen::Map<Eigen::MatrixXd> PDMs(this->onePDM[0],NB,NB);
			Eigen::Map<Eigen::MatrixXd> PDMz(this->onePDM[1],NB,NB);
			Eigen::Map<Eigen::MatrixXd> PDMorthos(this->onePDMOrtho[0],NB,NB);
			Eigen::Map<Eigen::MatrixXd> PDMorthoz(this->onePDMOrtho[1],NB,NB);
			Eigen::Map<Eigen::Matrix<double,-1,-1,Eigen::RowMajor>> ortho(this->ortho2,NB,NB);
			PDMs=ADM+BDM;
			PDMz=ADM-BDM;
			PDMorthos=ortho.transpose()*PDMs*ortho;
			PDMorthoz=ortho.transpose()*PDMz*ortho;
		}
	}
	template void SingleSlater<double,double>::fliporbits();
}
