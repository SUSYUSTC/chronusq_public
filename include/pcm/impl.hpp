#include <pcm.hpp>
#include <fields.hpp>
#include <basisset.hpp>
#include <molecule.hpp>
#include <chronusq_sys.hpp>


template<class MatsT>
void PCMBase::formpotential(CQMemManager& mem, MatsT* PDM, EMPerturbation& perb, BasisSet& basisset)
{
	//FIXME: memory leak
	Eigen::Matrix<MatsT,-1,1> Potential(this->grid_size);
	Eigen::Map<Eigen::Matrix<MatsT,-1,-1>> Matrix_PDM(PDM,this->nB,this->nB);
	Eigen::Matrix<MatsT,-1,-1> PDM_tran=Matrix_PDM.transpose();
	Eigen::Map<Eigen::Matrix<MatsT,-1,1>> DM(PDM_tran.data(),this->num_ele,1);
	if (this->store)
	{
		assert(this->is_stored);
		Potential=this->ints*DM;
	}
	else
	{
		for(int i=0;i!=this->grid_size;++i)
		{
			std::array<double,3> center={grid(0,i),grid(1,i),grid(2,i)};
			double* new_ints=PointFock(mem, perb,basisset,center);
			Eigen::Map<Eigen::RowVectorXd> V(new_ints,num_ele);
			Potential(i)=V*DM;
		}
	}
	//FIXME:A weird factor
	this->surp=convert_double(Potential)*8;
}
template<class MatsT>
MatsT* PCMBase::addFock(double* fock1, MatsT* fock2)
{
	Eigen::Map<Eigen::MatrixXd> Fock1(fock1,this->nB,this->nB);
	Eigen::Map<Eigen::Matrix<MatsT,-1,-1>> Fock2(fock2,this->nB,this->nB);
	Eigen::Matrix<MatsT,-1,-1> NewFock=Fock1+Fock2;
	return NewFock.data();
}


