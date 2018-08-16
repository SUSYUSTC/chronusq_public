/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2018 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */

#ifndef __PCM_HPP__
#define __PCM_HPP__
#include <PCMSolver/PCMInput.h>
#include <PCMSolver/pcmsolver.h>
#include <basisset.hpp>
#include <molecule.hpp>
#include <cxxapi/input.hpp>
#include <fields.hpp>
#include <iostream>
//#include <cxxapi/options.hpp>
#include <Eigen/Core>
namespace ChronusQ
{
	class PCMBase
	{
	public:
		bool use_PCM=false;
		bool store;//whether to store the matrix
		bool is_stored;
		PCMInput host_input;
		pcmsolver_context_t* pcm_context;
		size_t nB;
		size_t num_ele=nB;
		int grid_size;
		Eigen::Matrix3Xd grid;
		Eigen::VectorXd nucp;//nuclear potential
		Eigen::VectorXd surp;//surface potential
		Eigen::VectorXd surc;//surface charge
		Eigen::MatrixXd ints;
		Eigen::MatrixXd pcmfock;
		PCMBase();
		PCMBase(CQInputFile& input, BasisSet& basisset);
		void initialize(const Molecule& molecule);
		friend std::ostream& operator<< (std::ostream& os, const PCMBase pcmbase);
		double* PointFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset, std::array<double,3>& center);
		void formFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset);
		void storeFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset);
		Eigen::RowVectorXd convert_double(Eigen::RowVectorXcd vec);
		Eigen::RowVectorXd convert_double(Eigen::RowVectorXd vec);
		//template<class MatsT>
		//void formpotential(CQMemManager& mem, MatsT* PDM, EMPerturbation& perb, BasisSet& basisset);
		void formcharge();
		//template<class MatsT>
		//MatsT* addFock(double* fock1, MatsT* fock2);
		template<class MatsT>
		void formpotential(CQMemManager& mem, MatsT* PDM, EMPerturbation& perb, BasisSet& basisset)
		{
			//FIXME: memory leak
			Eigen::Matrix<MatsT,1,-1> Potential(1,this->grid_size);
			Eigen::Map<Eigen::Matrix<MatsT,-1,-1>> Matrix_PDM(PDM,this->nB,this->nB);
			Eigen::Matrix<MatsT,-1,-1> PDM_tran=Matrix_PDM.transpose();
			Eigen::Map<Eigen::Matrix<MatsT,1,-1>> DM(PDM_tran.data(),1,this->num_ele);
			std::cout << "Matrices created" << std::endl;
			if (this->store)
			{
				std::cout << "Store detected" << std::endl;
				assert(this->is_stored);
				Potential=DM*this->ints;
				std::cout << "Potential calculated" << std::endl;
			}
			else
			{
				std::cout << "No store detected" << std::endl;
				for(int i=0;i!=this->grid_size;++i)
				{
					std::cout << " Step " << i;
					std::array<double,3> center={grid(0,i),grid(1,i),grid(2,i)};
					double* new_ints=PointFock(mem, perb,basisset,center);
					Eigen::Map<Eigen::VectorXd> V(new_ints,num_ele);
					Potential(i)=DM*V;
				}
				std::cout << "Potential calculated" << std::endl;
			}
			this->surp=this->convert_double(Potential).transpose()+this->nucp;
		}
		template<class MatsT>
		void addFock(MatsT* fock)
		{
			Eigen::Map<Eigen::Matrix<MatsT,-1,-1>> Fock(fock,this->nB,this->nB);
			std::cout << "Matrices created" << std::endl;
			Eigen::Matrix<MatsT,-1,-1> NewFock=this->pcmfock+Fock;
			std::cout << "Addition finished" << std::endl;
			fock=NewFock.data();;
		}
	};
}


#endif
