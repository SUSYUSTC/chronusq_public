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
		bool use_PCM;
		bool store;//whether to store the matrix
		bool is_stored;
		PCMInput host_input;
		pcmsolver_context_t* pcm_context;
		size_t nB;
		size_t num_ele=nB;
		int grid_size;
		Eigen::Matrix3Xd grid;
		Eigen::VectorXd surp;
		Eigen::VectorXd surc;//surface charge
		Eigen::MatrixXd ints;
		PCMBase();
		PCMBase(CQInputFile& input, BasisSet& basisset);
		void initialize(const Molecule& molecule);
		friend std::ostream& operator<< (std::ostream& os, const PCMBase pcmbase);
		double* PointFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset, std::array<double,3>& center);
		double* formFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset);
		void storeFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset);
		//template<class MatsT>
		//void formpotential(CQMemManager& mem, MatsT* PDM, EMPerturbation& perb, BasisSet& basisset);
		void formcharge();
		//template<class MatsT>
		//MatsT* addFock(double* fock1, MatsT* fock2);
		template<class MatsT>
		void formpotential(CQMemManager& mem, MatsT* PDM, EMPerturbation& perb, BasisSet& basisset)
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
		MatsT* addFock(double* fock1, MatsT* fock2)
		{
			Eigen::Map<Eigen::MatrixXd> Fock1(fock1,this->nB,this->nB);
			Eigen::Map<Eigen::Matrix<MatsT,-1,-1>> Fock2(fock2,this->nB,this->nB);
			Eigen::Matrix<MatsT,-1,-1> NewFock=Fock1+Fock2;
			return NewFock.data();
		}
	};
}


#endif
