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
		bool store=false;//whether to store the matrix from surface charge to fock matrix, which is also the matrix from density matrix to electronic potential
		bool is_stored=false;
		std::string int_path;
		PCMInput host_input;//Interface of pcmsolver
		pcmsolver_context_t* pcm_context;//Interface of pcmsolver
		size_t DebugLevel=0;
		size_t DebugDepth=0;
		size_t nB;
		size_t num_ele=nB;//number of elements, seems that it's wrong, but this variable will be updated during initialize
		bool start_save=false;
		size_t times=0;//when times%savetep==0, surface potential and surface charge will be saved
		size_t savestep=100;//default is to save every 100 steps
		int grid_size;//num of grids
		Eigen::Matrix3Xd grid;//coordinates of grids
		Eigen::VectorXd nucp;//surface potential contributed by nuclei
		Eigen::VectorXd surp;//surface potential
		Eigen::VectorXd surc;//surface charge
		double* ints;//the matrix from potential to fock matrix
		double* pcmfock;
		PCMBase();
		PCMBase(CQInputFile& input, BasisSet& basisset);
		void initialize(CQMemManager& mem,const Molecule& molecule);
		friend std::ostream& operator<< (std::ostream& os, const PCMBase pcmbase);
		//calculate the fock matrix contributed by point charge
		double* PointFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset, std::array<double,3>& center);
		//calcualte the pcm fock matrix
		void formFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset);
		//calculate the matrix from surface potential to fock matrix
		void storeFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset);
		Eigen::RowVectorXd convert_double(Eigen::RowVectorXcd vec);
		Eigen::RowVectorXd convert_double(Eigen::RowVectorXd vec);
		//calculate the surface charge
		void formcharge();
		//calculate the surface potential
		template<class MatsT>
		void formpotential(CQMemManager& mem, MatsT* PDM, EMPerturbation& perb, BasisSet& basisset);
		//add pcm fock to total fock
		template<class MatsT>
		void addFock(MatsT* fock);
		//calcualte the energy contributed by pcm
		double computeEnergy();
	};
}


#endif
