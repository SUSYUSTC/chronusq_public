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
#include <singleslater.hpp>
#include <cxxapi/options.hpp>
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
		double* PointFock(std::shared_ptr<CQMemManager> mem, EMPerturbation& perb, BasisSet& basisset, std::array<double,3>& center);
		double* formFock(std::shared_ptr<CQMemManager> mem, EMPerturbation& perb, BasisSet& basisset);
		void storeFock(std::shared_ptr<CQMemManager> mem, EMPerturbation& perb, BasisSet& basisset);
		template<class MatsT>
		void formpotential(std::shared_ptr<CQMemManager> mem, MatsT* PDM, EMPerturbation& perb, BasisSet& basisset);
		void formcharge();
		template<class MatsT>
		MatsT* addFock(double* fock1, MatsT* fock2);
	};
	/*
	double* convert_double(double* array);
	double* convert_double(dcomplex* array);
	*/
}


#endif
