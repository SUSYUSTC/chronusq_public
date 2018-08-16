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


#include <pcm.hpp>
#include <chronusq_sys.hpp>
#include <libint2.hpp>
#include <aointegrals.hpp>
#include <fields.hpp>
#include <basisset.hpp>
#include <molecule.hpp>
#include <iostream>
#include <cstring>
#include <cnpy.h>
#include <cxxapi/options.hpp>


namespace ChronusQ
{
	std::ostream& operator<< (std::ostream& os, const PCMBase pcmbase)
	{
		os << "Use PCM?\t" << pcmbase.use_PCM << std::endl;
		os << "Store all the integrations?\t" << pcmbase.store << std::endl;
		os << "Is integrations calculated?\t" << pcmbase.is_stored << std::endl;
		os << "Number of grids:\t" << pcmbase.grid_size << std::endl;
		os << "Number of basis:\t" << pcmbase.nB << std::endl;
		if (pcmbase.grid.size()!=0 && pcmbase.surc.size()!=0)
		{
			os << "Coordinates of grids and surface charge:" << std::endl;
			for(int i=0;i!=pcmbase.grid_size;++i)
			{
				os << pcmbase.grid.col(i);
				os << pcmbase.surc[i];
				os << std::endl;
			}
		}
		else if (pcmbase.grid.size()!=0)
		{
			os << "Surface charge has not been calculated" << std::endl;
			os << "Coordinates of grids:" << std::endl;
			os << pcmbase.grid.transpose() << std::endl;
		}
		else if (pcmbase.surc.size()!=0)
		{
			std::cerr << "Hey stupid guy, your program curshes" << std::endl;
		}
		else
		{
			os << "Grids are not calculated" << std::endl;
		}
		static bool is_pcmstore_printed=false;
		if (pcmbase.ints.size()!=0 and not is_pcmstore_printed)
		{
			os << "Saving information of integrations to pcmint.npy" << std::endl;
			std::vector<size_t> size={size_t(pcmbase.grid_size),pcmbase.nB,pcmbase.nB};
			cnpy::npy_save("pcmint.npy",pcmbase.ints.data(),size,"w");
			os << "Saving finished" << std::endl;
			is_pcmstore_printed=true;
		}
		else if(is_pcmstore_printed)
		{
			os << "Information of integrations is stored, won't store again" << std::endl;
		}
		else
		{
			os << "Information of integrations will be calculated each step" << std::endl;
		}
		return os;
	}
	PCMBase::PCMBase()
	{
		this->use_PCM=false;
	}
	PCMBase::PCMBase(CQInputFile& input, BasisSet& basisset)
	{
		try
		{
			double area=input.getData<double>("PCM.AREA");
			std::string solver_type=input.getData<std::string>("PCM.SOLVER_TYPE");
			std::string solvent=input.getData<std::string>("PCM.SOLVENT");
			this->store=input.getData<bool>("PCM.STORE");
			this->host_input.area=area;
			std::strcpy(this->host_input.solver_type,solver_type.c_str());
			std::strcpy(this->host_input.solvent,solvent.c_str());
			std::strcpy(this->host_input.cavity_type, "Gepol");
			std::strcpy(this->host_input.radii_set, "Bondi");
			this->nB=basisset.nBasis;
			this->num_ele=nB*nB;
			this->host_input.scaling=true;
			this->use_PCM=true;
			this->grid_size=0;
		}
		catch(...)
		{
			std::cout << "Wrong type of PCM, PCM won't be used" << std::endl;
			this->use_PCM=false;
		}
	}
	//determine num of grids and locations of grids
	void PCMBase::initialize(const Molecule& molecule)
	{
		int nAtoms=molecule.nAtoms;
		double* charges=new double[nAtoms];
		for(int i=0;i!=nAtoms;++i)
		{
			charges[i]=double(molecule.atoms[i].atomicNumber);
		}
		double* coordinates=new double[3*nAtoms];
		for(int i=0;i!=nAtoms;++i)
		{
			for(int j=0;j!=3;++j)
				coordinates[3*i+j]=molecule.atoms[i].coord[j];
		}
		int symmetry_info[1]={0};
		this->pcm_context=pcmsolver_new(PCMSOLVER_READER_HOST,nAtoms,charges,coordinates,symmetry_info,&(this->host_input),[](const char* message){;});
		this->grid_size= pcmsolver_get_cavity_size(this->pcm_context); 
		this->grid = Eigen::Matrix3Xd(3,grid_size);
		pcmsolver_get_centers(pcm_context, grid.data());
	}
	double* PCMBase::PointFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset, std::array<double,3>& center)
	{
		std::cout << "PointFock started, ";
		Atom atom("H",center);
		atom.atomicNumber=1;
		std::vector<Atom> atoms(1);
		atoms[0]=atom;
		Molecule molecule(0,2,atoms);
		std::cout << "Molecule created, ";
		std::array<double,3> magAmp = perb.getDipoleAmp(Magnetic);
		std::vector<double*> Potential;
		AOIntegrals<double> aoints(mem, molecule, basisset);
		std::cout << "AOIntegrals created, ";
		Potential = aoints.OneEDriverLibint(libint2::Operator::nuclear,basisset.shells);
		std::cout << "Potential calculated, ";
		double* V=reinterpret_cast<double*>(Potential[0]);
		//this is because the structure of DM
		Eigen::Map<Eigen::VectorXd> EigenV(V,this->num_ele);
		EigenV*=2;
		std::cout << "Times Potential by 2";
		return EigenV.data();
	}
	Eigen::RowVectorXd PCMBase::convert_double(Eigen::RowVectorXd vec)
	{
		return vec;
	}
	Eigen::RowVectorXd PCMBase::convert_double(Eigen::RowVectorXcd vec)
	{
		return vec.real();
	}
	//input surface potential, store surface charge inside
	void PCMBase::formcharge()//surface potential
	{
		pcmsolver_set_surface_function(this->pcm_context, this->grid_size, this->surp.data(), "mep");
		std::cout << "Potential set" << std::endl;
		pcmsolver_compute_asc(pcm_context,"mep","asc",0);
		std::cout << "Charge calculated" << std::endl;
		double* surcdata=new double[this->grid_size];
		pcmsolver_get_surface_function(this->pcm_context, this->grid_size, surcdata, "asc");
		this->surc=Eigen::Map<Eigen::VectorXd>(surcdata,this->grid_size,1);
		delete surcdata;
		std::cout << "Charge stored" << std::endl;
	}
	//calculate Fock matrix generated by each grid and store it
	void PCMBase::storeFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset)
	{
		this->ints=Eigen::MatrixXd(this->num_ele, this->grid_size);
		for(int i=0;i!=this->grid_size;++i)
		{
			std::array<double,3> center={grid(0,i),grid(1,i),grid(2,i)};
			double* Vp=PointFock(mem, perb, basisset, center);
			//FIXME: memory leak of V
			Eigen::Map<Eigen::VectorXd> V(Vp,this->num_ele);
			ints.col(i)=V;
		}
		this->is_stored=true;
	}
	//TODO :optimize this part and make it parallel
	//If Fock matrix is stored, sum them up, else sum up one by one
	void PCMBase::formFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset)
	{
		Eigen::VectorXd sum_Fock(this->num_ele);
		if(this->store)
		{
			std::cout << "Store detected" << std::endl;
			assert(this->is_stored);
			sum_Fock=this->ints*this->surc;
			this->pcmfock=Eigen::Map<Eigen::MatrixXd>(sum_Fock.data(),this->nB,this->nB);
			std::cout << "Fock calculated" << std::endl;
		}
		else
		{
			std::cout << "No store detected" << std::endl;
			sum_Fock=Eigen::VectorXd::Zero(num_ele);
			for(int i=0;i!=this->grid_size;++i)
			{
				std::array<double,3> center={grid(0,i),grid(1,i),grid(2,i)};
				//get Fock matrix at each grid
				double* new_ints=PointFock(mem, perb,basisset,center);
				Eigen::Map<Eigen::VectorXd> V(new_ints,num_ele);
				sum_Fock+=V;
			}
			this->pcmfock=Eigen::Map<Eigen::MatrixXd>(sum_Fock.data(),this->nB,this->nB);
			std::cout << "Fock calculated" << std::endl;
		}
	}
}

