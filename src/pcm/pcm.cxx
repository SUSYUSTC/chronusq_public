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
#include <sjc_debug.hpp>


namespace ChronusQ
{
	std::ostream& operator<< (std::ostream& os, const PCMBase pcmbase)
	{
		os << "Use PCM?\t" << std::boolalpha << pcmbase.use_PCM << std::endl;
		os << "Store all the integrations?\t" << std::boolalpha << pcmbase.store << std::endl;
		os << "Is integrations calculated?\t" << std::boolalpha << pcmbase.is_stored << std::endl;
		os << "Number of grids:\t" << pcmbase.grid_size << std::endl;
		os << "Number of basis:\t" << pcmbase.nB << std::endl;
		os << "Save data every " << pcmbase.savestep << " steps" << std::endl;
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
			//All of these must be specified
			double area=input.getData<double>("PCM.AREA");
			std::string solver_type=input.getData<std::string>("PCM.SOLVER_TYPE");
			std::string solvent=input.getData<std::string>("PCM.SOLVENT");
			this->store=input.getData<bool>("PCM.STORE");
			this->host_input.outside_epsilon=input.getData<double>("PCM.EPSILON");
			this->savestep=input.getData<size_t>("QM.SAVESTEP");
			this->host_input.area=area;
			std::strcpy(this->host_input.solver_type,solver_type.c_str());
			std::strcpy(this->host_input.solvent,solvent.c_str());
			//The following four keywords is fixed
			std::strcpy(this->host_input.cavity_type, "Gepol");
			std::strcpy(this->host_input.radii_set, "UFF");
			std::strcpy(this->host_input.inside_type, "Vacuum");
			std::strcpy(this->host_input.outside_type, "UniformDielectric");
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
	//determine num of grids, locations of grids and the matrix from surface potential to surface charge
	void PCMBase::initialize(CQMemManager& mem,const Molecule& molecule)
	{
		int nAtoms=molecule.nAtoms;
		double* charges=mem.malloc<double>(nAtoms);
		for(int i=0;i!=nAtoms;++i)
		{
			charges[i]=double(molecule.atoms[i].atomicNumber);
		}
		double* coordinates=mem.malloc<double>(3*nAtoms);
		for(int i=0;i!=nAtoms;++i)
		{
			for(int j=0;j!=3;++j)
				coordinates[3*i+j]=molecule.atoms[i].coord[j];
		}
		int symmetry_info[1]={0};
		this->pcm_context=pcmsolver_new(PCMSOLVER_READER_HOST,nAtoms,charges,coordinates,symmetry_info,&(this->host_input),[](const char* message){std::cout << message << std::endl;});
		this->grid_size= pcmsolver_get_cavity_size(this->pcm_context); 
		this->grid = Eigen::Matrix3Xd(3,grid_size);
		std::vector<size_t> npy_size={(size_t)grid_size,3};
		pcmsolver_get_centers(pcm_context, grid.data());
		this->nucp=Eigen::VectorXd::Zero(this->grid_size,1);
		Eigen::Map<Eigen::Matrix3Xd> coor(coordinates,3,nAtoms);
		for(int i=0;i!=this->grid_size;++i)
		{
			for(int j=0;j!=nAtoms;++j)
			{
				Eigen::Vector3d vec=this->grid.col(i)-coor.col(j);
				nucp(i)+=charges[j]/vec.norm();
			}
		}
		for(int i=0;i!=this->grid_size;++i)
		this->pcmfock=mem.malloc<double>(this->nB*this->nB);
		mem.free(charges);
		mem.free(coordinates);
		Eigen::Matrix3Xd grid_angstrom=grid*0.52917721067;
		cnpy::npy_save("grid.npy",grid_angstrom.data(),npy_size,"w");
	}
	//perb is not used
	double* PCMBase::PointFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset, std::array<double,3>& center)
	{
		//std::cout << "PointFock started, ";
		//set the charge to 1
		Atom atom("H",center);
		atom.atomicNumber=1;
		std::vector<Atom> atoms(1);
		atoms[0]=atom;
		Molecule molecule(0,2,atoms);
		//std::cout << "Molecule created, ";
		std::array<double,3> magAmp = perb.getDipoleAmp(Magnetic);
		std::vector<double*> Potential;
		AOIntegrals<double> aoints(mem, molecule, basisset);
		//std::cout << "AOIntegrals created, ";
		Potential = aoints.OneEDriverLibint(libint2::Operator::nuclear,basisset.shells);
		//std::cout << "Potential calculated, ";
		double* V=reinterpret_cast<double*>(Potential[0]);
		for(int i=1;i!=int(Potential.size());++i)
			mem.free(Potential[i]);
		return V;
	}
	Eigen::RowVectorXd PCMBase::convert_double(Eigen::RowVectorXd vec)
	{
		return vec;
	}
	Eigen::RowVectorXd PCMBase::convert_double(Eigen::RowVectorXcd vec)
	{
		return vec.real();
	}
	//input surface potential, store surface charge inside pcmbase
	void PCMBase::formcharge()//surface potential
	{
		pcmsolver_set_surface_function(this->pcm_context, this->grid_size, this->surp.data(), "mep");
		//std::cout << "Potential set" << std::endl;
		pcmsolver_compute_asc(pcm_context,"mep","asc",0);
		//std::cout << "Charge calculated" << std::endl;
		double* surcdata=new double[this->grid_size];
		pcmsolver_get_surface_function(this->pcm_context, this->grid_size, surcdata, "asc");
		this->surc=Eigen::Map<Eigen::VectorXd>(surcdata,this->grid_size,1);
		if(this->start_save and this->times%this->savestep==0)
		{
			std::vector<size_t> npy_size={size_t(this->grid_size)};
			cnpy::npy_save("asc_"+std::to_string(this->times)+".npy",this->surc.data(),npy_size,"w");
			if (this->DebugLevel>=0)
				sjc_debug::debug0(this->DebugDepth,"Charge stored");
		}
	}
	//calculate Fock matrix generated by each grid and store it
	void PCMBase::storeFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset)
	{
		//std::cout << "Start store" << std::endl;
		this->ints=mem.malloc<double>(this->num_ele*this->grid_size);
		Eigen::Map<Eigen::MatrixXd> matints(this->ints,this->num_ele, this->grid_size);
		for(int i=0;i!=this->grid_size;++i)
		{
			std::array<double,3> center={grid(0,i),grid(1,i),grid(2,i)};
			double* Vp=PointFock(mem, perb, basisset, center);
			//FIXME: memory leak of V
			Eigen::Map<Eigen::VectorXd> V(Vp,this->num_ele);
			matints.col(i)=V;
			mem.free(Vp);
		}
		std::vector<size_t> npy_size={size_t(this->grid_size),this->nB,this->nB};
		cnpy::npy_save("pcmints.npy",this->ints,npy_size,"w");
		this->is_stored=true;
	}
	//If Fock matrix is stored, sum them up, else sum up one by one
	void PCMBase::formFock(CQMemManager& mem, EMPerturbation& perb, BasisSet& basisset)
	{
		Eigen::VectorXd sum_Fock(this->num_ele);
		Eigen::Map<Eigen::MatrixXd> matints(this->ints,this->num_ele,this->grid_size);
		Eigen::Map<Eigen::MatrixXd> matpcmfock(this->pcmfock,this->nB,this->nB);
		if(this->store)
		{
			//std::cout << "Store detected" << std::endl;
			assert(this->is_stored);
			sum_Fock=matints*this->surc;
		}
		else
		{
			//std::cout << "No store detected" << std::endl;
			sum_Fock=Eigen::VectorXd::Zero(num_ele);
			for(int i=0;i!=this->grid_size;++i)
			{
				std::array<double,3> center={grid(0,i),grid(1,i),grid(2,i)};
				//get Fock matrix at each grid
				double* new_ints=PointFock(mem, perb,basisset,center);
				Eigen::Map<Eigen::VectorXd> V(new_ints,num_ele);
				sum_Fock+=V;
			}
		}
		sum_Fock*=2;
		matpcmfock=Eigen::Map<Eigen::MatrixXd>(sum_Fock.data(),this->nB,this->nB);
		//std::cout << "Fock calculated" << std::endl;
	}
	template<class MatsT>
	void PCMBase::formpotential(CQMemManager& mem, MatsT* PDM, EMPerturbation& perb, BasisSet& basisset)
	{
		//FIXME: memory leak
		Eigen::Matrix<MatsT,1,-1> Potential(1,this->grid_size);
		Eigen::Map<Eigen::Matrix<MatsT,-1,-1>> Matrix_PDM(PDM,this->nB,this->nB);
		Eigen::Matrix<MatsT,-1,-1> PDM_tran=Matrix_PDM.transpose();
		Eigen::Map<Eigen::Matrix<MatsT,1,-1>> DM(PDM_tran.data(),1,this->num_ele);
		Eigen::Map<Eigen::MatrixXd> matints(this->ints,this->num_ele,this->grid_size);
		//std::cout << "Matrices created" << std::endl;
		if (this->store)
		{
			//std::cout << "Store detected" << std::endl;
			assert(this->is_stored);
			Potential=DM*matints;
			//std::cout << "Potential calculated" << std::endl;
		}
		else
		{
			//std::cout << "No store detected" << std::endl;
			for(int i=0;i!=this->grid_size;++i)
			{
				//if(i%50==0)
					//std::cout << "Step " << i;
				std::array<double,3> center={grid(0,i),grid(1,i),grid(2,i)};
				double* new_ints=PointFock(mem, perb,basisset,center);
				Eigen::Map<Eigen::VectorXd> V(new_ints,num_ele);
				Potential(i)=DM*V;
			}
			//std::cout << "Potential calculated" << std::endl;
		}
		this->surp=this->convert_double(Potential).transpose()+this->nucp;
		if (this->start_save and this->times%this->savestep==0)
		{
			std::vector<size_t> npy_size={size_t(this->grid_size)};
			cnpy::npy_save("potential_"+std::to_string(this->times)+".npy",this->surp.data(),npy_size,"w");
			if (this->DebugLevel>=0)
				sjc_debug::debug0(this->DebugDepth,"Potential stored");
		}
	}
	template void PCMBase::formpotential(CQMemManager& mem, double* PDM, EMPerturbation& perb, BasisSet& basisset);
	template void PCMBase::formpotential(CQMemManager& mem, dcomplex* PDM, EMPerturbation& perb, BasisSet& basisset);
	template<class MatsT>
	void PCMBase::addFock(MatsT* fock)
	{
		Eigen::Map<Eigen::Matrix<MatsT,-1,-1>> Fock(fock,this->nB,this->nB);
		Eigen::Map<Eigen::MatrixXd> matpcmfock(this->pcmfock,this->nB,this->nB);
		//std::cout << "Matrices created" << std::endl;
		Fock=matpcmfock+Fock;
		//std::cout << "Addition finished" << std::endl;
	}
	template void PCMBase::addFock(dcomplex* fock);
	template void PCMBase::addFock(double* fock);
	double PCMBase::computeEnergy()
	{
		double energy=this->surp.transpose().dot(this->surc)/2;
		return energy;
	}
}

