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
#ifndef __INCLUDED_SINGLESLATER_FOCK_HPP__
#define __INCLUDED_SINGLESLATER_FOCK_HPP__

#include <singleslater.hpp>

#include <util/time.hpp>
#include <cqlinalg/blasext.hpp>

#include <cqlinalg.hpp>
#include <cqlinalg/svd.hpp>
#include <cqlinalg/blasutil.hpp>
#include <physcon.hpp>
#include <util/matout.hpp>
#include <util/threads.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <cnpy.h>
#include <pcm.hpp>

//#define _DEBUGORTHO

namespace ChronusQ {

  /**
   *  \brief Forms the Fock matrix for a single slater determinant using
   *  the 1PDM.
   *
   *  \param [in] increment Whether or not the Fock matrix is being 
   *  incremented using a previous density
   *
   *  Populates / overwrites fock strorage
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::formFock(
    EMPerturbation &pert, bool increment, double xHFX) {

    size_t NB = this->aoints.basisSet().nBasis;
    size_t NB2 = NB*NB;

    auto GDStart = tick(); // Start time for G[D]

    // Form G[D]
    formGD(pert,increment,xHFX);

    GDDur = tock(GDStart); // G[D] Duraction

    ROOT_ONLY(comm); 


    // Zero out the Fock
    for(auto &F : fockMatrix) std::fill_n(F,NB2,MatsT(0.));

    // Copy over the Core Hamiltonian
    for(auto i = 0; i < coreH.size(); i++)
      SetMat('N',NB,NB,MatsT(1.),coreH[i],NB,fockMatrix[i],NB);

    // Add in the two electron term
    for(auto i = 0ul; i < fockMatrix.size(); i++)
      MatAdd('N','N', NB, NB, MatsT(1.), fockMatrix[i], NB, MatsT(1.), twoeH[i], NB, fockMatrix[i], NB);




    // Add in the electric field contributions
    // FIXME: the magnetic field contribution should go here as well to allow for RT
    // manipulation

    if( pert_has_type(pert,Electric) ) {

      auto dipAmp = pert.getDipoleAmp(Electric);

      // Because IntsT and MatsT need not be the same: MatAdd not possible
      for(auto i = 0;    i < 3;     i++)
      for(auto k = 0ul;  k < NB*NB; k++)
        fockMatrix[SCALAR][k] -= 
          2. * dipAmp[i] * this->aoints.lenElecDipole[i][k];

    }
	  //PCM
	  std::vector<size_t> npy_size={NB,NB};
	  if(this->pcm!=nullptr and this->pcm->use_PCM)
	  {
		  std::cout << std::setprecision(4);
		  std::cout << "PCM iteration begins" << std::endl;
		  //std::cout << "formpotential starts" << std::endl;
		  this->pcm->formpotential(this->aoints.memManager_,this->onePDM[0],pert,this->aoints.basisSet_);
		  //std::cout << "Total potential" << std::endl;
		  //std::cout << this->pcm->surp.transpose() << std::endl;
		  //std::cout << "Nuclear potential" << std::endl;
		  //std::cout << this->pcm->nucp.transpose() << std::endl;
		  //std::cout << "formcharge starts" << std::endl;
		  this->pcm->formcharge();
		  //std::cout << "charge" << std::endl;
		  //std::cout << this->pcm->surc.transpose() << std::endl;
		  //std::cout << "total charge: " << this->pcm->surc.sum() << std::endl;
		  //std::cout << "formFock starts" << std::endl;
		  this->pcm->formFock(this->aoints.memManager_,pert,this->aoints.basisSet_);
		  //std::cout << "addFock starts" << std::endl;
		  this->pcm->addFock(this->fockMatrix[0]);
		  std::cout << "PCM iteration ends" << std::endl;
		  //std::cout << Eigen::Map<Eigen::MatrixXd>(this->pcm->pcmfock,this->pcm->nB,this->pcm->nB) << std::endl;
		  //cnpy::npy_save("pcmFock"+std::to_string(times)+".npy",this->pcm->pcmfock,npy_size,"w");
	  }
	  if (pcm->start_save and (pcm->times%pcm->savestep==0))
	  {
		  for(int i=0;i!=this->fockMatrix.size();++i)
		  {
			  //std::cout << "The " << i << "th component of Fock" << std::endl;
			  //std::cout << Eigen::Map<Eigen::Matrix<MatsT,-1,-1>>(this->fockMatrix[i],this->pcm->nB,this->pcm->nB) << std::endl;
			  cnpy::npy_save("Fock"+std::to_string(i)+"_"+std::to_string(pcm->times)+".npy",this->fockMatrix[i],npy_size,"w");
		  }
		  for(int i=0;i!=this->onePDM.size();++i)
		  {
			  //std::cout << "The " << i << "th component of Density Matrix" << std::endl;
			  //std::cout << Eigen::Map<Eigen::Matrix<MatsT,-1,-1>>(this->onePDM[i],this->pcm->nB,this->pcm->nB) << std::endl;
			  cnpy::npy_save("DM"+std::to_string(i)+"_"+std::to_string(pcm->times)+".npy",this->onePDM[i],npy_size,"w");
		  }
		  std::cout << "Fock and DM Saved" << std::endl;
	  }
	  if (pcm->start_save)
		  pcm->times++;

  
#if 0
    printFock(std::cout);
#endif

  }; // SingleSlater::fockFock


  /**
   *  \brief Forms the Hartree-Fock perturbation tensor
   *
   *  Populates / overwrites GD storage (and JScalar and K storage)
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::formGD(EMPerturbation &pert, bool increment, double xHFX) {

    // Decide list of onePDMs to use
    oper_t_coll &contract1PDM  = increment ? deltaOnePDM : this->onePDM;

    size_t NB = this->aoints.basisSet().nBasis;
    size_t NB2 = NB*NB;

    size_t mpiRank   = MPIRank(comm);
    bool   isNotRoot = mpiRank != 0;


    // Zero out J
    if(not increment) memset(coulombMatrix,0.,NB2*sizeof(MatsT));

    std::vector<TwoBodyContraction<MatsT>> contract =
      { {true, contract1PDM[SCALAR], coulombMatrix, true, COULOMB} };

    // Determine how many (if any) exchange terms to calculate
    if( std::abs(xHFX) > 1e-12 )
    for(auto i = 0; i < exchangeMatrix.size(); i++) {
      contract.push_back(
          {true, contract1PDM[i], exchangeMatrix[i], true, EXCHANGE}
      );

      // Zero out K[i]
      if(not increment) memset(exchangeMatrix[i],0,NB2*sizeof(MatsT));
    }

    this->aoints.twoBodyContract(comm, contract, pert);

    ROOT_ONLY(comm); // Return if not root (J/K only valid on root process)


    // Form GD: G[D] = 2.0*J[D] - K[D]
    if( std::abs(xHFX) > 1e-12 ) {
      for(auto i = 0; i < exchangeMatrix.size(); i++)
        MatAdd('N','N', NB, NB, MatsT(0.), twoeH[i], NB, MatsT(-xHFX), exchangeMatrix[i], NB, twoeH[i], NB);
    } else {
      for(auto i = 0; i < fockMatrix.size(); i++) memset(twoeH[i],0,NB2*sizeof(MatsT));
    }
    // G[D] += 2*J[D]
    MatAdd('N','N', NB, NB, MatsT(1.), twoeH[SCALAR], NB, MatsT(2.), coulombMatrix, NB, twoeH[SCALAR], NB);
      
#if 0
  //printJ(std::cout);
    printK(std::cout);
  //printGD(std::cout);
#endif

  }; // SingleSlater::formGD


  /**
   *  \brief Compute the Core Hamiltonian.
   *
   *  \param [in] typ Which Hamiltonian to build
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::formCoreH(EMPerturbation& emPert, CORE_HAMILTONIAN_TYPE typ) {

    ROOT_ONLY(comm);

    if( coreH.size() != 0 )
      CErr("Recomputing the CoreH is not well-defined behaviour",std::cout);

    size_t NB = this->aoints.basisSet().nBasis;
    size_t nSQ  = NB*NB;

    coreH.emplace_back(memManager.malloc<MatsT>(nSQ));
    if(not iCS and nC == 1 and this->aoints.basisSet().basisType == COMPLEX_GIAO) 
      coreH.emplace_back(memManager.malloc<MatsT>(nSQ));
    else if(nC == 2) {
      coreH.emplace_back(memManager.malloc<MatsT>(nSQ));
      coreH.emplace_back(memManager.malloc<MatsT>(nSQ));
      coreH.emplace_back(memManager.malloc<MatsT>(nSQ));
    }

    this->aoints.computeAOOneE(emPert,oneETerms); // compute the necessary 1e ints

    if(typ == NON_RELATIVISTIC) computeNRCH(emPert,coreH);
    if(std::is_same<MatsT,dcomplex>::value && typ == RELATIVISTIC_X2C_1E) 
      computeX2CCH(emPert,coreH);





    // Compute Orthonormalization trasformations
    computeOrtho();


    // Save the Core Hamiltonian
    if( savFile.exists() ) {

      size_t NB = this->aoints.basisSet().nBasis;
      const std::array<std::string,4> spinLabel =
        { "SCALAR", "MZ", "MY", "MX" };

      for(auto i = 0; i < coreH.size(); i++)
        savFile.safeWriteData("INTS/CORE_HAMILTONIAN_" +
          spinLabel[i], coreH[i], {NB,NB});

    }


  }; // SingleSlater<MatsT,IntsT>::computeCoreH


  /**
   *  \brief Allocate, compute and store the orthonormalization matricies 
   *  over the CGTO basis.
   *
   *  Computes either the Lowdin or Cholesky transformation matricies based
   *  on orthoType
   */ 
  template <typename MatsT, typename IntsT> 
  void SingleSlater<MatsT,IntsT>::computeOrtho() {

    size_t NB = this->aoints.basisSet().nBasis;
    size_t nSQ  = NB*NB;

    // Allocate orthogonalization matricies
    ortho1 = memManager.malloc<MatsT>(nSQ);
    ortho2 = memManager.malloc<MatsT>(nSQ);

    std::fill_n(ortho1,nSQ,0.);
    std::fill_n(ortho2,nSQ,0.);

    // Allocate scratch
    MatsT* SCR1 = memManager.malloc<MatsT>(nSQ);

    // Copy the overlap over to scratch space
    std::copy_n(this->aoints.overlap,nSQ,SCR1);

    if(orthoType == LOWDIN) {

      // Allocate more scratch
      MatsT* sE   = memManager.malloc<MatsT>(NB);
      MatsT* SCR2 = memManager.malloc<MatsT>(nSQ);

      
      // Diagonalize the overlap in scratch S = V * s * V**T
      HermetianEigen('V','U',NB,SCR1,NB,sE,memManager);


      if( std::abs( sE[0] ) < 1e-10 )
        CErr("Contracted Basis Set is Linearly Dependent!");

      // Compute X = V * s^{-1/2} 
      for(auto j = 0; j < NB; j++)
      for(auto i = 0; i < NB; i++)
        SCR2[i + j*NB] = 
          SCR1[i + j*NB] / std::sqrt(sE[j]);

      // Compute O1 = X * V**T
      Gemm('N','C',NB,NB,NB,
        static_cast<MatsT>(1.),SCR2,NB,SCR1,NB,
        static_cast<MatsT>(0.),ortho1,NB);


      // Compute X = V * s^{1/2} in place (by multiplying by s)
      for(auto j = 0; j < NB; j++)
      for(auto i = 0; i < NB; i++)
        SCR2[i + j*NB] = 
          SCR2[i + j*NB] * sE[j];

      // Compute O2 = X * V**T
      Gemm('N','C',NB,NB,NB,
        static_cast<MatsT>(1.),SCR2,NB,SCR1,NB,
        static_cast<MatsT>(0.),ortho2,NB);

#ifdef _DEBUGORTHO
      // Debug code to validate the Lowdin orthogonalization

      std::cerr << "Debugging Lowdin Orthogonalization" << std::endl;
      double maxDiff(-10000000);

      // Check that ortho1 and ortho2 are inverses of eachother
      Gemm('N','N',NB,NB,NB,
        static_cast<MatsT>(1.),ortho1,NB,ortho2,NB,
        static_cast<MatsT>(0.),SCR1,NB);
      
      for(auto j = 0; j < NB; j++)
      for(auto i = 0; i < NB; i++) {

        if( i == j ) maxDiff = 
          std::max(maxDiff, std::abs(1. - SCR1[i + j*NB]));
        else maxDiff = 
          std::max(maxDiff,std::abs(SCR1[i + j*NB])); 

      }

      std::cerr << "  Ortho1 * Ortho2 = I: " << maxDiff << std::endl;

      // Check that ortho2 * ortho2 is the overlap
      Gemm('N','N',NB,NB,NB,
        static_cast<MatsT>(1.),ortho2,NB,ortho2,NB,
        static_cast<MatsT>(0.),SCR1,NB);
      
      maxDiff = -100000;

      for(auto j = 0; j < NB; j++)
      for(auto i = 0; i < NB; i++) {

          maxDiff = std::max(maxDiff,
          std::abs(SCR1[i + j*NB] - 
            this->aoints.overlap[i + j*NB])); 

      }

      std::cerr << "  Ortho2 * Ortho2 = S: " << maxDiff << std::endl;

      // Check that ortho1 * ortho1 is the inverse of the overlap
      Gemm('N','N',NB,NB,NB,
        static_cast<MatsT>(1.),ortho1,NB,ortho1,NB,
        static_cast<MatsT>(0.),SCR1,NB);
      Gemm('N','N',NB,NB,NB,
        static_cast<MatsT>(1.),SCR1,NB,reinterpret_cast<MatsT*>(this->aoints.overlap),NB,
        static_cast<MatsT>(0.),SCR2,
        NB);
      
      maxDiff = -10000;
      for(auto j = 0; j < NB; j++)
      for(auto i = 0; i < NB; i++) {

        if( i == j ) maxDiff = 
          std::max(maxDiff, std::abs(1. - SCR2[i + j*NB]));
        else maxDiff = 
          std::max(maxDiff,std::abs(SCR2[i + j*NB])); 

      }

      std::cerr << "  Ortho1 * Ortho1 * S = I: " << maxDiff << std::endl;

#endif

      // Free Scratch Space
      memManager.free(sE,SCR2);

    } else if(orthoType == CHOLESKY) {

      std::cout << 
      "*** WARNING: Cholesky orthogonalization has not yet been confirmed ***" 
      << std::endl;

      // Compute the Cholesky factorization of the overlap S = L * L**T
      Cholesky('L',NB,SCR1,NB);

      // Copy the lower triangle to ortho2 (O2 = L)
      for(auto j = 0; j < NB; j++)
      for(auto i = j; i < NB; i++)
        ortho2[i + j*NB] = SCR1[i + j*NB];

      // Compute the inverse of the overlap using the Cholesky factors
      CholeskyInv('L',NB,SCR1,NB);

      // O1 = O2**T * S^{-1}
      Gemm('T','N',NB,NB,NB,
        static_cast<MatsT>(1.),ortho2,NB,SCR1,NB,
        static_cast<MatsT>(0.),ortho1,NB);

      // Remove upper triangle junk from O1
      for(auto j = 0; j < NB; j++)
      for(auto i = 0; i < j               ; i++)
        ortho1[i + j*NB] = 0.;

#ifdef _DEBUGORTHO
      // Debug code to validate the Lowdin orthogonalization

      std::cerr << "Debugging Cholesky Orthogonalization" << std::endl;

      // Debug code to validate the Cholesky orthogonalization
      MatsT* SCR2 = memManager.malloc<MatsT>(nSQ);
        
      double maxDiff = -1000;
      Gemm('T','N',NB,NB,NB,
        static_cast<MatsT>(1.),ortho1,NB,reinterpret_cast<MatsT*>(this->aoints.overlap),NB,
        static_cast<MatsT>(0.),SCR1,
        NB);
      Gemm('N','N',NB,NB,NB,
        static_cast<MatsT>(1.),SCR1,NB,ortho1,NB,
        static_cast<MatsT>(0.),SCR2,
        NB);

      for(auto j = 0; j < NB; j++)
      for(auto i = 0; i < NB; i++) {

        if( i == j ) maxDiff = 
          std::max(maxDiff, std::abs(1. - SCR2[i + j*NB]));
        else maxDiff = 
          std::max(maxDiff,std::abs(SCR2[i + j*NB])); 

      }

      std::cerr << "Ortho1**T * S ** Ortho1 = I: " << maxDiff << std::endl;

      memManager.free(SCR2); // Free SCR2
#endif
        

    }

    memManager.free(SCR1); // Free SCR1

  }; // computeOrtho

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::MOFOCK() {

    const size_t NB   = this->aoints.basisSet().nBasis;
    const size_t NB2  = NB * NB;
    const size_t NBC  = nC * NB;
    const size_t NBC2 = NBC * NBC;


    if( fockMO.empty() ) {
      fockMO.emplace_back(memManager.template malloc<MatsT>(NBC2));
      if( nC == 1 and not iCS )
        fockMO.emplace_back(memManager.template malloc<MatsT>(NBC2));
    }

    if( MPIRank(comm) == 0 ) {
      MatsT* SCR1 = memManager.template malloc<MatsT>(NBC2);

      MatsT* FSTORAGE = fockMatrix[0];

      if( nC == 1 and not iCS ) {

        MatAdd('N','N',NB,NB,MatsT(1.),fockMatrix[0],NB,MatsT(1.) ,fockMatrix[1],NB,fockMO[0],NB);
        MatAdd('N','N',NB,NB,MatsT(1.),fockMatrix[0],NB,MatsT(-1.),fockMatrix[1],NB,fockMO[1],NB);

        FSTORAGE = fockMO[0];

      } else if( nC == 2 ) {

        SpinGather(NB,fockMO[0],NBC,fockMatrix[SCALAR],NB,fockMatrix[MZ],NB,
          fockMatrix[MY],NB,fockMatrix[MX],NB);

        FSTORAGE = fockMO[0];

      }

      
      MatsT tfact = (nC == 2) ? 1. : 0.5;

      Gemm('C','N',NBC,NBC,NBC,MatsT(1.),this->mo1,NBC,FSTORAGE,NBC,MatsT(0.),SCR1,NBC);
      Gemm('N','N',NBC,NBC,NBC,tfact,SCR1,NBC,this->mo1,NBC,MatsT(0.),
          fockMO[0],NBC);

      if( nC == 1 and not iCS ) {
        Gemm('C','N',NB,NB,NB,MatsT(1.),this->mo2,NB,fockMO[1],NB,MatsT(0.),SCR1,NB);
        Gemm('N','N',NB,NB,NB,tfact,SCR1,NB,this->mo2,NB,MatsT(0.),fockMO[1],NB);
      }

      memManager.free(SCR1);
    }

  //prettyPrintSmart(std::cout,"MOF",fockMO[0],NBC,NBC,NBC);

#ifdef CQ_ENABLE_MPI

    if(MPISize(comm) > 1) {
      std::cerr  << "  *** Scattering MO-FOCK ***\n";
      for(int k = 0; k < fockMO.size(); k++)
        MPIBCast(fockMO[k],NBC*NBC,0,comm);
    }

#endif
  };

}; // namespace ChronusQ

#endif
