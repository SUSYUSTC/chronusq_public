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
#ifndef __INCLUDED_SINGLESLATER_SCF_HPP__
#define __INCLUDED_SINGLESLATER_SCF_HPP__

#include <singleslater.hpp>
#include <util/matout.hpp>
#include <util/math.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <vector>

// SCF definitions for SingleSlaterBase
#include <singleslater/base/scf.hpp> 

#include <ctime>
#include <sjc_debug.hpp>
#include <cnpy.h>

namespace ChronusQ {

  /**
   *  \brief Saves the current state of wave function
   *
   *  Saves a copy of the current AO 1PDM and orthonormal Fock
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::saveCurrentState() {

    ROOT_ONLY(comm); 


    // Checkpoint if file exists
    if( savFile.exists() ) {

      size_t NB = this->aoints.basisSet().nBasis;
      size_t NBC = this->nC * NB;

      auto t_type_indx = std::type_index(typeid(MatsT));
      size_t t_hash = t_type_indx.hash_code();

      // Save Field type
      savFile.safeWriteData("SCF/FIELD_TYPE",&t_hash,{1});


      const std::array<std::string,4> spinLabel =
        { "SCALAR", "MZ", "MY", "MX" };

      // Save Matricies
      for(auto i = 0; i < this->fockMatrix.size(); i++) {

        savFile.safeWriteData("SCF/1PDM_" + spinLabel[i],
          this->onePDM[i],{NB,NB});

        savFile.safeWriteData("SCF/FOCK_" + spinLabel[i],
          this->fockMatrix[i],{NB,NB});

        savFile.safeWriteData("SCF/1PDM_ORTHO_" + spinLabel[i],
          this->onePDMOrtho[i],{NB,NB});

        savFile.safeWriteData("SCF/FOCK_ORTHO_" + spinLabel[i],
          this->fockMatrixOrtho[i],{NB,NB});

      }

	  if(this->DebugLevel>=2)
		  sjc_debug::debug0(this->DebugDepth,"start saving ortho");
	  try
	  {
		  savFile.safeWriteData("SCF/ORTHO1", this->ortho1,{NB,NB});
	  }
	  catch(...)
	  {
		  std::cout << "No Ortho1" << std::endl;
	  }

	  try
	  {
		  savFile.safeWriteData("SCF/ORTHO2", this->ortho2,{NB,NB});
	  }
	  catch(...)
	  {
		  std::cout << "No Ortho2" << std::endl;
	  }


	  // Save MOs
	  savFile.safeWriteData("SCF/MO1", this->mo1, {NBC,NBC});
	  if( this->nC == 1 and not this->iCS )
		  savFile.safeWriteData("SCF/MO2", this->mo2, {NBC,NBC});

      // Save Energies
      savFile.safeWriteData("SCF/TOTAL_ENERGY",&this->totalEnergy,
        {1});
      savFile.safeWriteData("SCF/ONE_BODY_ENERGY",&this->OBEnergy,
        {1});
      savFile.safeWriteData("SCF/MANY_BODY_ENERGY",&this->MBEnergy,
        {1});

      // Save Multipoles
      savFile.safeWriteData("SCF/LEN_ELECTRIC_DIPOLE",&this->elecDipole[0],
        {3});
      savFile.safeWriteData("SCF/LEN_ELECTRIC_QUADRUPOLE",
        &this->elecQuadrupole[0][0], {3,3});
      savFile.safeWriteData("SCF/LEN_ELECTRIC_OCTUPOLE",
        &this->elecOctupole[0][0][0], {3,3,3});

      // Save Spin
      savFile.safeWriteData("SCF/S_EXPECT",&this->SExpect[0],{3});
      savFile.safeWriteData("SCF/S_SQUARED",&this->SSq,{1});
      

    // If file doesnt exist, checkpoint important bits in core
    } else {

      size_t OSize = memManager.template getSize(fockMatrix[SCALAR]);

      // Copy over current AO density matrix
      for(auto i = 0; i < this->onePDM.size(); i++)
        std::copy_n(this->onePDM[i],OSize,curOnePDM[i]);

      // Copy the previous orthonormal Fock matrix for damping. It's the 
      // previous Fock since saveCurrentState is called at the beginning 
      // of the SCF loop. 
      if ( scfControls.doExtrap and scfControls.doDamp) {
          
        // Avoid saving the guess Fock for extrapolation
        if (scfConv.nSCFIter > 0) {
          for(auto i = 0; i < this->fockMatrixOrtho.size(); i++)
            std::copy_n(this->fockMatrixOrtho[i],OSize,prevFock[i]);
        }
      }

    }

  }; // SingleSlater<MatsT>::saveCurrentState()


  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT, IntsT>::ConventionalSCF(bool modF) {

    // Transform AO fock into the orthonormal basis (on root MPI process)
    ao2orthoFock();

    // Modify fock matrix if requested (on root MPI process)
    if( modF ) modifyFock();

    // Diagonalize the orthonormal fock Matrix (on root MPI process)
    diagOrthoFock();

    // Form the orthonormal density (in the AO storage, 
    // updates all MPI processes)
    formDensity();

    // Copy the AO storage to orthonormal storage and back transform
    // the density into the AO basis. This is because ortho2aoDen
    // requires the onePDMOrtho storage is populated.
    //
    // *** Replicated on all MPI processes ***
    for(auto i = 0; i < this->onePDM.size(); i++)
      std::copy_n(this->onePDM[i],
        memManager.template getSize(onePDMOrtho[i]),
        onePDMOrtho[i]);

    // Transform the orthonormal density to the AO basis 
    // (on root MPI process)
    ortho2aoDen();

    ortho2aoMOs();

  }

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::NewtonRaphsonSCF() {

    MatsT* C = getNRCoeffs();

    // MO(:,i) = MO(:,i) + \sum_a C(a,i) MO(:,a)
    const size_t NB   = this->aoints.basisSet().nBasis;
    const size_t NB2  = NB * NB;
    const size_t NBC  = nC * NB;
    const size_t NBC2 = NBC * NBC;

    const size_t NO    = (nC == 2) ? this->nO : this->nOA;
    const size_t NV    = (nC == 2) ? this->nV : this->nVA;
    const size_t nOAVA = this->nOA * this->nVA;
    const size_t nOBVB = this->nOB * this->nVB;

    for(auto i = 0ul, ai = 0ul; i < NO;  i++)
    for(auto a = NO           ; a < NBC; a++, ai++) 
      AXPY(NBC,-C[ai],this->mo1 + a*NBC,1,this->mo1 + i*NBC,1);

    if( nC == 1 and not iCS )
      for(auto i = 0ul, ai = 0ul; i < this->nOB;  i++)
      for(auto a = this->nOB   ; a < NB; a++, ai++) 
        AXPY(NB,-C[ai + nOAVA],this->mo2 + a*NB,1,this->mo2 + i*NB,1);

    this->memManager.free(C);


    orthoAOMO();
    formDensity();

  }



  /**
   *  \brief Obtain a new set of orbitals given a Fock matrix.
   *
   *  Currently implements the fixed-point SCF procedure.
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::getNewOrbitals(EMPerturbation &pert, 
      bool frmFock) {

    bool increment = scfControls.doIncFock and 
                     scfConv.nSCFIter % scfControls.nIncFock != 0 and
                     scfControls.guess != RANDOM;

	if(this->DebugLevel>=2)
		sjc_debug::debugP(this->DebugDepth,"getNewOrbitals","formFock");
    // Form the Fock matrix D(k) -> F(k)
    if( frmFock ) {
      formFock(pert,increment);
      //if( MPIRank(comm) == 0 )
      //  printFockTimings(std::cout);
    }
	if(this->DebugLevel>=2)
		sjc_debug::debugN(this->DebugDepth,"getNewOrbitals","formFock");

    if( scfControls.scfAlg == _NEWTON_RAPHSON_SCF and scfConv.nSCFIter > 0 )
      scfControls.scfStep = _NEWTON_RAPHSON_STEP;

    if( scfControls.scfStep == _CONVENTIONAL_SCF_STEP )
      ConventionalSCF(scfControls.doExtrap and frmFock);
    else
      NewtonRaphsonSCF();


#ifdef CQ_ENABLE_MPI
    // Broadcast the AO 1PDM to all MPI processes
    if( MPISize(comm) > 1 ) {
      std::cerr  << "  *** Scattering the 1PDM ***\n";
      for(auto k = 0; k < this->onePDM.size(); k++)
        MPIBCast(this->onePDM[k],
          memManager.template getSize(this->onePDM[k]),0,comm);
    }
#endif

  }; // SingleSlater<T>::getNewOrbitals



  template <typename MatsT, typename IntsT>
  bool SingleSlater<MatsT, IntsT>::checkStability() {

    double W;
    MatsT* J;
    std::tie(W,J) = this->getStab();
    std::cout << "  * LOWEST STABILITY EIGENVALUE " << 
      std::scientific << W << "\n";

    if( W < 0. and std::abs(W) > 1e-08 )
      std::cout << "  * LOWEST STABILITY EIGENVALUE NEGATIVE: " 
        << "PERFORMING THOULESS ROTATION\n";
    else { 

      std::cout << "  * LOWEST STABILITY EIGENVALUE POSITIVE: " 
        << "WAVE FUNCTION 2nd ORDER STABLE\n";

      this->memManager.free(J); return true; 

    }

    const size_t NB   = this->aoints.basisSet().nBasis;
    const size_t NB2  = NB * NB;
    const size_t NBC  = nC * NB;
    const size_t NBC2 = NBC * NBC;

    const size_t NO    = (nC == 2) ? this->nO : this->nOA;
    const size_t NV    = (nC == 2) ? this->nV : this->nVA;
    const size_t nOAVA = this->nOA * this->nVA;
    const size_t nOBVB = this->nOB * this->nVB;


    MatsT* ROT    = this->memManager.template malloc<MatsT>(NBC2);
    MatsT* EXPROT = this->memManager.template malloc<MatsT>(NBC2);
    std::fill_n(ROT,NBC2,0.);

    for(auto i = 0ul, ai = 0ul; i < NO;  i++)
    for(auto a = NO;            a < NBC; a++, ai++) {

      ROT[a + i*NBC] =  J[ai];
      ROT[i + a*NBC] = -SmartConj(J[ai]);

    }      

    // FIXME: need to generalize MatExp to take non-hermetian and real
    // matricies
    //MatExp('D',NBC,T(-1.),ROT,NBC,EXPROT,NBC,this->memManager);

    // Taylor
    MatsT s = 1.;
    std::copy_n(ROT,NBC2,EXPROT); // n = 1
    Scale(NBC2,-s,EXPROT,1);
    for(auto j = 0; j < NBC; j++) EXPROT[j*(NBC+1)] += 1.; // n = 0

    MatsT* SCR  = this->memManager.template malloc<MatsT>(NBC2);
    MatsT* SCR2 = this->memManager.template malloc<MatsT>(NBC2);
    std::copy_n(ROT,NBC2,SCR);

    size_t tayMax = 30; 
    for(auto n = 2; n <= tayMax; n++) {

      MatsT* M = nullptr;
      if( n % 2 ) {
        Gemm('N','N',NBC,NBC,NBC,MatsT(1.),ROT,NBC,SCR2,NBC,MatsT(0.),SCR,NBC);
        M = SCR;
      } else {
        Gemm('N','N',NBC,NBC,NBC,MatsT(1.),ROT,NBC,SCR,NBC,MatsT(0.),SCR2,NBC);
        M = SCR2;
      }

      MatsT fact = std::pow(-s,n)/factorial(n);
      MatAdd('N','N',NBC,NBC,MatsT(1.),EXPROT,NBC,fact,M,NBC,
        EXPROT,NBC);

    }

    /*
    Gemm('C','N',NBC,NBC,NBC,T(1.),EXPROT,NBC,EXPROT,NBC,T(0.),SCR,NBC);
   // prettyPrintSmart(std::cerr,"ROT",ROT,NBC,NBC,NBC);
   // prettyPrintSmart(std::cerr,"EXPROT",EXPROT,NBC,NBC,NBC);
    prettyPrintSmart(std::cout,"SCR",SCR,NBC,NBC,NBC);
   // CErr();
   */


    // MO1 = MO1 * EXPROT
    Gemm('N','N',NBC,NBC,NBC,MatsT(1.),this->mo1,NBC,EXPROT,NBC,MatsT(0.),
      ROT,NBC);
    std::copy_n(ROT,NBC2,this->mo1);




    orthoAOMO();
    this->formDensity();
    this->formDelta();
    this->memManager.free(J,ROT,EXPROT);

    return false;

  }





  /**
   *  \brief Evaluate SCF convergence based on various criteria.
   *
   *  Checks the norm of [F,D], if converged -> SCF converged.
   *
   *  Checks change in energy and density between SCF iterations,
   *    if *both* converged -> SCF converged.
   */ 
  template <typename MatsT, typename IntsT>
  bool SingleSlater<MatsT, IntsT>::evalConver(EMPerturbation &pert) {

    bool isConverged;

    formDelta(); // Get change in density on all MPI processes

    // Compute all SCF convergence information on root process
    if( MPIRank(comm) == 0 ) {
      
      // Check energy convergence
        
      // Save copy of old Energy
      double oldEnergy = this->totalEnergy;

      // Compute new energy (with new Density)
      this->computeProperties(pert);
	if(this->pcm!=nullptr and this->pcm->use_PCM)
		this->totalEnergy+=this->pcm->computeEnergy();
      scfConv.deltaEnergy = this->totalEnergy - oldEnergy;

      bool energyConv = std::abs(scfConv.deltaEnergy) < 
                        scfControls.eneConvTol;

      /*
      bool energySuperConv = std::abs(scfConv.deltaEnergy) < 
                             1e-2 * scfControls.eneConvTol;
                             */
      bool energySuperConv = false;



      // Check density convergence
 
      size_t DSize = memManager. template getSize(fockMatrix[SCALAR]);
      size_t NB    = this->aoints.basisSet().nBasis;
      scfConv.RMSDenScalar = 
        TwoNorm<double>(DSize,deltaOnePDM[SCALAR],1) / NB;
      scfConv.RMSDenMag = 0.;
      for(auto i = 1; i < deltaOnePDM.size(); i++)
        scfConv.RMSDenMag += 
          std::pow(TwoNorm<double>(DSize,deltaOnePDM[i],1),2.);
 
      scfConv.RMSDenMag = std::sqrt(scfConv.RMSDenMag) / NB;
      
 
      bool denConv = scfConv.RMSDenScalar < scfControls.denConvTol;
 
      // Check FP convergence
      bool FDConv(false);
 
      isConverged = FDConv or (energyConv and denConv) or energySuperConv;
 
      // Toggle damping based on energy difference
      if( scfControls.doExtrap ) {
        // TODO: should enable print statements only when print flag is high 
        // enough 
        bool largeEDiff = 
          std::abs(scfConv.deltaEnergy) > scfControls.dampError;
 
        if( scfControls.doDamp and not largeEDiff and 
            scfControls.dampParam > 0.) {
 
          if( printLevel > 0 )
            std::cout << 
              "    *** Damping Disabled - Energy Difference Fell Below " <<
              scfControls.dampError << " ***" << std::endl;
 
          scfControls.dampParam = 0.;
 
        } else if( scfControls.doDamp and largeEDiff and 
                   scfControls.dampParam <= 0.) {
 
          if( printLevel > 0 )
            std::cout << "    *** Damping Enabled Due to "<<
              scfControls.dampError << " Oscillation in Energy ***" << std::endl;
 
          scfControls.dampParam = scfControls.dampStartParam;
 
        }
      }

    }


    // If converged and NR, check for saddle point and poke along
    // that direction ala https://doi.org/10.1063/1.4918561
    if( scfControls.scfAlg == _NEWTON_RAPHSON_SCF and isConverged )
      isConverged = checkStability();


#ifdef CQ_ENABLE_MPI
    // Broadcast whether or not we're converged to ensure that all
    // MPI processes exit the SCF simultaneously
    if( MPISize(comm) > 1 ) MPIBCast(isConverged,0,comm);
#endif

    return isConverged;

  }; // SingleSlater<MatsT>::evalConver





  /**
   *  \brief Computes the change in the current wave function
   *
   *  Saves onePDM - curOnePDM in deltaOnePDM
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::formDelta() {

    size_t NB = this->aoints.basisSet().nBasis;

    // Compute difference on root MPI process
    if( MPIRank(comm) == 0 ) {
      if( not savFile.exists() )
        for(auto i = 0; i < this->onePDM.size(); i++)
          MatAdd('N','N',NB,NB,MatsT(1.),this->onePDM[i],NB,MatsT(-1.),
            curOnePDM[i],NB,deltaOnePDM[i],NB);
      else {

        MatsT* DENSCR = this->memManager.template malloc<MatsT>(NB*NB);
        const std::array<std::string,4> spinLabel =
          { "SCALAR", "MZ", "MY", "MX" };

        for(auto i = 0; i < this->onePDM.size(); i++) {

          savFile.readData("/SCF/1PDM_" + spinLabel[i],DENSCR);

          MatAdd('N','N',NB,NB,MatsT(1.),this->onePDM[i],NB,MatsT(-1.),
            DENSCR,NB,deltaOnePDM[i],NB);
        }

        this->memManager.free(DENSCR);

      }
    }

#ifdef CQ_ENABLE_MPI
    // Broadcast the change in the 1PDM
    if( MPISize(comm) > 1 ) 
      for(auto &X : deltaOnePDM) {
        MPIBCast(X,NB*NB,0,comm);
      }
#endif

  }; // SingleSlater<MatsT>:formDelta

  



  /**
   *  \brief Diagonalize the orthonormal fock matrix
   *
   *  General purpose routine which diagonalizes the orthonormal
   *  fock matrix and stores a set of orthonormal MO coefficients
   *  (in WaveFunction::mo1 and possibly WaveFunction::mo2) and
   *  orbital energies. General for both 1 and 2 spin components
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::diagOrthoFock() {

    ROOT_ONLY(comm); 
    size_t NB = this->aoints.basisSet().nBasis * nC;
    size_t NB2 = NB*NB;

    // Copy over the fockMatrixOrtho into MO storage
    if(nC == 1 and iCS) 
      std::transform(fockMatrixOrtho[SCALAR],fockMatrixOrtho[SCALAR] + NB2,this->mo1,
        [](MatsT a){ return a / 2.; }
      );
    else if(nC == 1)
      for(auto j = 0; j < NB2; j++) {
        this->mo1[j] = 0.5 * (fockMatrixOrtho[SCALAR][j] + fockMatrixOrtho[MZ][j]); 
        this->mo2[j] = 0.5 * (fockMatrixOrtho[SCALAR][j] - fockMatrixOrtho[MZ][j]); 
      }
    else { 

      SpinGather(NB/2,this->mo1,NB,fockMatrixOrtho[SCALAR],NB/2,fockMatrixOrtho[MZ],
        NB/2,fockMatrixOrtho[MY],NB/2,fockMatrixOrtho[MX],NB/2);

    }

    // Diagonalize the Fock Matrix
    int INFO = HermetianEigen('V', 'L', NB, this->mo1, NB, this->eps1, 
      memManager );
    if( INFO != 0 ) CErr("HermetianEigen failed in Fock1",std::cout);

    if(nC == 1 and not iCS) {
      INFO = HermetianEigen('V', 'L', NB, this->mo2, NB, this->eps2, 
        memManager );
      if( INFO != 0 ) CErr("HermetianEigen failed in Fock2",std::cout);
    }

#if 0
    printMO(std::cout);
#endif

  }; // SingleSlater<MatsT>::diagOrthoFock


  /**
   *  \brief Transforms all of the spin components of the AO fock
   *  matrix to the orthonormal basis.
   *
   *  Populates / overwrites fockMatrixOrtho storage
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::ao2orthoFock() {

    ROOT_ONLY(comm); 
    for(auto i = 0; i < fockMatrix.size(); i++)
      Ortho1TransT(fockMatrix[i],fockMatrixOrtho[i]);

  }; // SingleSlater<MatsT>::ao2orthoFock



  /**
   *  \brief Transforms all of the spin compoenents of the orthonormal
   *  1PDM to the AO basis.
   *
   *  Populates / overwrites onePDM storage
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::ortho2aoDen() {

    ROOT_ONLY(comm);

    for(auto i = 0; i < onePDMOrtho.size(); i++)
      Ortho1Trans(onePDMOrtho[i],this->onePDM[i]);

#if 0
    print1PDMOrtho(std::cout);
#endif

  }; // SingleSlater<MatsT>::ao2orthoFock

  
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::ortho2aoMOs() {

    size_t NB = this->aoints.basisSet().nBasis;

    // Transform MOs on MPI root as slave processes do not have
    // updated MO coefficients
    if( MPIRank(comm) == 0 ) {
      MatsT* SCR = this->memManager.template malloc<MatsT>(this->nC*NB*NB);

      // Transform the (top half) of MO1
      Gemm('N','N',NB,this->nC*NB,NB,MatsT(1.),this->ortho1,NB,
        this->mo1,this->nC*NB,MatsT(0.),SCR,NB);
      SetMat('N',NB,this->nC*NB,MatsT(1.),SCR,NB,this->mo1,this->nC*NB);
      
      if( this->nC == 2 ) {

        // Transform the bottom half of MO1
        Gemm('N','N',NB,this->nC*NB,NB,MatsT(1.),this->ortho1,NB,
          this->mo1 + NB,this->nC*NB,MatsT(0.),SCR,NB);
        SetMat('N',NB,this->nC*NB,MatsT(1.),SCR,NB,this->mo1 + NB,this->nC*NB);

      } else if( not this->iCS ) {

        // Transform MO2
        Gemm('N','N',NB,NB,NB,MatsT(1.),this->ortho1,NB,this->mo2,NB,MatsT(0.),
          SCR,NB);
        SetMat('N',NB,NB,MatsT(1.),SCR,NB,this->mo2,NB);

      }

      this->memManager.free(SCR);

    }

#ifdef CQ_ENABLE_MPI


    // Broadcast the updated MOs to all MPI processes
    if( MPISize(comm) > 1 ) {

      std::cerr  << "  *** Scattering the AO-MOs ***\n";
      MPIBCast(this->mo1,nC*nC*NB*NB,0,comm);
      if( nC == 1 and not iCS )
        MPIBCast(this->mo2,nC*nC*NB*NB,0,comm);

      std::cerr  << "  *** Scattering EPS ***\n";
      MPIBCast(this->eps1,nC*NB,0,comm);
      if( nC == 1 and not iCS )
        MPIBCast(this->eps2,nC*NB,0,comm);

      std::cerr  << "  *** Scattering FOCK ***\n";
      for(int k = 0; k < fockMatrix.size(); k++)
        MPIBCast(fockMatrix[k],NB*NB,0,comm);

    }

#endif

    MOFOCK(); // Form the MO fock matrix

#if 0
    // Check proper orthonormalized wrt overlap

    T* SCR2 = this->memManager.template malloc<T>(this->nC*this->nC*NB*NB);
    T* SCR3 = this->memManager.template malloc<T>(this->nC*this->nC*NB*NB);


    // MO1 inner product
    Gemm('N','N',NB,this->nC*NB,NB,T(1.),this->aoints.overlap,NB,
      this->mo1,this->nC*NB,T(0.),SCR2,this->nC*NB);
    if(this->nC == 2)
      Gemm('N','N',NB,this->nC*NB,NB,T(1.),this->aoints.overlap,NB,
        this->mo1+NB,this->nC*NB,T(0.),SCR2+NB,this->nC*NB);
   
    Gemm('C','N',this->nC*NB,this->nC*NB,this->nC*NB,T(1.),this->mo1,
      this->nC*NB,SCR2,this->nC*NB,T(0.),SCR3,this->nC*NB);

    for(auto i = 0; i < this->nC*NB; i++)
      SCR3[i*(this->nC*NB + 1)] -= 1.;


    std::cerr << "Error in orthonormazation of MO1 = " 
      << MatNorm<double>('F',this->nC*NB,this->nC*NB,SCR3,this->nC*NB)
      << std::endl;
             


    if(this->nC == 1 and not this->iCS) {
      Gemm('N','N',NB,NB,NB,T(1.),this->aoints.overlap,NB,this->mo2,NB,T(0.),SCR2,NB);
      Gemm('C','N',NB,NB,NB,T(1.),this->mo2,NB,SCR2,NB,T(0.),SCR3,NB);

      for(auto i = 0; i < this->nC*NB; i++)
        SCR3[i*(this->nC*NB + 1)] -= 1.;

      std::cerr << "Error in orthonormazation of MO2 = " 
        << MatNorm<double>('F',NB,NB,SCR3,NB) << std::endl;
    }

    this->memManager.free(SCR2,SCR3);

#endif

  }; // SingleSlater<MatsT>::ortho2aoMOs


  /**
   *  \brief Initializes the environment for the SCF caluclation.
   *
   *  Allocate memory for extrapolation and compute the energy
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::SCFInit() {

    // Allocate additional storage if doing some type of 
    // extrapolation during the SCF procedure
    if ( scfControls.doExtrap ) allocExtrapStorage();

  }; // SingleSlater<MatsT>::SCFInit




  /**
   *  \brief Finalizes the environment for the SCF caluclation.
   *
   *  Deallocate the memory allocated for extrapolation.
   */ 
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::SCFFin() {

    // Deallocate extrapolation storage
    if ( scfControls.doExtrap ) deallocExtrapStorage();

  }; // SingleSlater<MatsT>::SCFFin




  /**
   *  \brief Reorthogonalize the MOs wrt overlap
   */
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT, IntsT>::orthoAOMO() {

    const size_t NB   = this->aoints.basisSet().nBasis;
    const size_t NB2  = NB * NB;
    const size_t NBC  = this->nC * NB;
    const size_t NBC2 = NBC * NBC;

    // Reorthogonalize MOs wrt S
    MatsT* SCR  = this->memManager.template malloc<MatsT>(NBC*NBC);
    MatsT* SCR2 = this->memManager.template malloc<MatsT>(NBC*NBC);

    // SCR2 = C**H S C
    Gemm('N','N',NB,NBC,NB,MatsT(1.),this->aoints.overlap,NB,this->mo1,NBC,MatsT(0.),
      SCR,NBC);
    if( this->nC == 2 )
      Gemm('N','N',NB,NBC,NB,MatsT(1.),this->aoints.overlap,NB,this->mo1+NB,NBC,
        MatsT(0.),SCR+NB,NBC);

    Gemm('C','N',NBC,NBC,NBC,MatsT(1.),this->mo1,NBC,SCR,NBC,MatsT(0.),SCR2,NBC);

    // SCR2 = L L**H -> L
    int INFO = Cholesky('L',NBC,SCR2,NBC);

    // SCR2 = L^-1
    INFO = TriInv('L','N',NBC,SCR2,NBC);

    // MO1 = MO1 * L^-H
    Trmm('R','L','C','N',NBC,NBC,MatsT(1.),SCR2,NBC,this->mo1,NBC);

    // Reorthogonalize MOB
    if( this->nC == 1 and not this->iCS ) {
      Gemm('N','N',NB,NB,NB,MatsT(1.),this->aoints.overlap,NB,this->mo2,NB,MatsT(0.),
        SCR,NB);
      Gemm('C','N',NB,NB,NB,MatsT(1.),this->mo2,NBC,SCR,NB,MatsT(0.),SCR2,NB);

      // SCR2 = L L**H -> L
      INFO = Cholesky('L',NB,SCR2,NB);

      // SCR2 = L^-1
      INFO = TriInv('L','N',NB,SCR2,NB);

      // MO2 = MO2 * L^-H
      Trmm('R','L','C','N',NB,NB,MatsT(1.),SCR2,NB,this->mo2,NB);

    }

    this->memManager.free(SCR,SCR2);
  }

	template<typename MatsT, typename IntsT>
	void SingleSlater<MatsT,IntsT>::swaporbit(bool is_swap_alpha, bool is_swap_beta, std::pair<int,int> swap_alpha, std::pair<int,int> swap_beta) {
		if(this->nC==1)
		{
			typedef Eigen::Matrix<MatsT,-1,-1> Matrix;
			typedef Eigen::Matrix<MatsT,-1,-1,Eigen::RowMajor> RowMatrix;
			typedef Eigen::Map<Eigen::Matrix<MatsT,-1,-1>> MapMatrix;
			typedef Eigen::Map<Eigen::Matrix<MatsT,-1,-1,Eigen::RowMajor>> MapRowMatrix;
			size_t NB = this->aoints.basisSet().nBasis;
			MapRowMatrix AMO(this->mo1,NB,NB);
			MapRowMatrix BMO(this->mo2,NB,NB);
			if(this->DebugLevel>=2)
			{
				std::cout << "before change" << std::endl;
				std::cout << "AMO" << std::endl;
				std::cout << AMO << std::endl;
				std::cout << "BMO" << std::endl;
				std::cout << BMO << std::endl;
				for(size_t i=0;i<this->onePDM.size();++i)
				{
					std::cout << "onePDM " << i << std::endl;
					std::cout << MapRowMatrix(this->onePDM[i],NB,NB) << std::endl;
				}
				for(size_t i=0;i<this->onePDMOrtho.size();++i)
				{
					std::cout << "onePDMOrtho " << i << std::endl;
					std::cout << MapRowMatrix(this->onePDMOrtho[i],NB,NB) << std::endl;
				}
			}
			Eigen::Matrix<MatsT,1,-1> temprow;
			if(is_swap_alpha)
			{
				temprow=AMO.row(this->nOA-1-swap_alpha.first);
				AMO.row(this->nOA-1-swap_alpha.first)=AMO.row(this->nOA+swap_alpha.second);
				AMO.row(this->nOA+swap_alpha.second)=temprow;
			}
			if(is_swap_beta)
			{
				temprow=BMO.row(this->nOB-1-swap_beta.first);
				BMO.row(this->nOB-1-swap_beta.first)=BMO.row(this->nOB+swap_beta.second);
				BMO.row(this->nOB+swap_beta.second)=temprow;
			}
			this->formDensity();
			for(size_t i=0;i<this->onePDM.size();++i)
				this->Ortho2TransT(this->onePDM[i],this->onePDMOrtho[i]);
			if(this->DebugLevel>=2)
			{
				std::cout << "after change" << std::endl;
				std::cout << "AMO" << std::endl;
				std::cout << AMO << std::endl;
				std::cout << "BMO" << std::endl;
				std::cout << BMO << std::endl;
				for(size_t i=0;i<this->onePDM.size();++i)
				{
					std::cout << "onePDM " << i << std::endl;
					std::cout << MapRowMatrix(this->onePDM[i],NB,NB) << std::endl;
				}
				for(size_t i=0;i<this->onePDMOrtho.size();++i)
				{
					std::cout << "onePDMOrtho " << i << std::endl;
					std::cout << MapRowMatrix(this->onePDMOrtho[i],NB,NB) << std::endl;
				}
			}
		}
	}

	template<typename MatsT, typename IntsT>
	void SingleSlater<MatsT,IntsT>::readpdm(){
		size_t NB = this->aoints.basisSet().nBasis;
		size_t size=this->onePDM.size();
		for(size_t i=0;i<size;++i)
		{
			if(this->DebugLevel>=1)
				sjc_debug::debug0(this->DebugDepth,"Loading "+std::to_string(int(i)));
			std::string path="restart_"+std::to_string(int(i))+".npy";
			cnpy::NpyArray arr=cnpy::npy_load(path);
			MatsT* data=arr.data<MatsT>();
			for(size_t j=0;j<NB*NB;++j)
				this->onePDM[i][j]=data[j];
			this->Ortho2TransT(this->onePDM[i],this->onePDMOrtho[i]);
			//TODO: memery leak here
			if(this->DebugLevel>=1)
				sjc_debug::debug0(this->DebugDepth,"Loading "+std::to_string(i)+" Done");
		}
	}
	template<typename MatsT, typename IntsT>
	void SingleSlater<MatsT,IntsT>::savenpy(std::string word){

		size_t NB = this->aoints.basisSet().nBasis;
		std::vector<size_t> npy_size={NB,NB};
		for(int i=0;i!=this->fockMatrix.size();++i)
			cnpy::npy_save("Fock"+std::to_string(i)+"_"+word+".npy",this->fockMatrix[i],npy_size,"w");
		for(int i=0;i!=this->onePDM.size();++i)
			cnpy::npy_save("DM"+std::to_string(i)+"_"+word+".npy",this->onePDM[i],npy_size,"w");
		if (this->DebugLevel>=0)
			sjc_debug::debug0(this->DebugDepth,"Fock and DM Saved");
	}

}; // namespace ChronusQ

#endif
