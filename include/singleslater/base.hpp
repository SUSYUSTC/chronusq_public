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
#ifndef __INCLUDED_SINGLESLATER_BASE_HPP__
#define __INCLUDED_SINGLESLATER_BASE_HPP__

#include <chronusq_sys.hpp>
#include <wavefunction/base.hpp>
#include <aointegrals.hpp>

#include <fields.hpp>
#include <util/files.hpp>
#include <pcm.hpp>

namespace ChronusQ {

  enum DIIS_ALG {
    CDIIS,      ///< Commutator DIIS
    EDIIS,      ///< Energy DIIS
    CEDIIS,     ///< Commutator & Energy DIIS
    NONE = -1  
  };

  /**
   *  The Single Slater guess types
   */ 
  enum SS_GUESS {
    CORE,
    SAD,
    RANDOM,
    READMO,
    READDEN
  };

  /**
   *  The types of steps for the SCF
   */
  enum SCF_STEP {
    _CONVENTIONAL_SCF_STEP,
    _NEWTON_RAPHSON_STEP
  };

  /**
   *  SCF Algorithms
   */
  enum SCF_ALG {
    _CONVENTIONAL_SCF,
    _NEWTON_RAPHSON_SCF
  };

  /**
   *  \brief A struct to hold the information pertaining to
   *  the control of an SCF procedure.
   *
   *  Holds information like convergence critera, DIIS settings, 
   *  max iterations, etc.
   */ 
  struct SCFControls {

    // Convergence criteria
    double denConvTol = 1e-8;  ///< Density convergence criteria
    double eneConvTol = 1e-10; ///< Energy convergence criteria

    // TODO: need to add logic to set this
    // Extrapolation flag for DIIS and damping
    bool doExtrap = true;     ///< Whether to extrapolate Fock matrix


    // Algorithm and step
    SCF_STEP  scfStep = _CONVENTIONAL_SCF_STEP;
    SCF_ALG   scfAlg  = _CONVENTIONAL_SCF;

    // Guess Settings
    SS_GUESS guess = SAD;

    // DIIS settings 
    DIIS_ALG diisAlg = CDIIS; ///< Type of DIIS extrapolation 
    size_t nKeep     = 10;    ///< Number of matrices to use for DIIS

    // Static Damping settings
    bool   doDamp         = true;           ///< Flag for turning on damping
    double dampStartParam = 0.7;            ///< Starting damping parameter
    double dampParam      = dampStartParam; ///< Current Damp parameter 
    double dampError      = 1e-1; ///< Energy oscillation to turn off damp

    // Incremental Fock build settings
    bool   doIncFock = true; ///< Whether to perform an incremental fock build
    size_t nIncFock  = 20;   ///< Restart incremental fock build after n steps

    // Misc control
    size_t maxSCFIter = 128; ///< Maximum SCF iterations.



    // Printing
    bool printMOCoeffs = false;

  }; // SCFControls struct

  /**
   *  \brief A struct to hold the current status of an SCF procedure
   *
   *  Holds information like current density / energy changes, number of 
   *  iterations, etc.
   */ 
  struct SCFConvergence {

    double deltaEnergy;  ///< Convergence of Energy
    double RMSDenScalar; ///< RMS change in Scalar density
    double RMSDenMag;    ///< RMS change in magnetization (X,Y,Z) density
    double nrmFDC;       ///< 2-Norm of [F,D]

    size_t nSCFIter = 0; ///< Number of SCF Iterations

  }; // SCFConvergence struct


  /**
   *  \brief The SingleSlaterBase class. The abstraction of information
   *  relating to the SingleSlater class which are independent of storage
   *  type.
   *
   *  Specializes WaveFunctionBase interface.
   *
   *  See SingleSlater for further docs.
   */ 
  class SingleSlaterBase : virtual public WaveFunctionBase {

  protected:

    std::string refLongName_;  ///< Long form of the reference name
    std::string refShortName_; ///< Short form of the reference name

  private:
  public:
	size_t DebugLevel;
	size_t DebugDepth=0;
	bool save_status=false;
	size_t save_time=0;
	size_t save_step;

    // Save / Restart File
    SafeFile savFile;
       
    // Print Controls
    size_t printLevel; ///< Print Level

    // Current Timings
    double GDDur;
              
    // Integral variables
    CORE_HAMILTONIAN_TYPE coreType   = NON_RELATIVISTIC;  ///< Core Hamiltonian type
    ORTHO_TYPE            orthoType  = LOWDIN; ///< Orthogonalization scheme

    OneETerms oneETerms; ///< One electron terms to be computed

    // SCF Variables
    SCFControls    scfControls; ///< Controls for the SCF procedure
    SCFConvergence scfConv;     ///< Current status of SCF convergence
	//PCM
	std::shared_ptr<PCMBase> pcm=nullptr;

    // Constructors (all defaulted)
    SingleSlaterBase(const SingleSlaterBase &) = default;
    SingleSlaterBase(SingleSlaterBase &&)      = default;

    SingleSlaterBase() = delete;

    SingleSlaterBase(MPI_Comm c, CQMemManager &mem, size_t _nC, bool iCS) : 
      WaveFunctionBase(c, mem,_nC,iCS), QuantumBase(c, mem,_nC,iCS),
      printLevel((MPIRank(c) == 0) ? 1 : 0) { };

    
    // Setters
    inline void setCoreH(CORE_HAMILTONIAN_TYPE cType) {

      coreType = cType;

      if(coreType == NON_RELATIVISTIC) 
        oneETerms={false,true,false};

      else if(coreType == RELATIVISTIC_X2C_1E) 
        oneETerms={true,true,true};

    }
      


    // Procedural Functions to be defined in all derived classes
      
    // In essence, all derived classes should be able to:
    //   Form a Fock matrix with the ability to increment
    virtual void formFock(EMPerturbation &, bool increment = false, double xHFX = 1.) = 0;

    //   Form an initial Guess (which populates the Fock, Density 
    //   and energy)
    virtual void formGuess() = 0;

    //   Form the core Hamiltonian
    virtual void formCoreH(EMPerturbation&) = 0;

    //   Obtain a new set of orbitals / densities from current
    //   set of densities
    virtual void getNewOrbitals(EMPerturbation &, bool frmFock = true) = 0;

    //   Save the current state of the wave function
    virtual void saveCurrentState() = 0;

    //   Save some metric regarding the change in the wave function
    //   from the currently saved state (i.e. between SCF iterations)
    virtual void formDelta() = 0;

    //   Evaluate SCF convergence. This function should populate the
    //   SingleSlaterBase::scfConv variable and compare it to the 
    //   SingleSlaterBase::scfControls variable to evaluate convergence
    virtual bool evalConver(EMPerturbation &) = 0;

    //   Print SCF header, footer and progress
    void printSCFHeader(std::ostream &out, EMPerturbation &);
    void printSCFProg(std::ostream &out = std::cout,
      bool printDiff = true);

    //   Initialize and finalize the SCF environment
    virtual void SCFInit() = 0;
    virtual void SCFFin()  = 0;

    //   Print various matricies
    virtual void printFock(std::ostream& )     = 0;
    virtual void print1PDMOrtho(std::ostream&) = 0;
    virtual void printGD(std::ostream&)        = 0;
    virtual void printJ(std::ostream&)         = 0;
    virtual void printK(std::ostream&)         = 0;

    virtual void printFockTimings(std::ostream&) = 0;

    // Procedural Functions to be shared among all derived classes
      
    // Perform an SCF procedure (see include/singleslater/scf.hpp for docs)
    void SCF(EMPerturbation &);
	void initpcm(std::shared_ptr<PCMBase>);
	virtual void swaporbit(bool, bool, std::pair<int,int>, std::pair<int,int>)=0;
	virtual	void readpdm()=0;
	virtual void savenpy(std::string)=0;

  }; // class SingleSlaterBase

}; // namespace ChronusQ

#endif
