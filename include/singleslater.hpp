/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2017 Li Research Group (University of Washington)
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
#ifndef __INCLUDED_SINGLESLATER_HPP__
#define __INCLUDED_SINGLESLATER_HPP__

#include <chronusq_sys.hpp>
#include <wavefunction.hpp>

// Debug print triggered by Wavefunction
  
#ifdef _WaveFunctionDebug
  #define _SingleSlaterDebug
#endif

namespace ChronusQ {

  enum DIIS_ALG {
    CDIIS,      ///< Commutator DIIS
    EDIIS,      ///< Energy DIIS
    CEDIIS,     ///< Commutator & Energy DIIS
    NONE = -1  
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

    // DIIS settings 
    DIIS_ALG diisAlg = CDIIS; ///< Type of DIIS extrapolation 
    size_t nKeep     = 10;    ///< Number of matrices to use for DIIS

    // Static Damping settings
    bool   doDamp         = true;           ///< Flag for turning on damping
    double dampStartParam = 0.7;            ///< Starting damping parameter
    double dampParam      = dampStartParam; ///< Damp parameter current iteration
    double dampError      = 1e-1;           ///< Energy oscillation to turn off damp

    // Misc control
    size_t maxSCFIter = 128; ///< Maximum SCF iterations.

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

  template <typename T>
  class SingleSlater : public WaveFunction<T> {

  protected:
    // Useful typedefs
    typedef T*                        oper_t;
    typedef std::vector<oper_t>       oper_t_coll;
    typedef std::vector<oper_t_coll>  oper_t_coll2;

    std::string refLongName_;
    std::string refShortName_;

  private:
  public:
    
    // Operator storage

    // AO Fock Matrix
    oper_t_coll fock;       ///< List of populated AO Fock matricies

    // Orthonormal Fock
    oper_t_coll fockOrtho;   ///< List of populated orthonormal Fock matricies

    // Coulomb (J[D])
    double* JScalar; ///< Scalar Coulomb Matrix

    // Exchange (K[D])
    oper_t_coll K;       ///< List of populated exact (HF) exchange matricies

    // Exact Perturbation Tensor (G[D])
    oper_t_coll GD;      ///< List of populated HF perturbation tensors


    // Orthonormal density
    oper_t_coll onePDMOrtho;   ///< List of populated orthonormal 1PDM matricies


    // SCF Variables
    SCFControls    scfControls; ///< Controls for the SCF procedure
    SCFConvergence scfConv;     ///< Current status of SCF convergence


    // Current / change in state information (for use with SCF)
    oper_t_coll curOnePDM;    ///< List of the current 1PDMs
    oper_t_coll deltaOnePDM;  ///< List of the changes in the 1PDMs

    // Stores the previous Fock matrix to use for damping    
    oper_t_coll prevFock;     ///< AO Fock from the previous SCF iteration

    // Stores the previous matrices to use for DIIS 
    oper_t_coll2 diisFock;    ///< List of AO Fock matrices for DIIS extrap
    oper_t_coll2 diisOnePDM;  ///< List of AO Density matrices for DIIS extrap
    oper_t_coll2 diisError;   ///< List of orthonormal [F,D] for DIIS extrap

    // Constructors
      
    /**
     *  SingleSlater Constructor. Constructs a SingleSlater object
     *
     *  \param [in] aoi  AOIntegrals object (which handels the BasisSet, etc)
     *  \param [in] nC  Number of spin components (1 and 2 are supported)
     */ 
    SingleSlater(AOIntegrals &aoi, size_t nC) : 
      WaveFunction<T>(aoi,nC), JScalar(nullptr) {

      // Allocate SingleSlater Object
      alloc(); 

      // Determine method string
      if(std::is_same<T,double>::value) {
        refLongName_  = "Real ";
        refShortName_ = "R-";
      } else {
        refLongName_  = "Complex ";
        refShortName_ = "C-";
      }

      if(this->nC == 1) {
        if(this->iCS) {
          refLongName_  += "Restricted Hartree-Fock";
          refShortName_ += "RHF";
        } else {
          refLongName_  += "Unrestricted Hartree-Fock";
          refShortName_ += "UHF";
        }
      } else {
        refLongName_  += "Generalized Hartree-Fock";
        refShortName_ += "GHF";
      }

    }

    // See include/singleslater/impl.hpp for documentation 
    // on the following constructors

    // Different type
    template <typename U> SingleSlater(const SingleSlater<U> &, int dummy = 0);
    template <typename U> SingleSlater(SingleSlater<U> &&     , int dummy = 0);

    // Same type
    SingleSlater(const SingleSlater &);
    SingleSlater(SingleSlater &&);     

    /**
     *  Destructor.
     *
     *  Destructs a SingleSlater object
     */ 
    ~SingleSlater() { dealloc(); }



    // Public Member functions

    // Deallocation (see include/singleslater/impl.hpp for docs)
    void alloc();
    void dealloc();


    // Declarations from Quantum 
    // (see include/singleslater/quantum.hpp for docs)
    void formDensity();
    void computeEnergy();

    // Form a fock matrix (see include/singleslater/fock.hpp for docs)
    virtual void formFock(bool increment = false);
    void formGD();

    // Form initial guess orbitals
    // see include/singleslater/guess.hpp for docs)
    void formGuess();
    



    // SCF procedural functions (see include/singleslater/scf.hpp for docs)
      
    // Perform the SCF
    void SCF(); 

    // Transformation functions to and from the orthonormal basis
    void ao2orthoFock();
    void ortho2aoDen();

    // Evaluate convergence
    bool evalConver();

    // Obtain new orbitals
    void getNewOrbitals(bool frmFock = true);

    // Misc procedural
    void diagOrthoFock();
    void printSCFProg(std::ostream &out = std::cout);
    void printSCFHeader(std::ostream &out = std::cout);
    void FDCommutator(oper_t_coll &);
    virtual void saveCurrentState();
    virtual void formDelta();

    // SCF extrapolation functions (see include/singleslater/extrap.hpp for docs)
    void allocExtrapStorage();
    void deallocExtrapStorage();
    void modifyFock();
    void fockDamping();
    void scfDIIS(size_t);

  }; // class SingleSlater

}; // namespace ChronusQ

#endif