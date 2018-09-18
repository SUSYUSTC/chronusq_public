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
#ifndef __INCLUDED_SINGLESLATER_HPP__
#define __INCLUDED_SINGLESLATER_HPP__

#include <chronusq_sys.hpp>
#include <wavefunction.hpp>
#include <singleslater/base.hpp>
#include <sjc_debug.hpp>

// Debug print triggered by Wavefunction
  
#ifdef _WaveFunctionDebug
  #define _SingleSlaterDebug
#endif

namespace ChronusQ {


  /**
   *  \brief The SingleSlater class. The typed abstract interface for all
   *  classes for which the wave function is described by a single slater
   *  determinant (HF, KS, PHF, etc).
   *
   *  Adds knowledge of storage type to SingleSlaterBase
   *
   *  Specializes the WaveFunction class of the same type
   */ 
  template <typename MatsT, typename IntsT>
  class SingleSlater : public SingleSlaterBase, public WaveFunction<MatsT,IntsT> {

  protected:

    // Useful typedefs
    typedef MatsT*                    oper_t;
    typedef std::vector<oper_t>       oper_t_coll;
    typedef std::vector<oper_t_coll>  oper_t_coll2;

  private:
  public:

    typedef MatsT value_type;
    typedef IntsT ints_type;

    //CORE_HAMILTONIAN_TYPE coreType;  ///< Core Hamiltonian type
    //ORTHO_TYPE            orthoType; ///< Orthogonalization scheme

    //OneETerms oneETerms; ///< One electron terms to be computed

    // Operator storage

    // AO Fock Matrix
    oper_t_coll fockMatrix;       ///< List of populated AO Fock matricies
    oper_t_coll fockMO;     ///< Fock matrix in the MO basis

    // Orthonormal Fock
    oper_t_coll fockMatrixOrtho;   ///< List of populated orthonormal Fock matricies

    // Coulomb (J[D])
    oper_t  coulombMatrix = nullptr;         ///< scalar Coulomb Matrix

    // Exchange (K[D])
    oper_t_coll exchangeMatrix;    ///< List of populated exact (HF) exchange matricies

    // Two-electron Hamiltonian (G[D])
    oper_t_coll twoeH;      ///< List of populated HF perturbation tensors

    // Orthonormal density
    oper_t_coll onePDMOrtho;   ///< List of populated orthonormal 1PDM matricies


    // Current / change in state information (for use with SCF)
    oper_t_coll curOnePDM;    ///< List of the current 1PDMs
    oper_t_coll deltaOnePDM;  ///< List of the changes in the 1PDMs

    // Stores the previous Fock matrix to use for damping    
    oper_t_coll prevFock;     ///< AO Fock from the previous SCF iteration

    // Stores the previous matrices to use for DIIS 
    oper_t_coll2 diisFock;    ///< List of AO Fock matrices for DIIS extrap
    oper_t_coll2 diisOnePDM;  ///< List of AO Density matrices for DIIS extrap
    oper_t_coll2 diisError;   ///< List of orthonormal [F,D] for DIIS extrap

    // 1-e integrals
    oper_t ortho1 = nullptr;   ///< Orthogonalization matrix which S -> I
    oper_t ortho2 = nullptr;   ///< Inverse of ortho1

    oper_t_coll coreH;          ///< Core Hamiltonian (scalar and magnetization)
    oper_t_coll coreHPerturbed; ///< Perturbed Core Hamiltonian (scalar and magnetization)

    // Method specific propery storage
    std::vector<double> mullikenCharges;
    std::vector<double> lowdinCharges;



    // Constructors
      
    /**
     *  SingleSlater Constructor. Constructs a SingleSlater object
     *
     *  \param [in] aoi  AOIntegrals object (which handels the BasisSet, etc)
     *  \param [in] args Parameter pack for the remaining parameters of the
     *                   WaveFunction constructor. See include/wavefunction.hpp
     *                   for details. 
     */ 
    template <typename... Args>
    SingleSlater(MPI_Comm c, AOIntegrals<IntsT> &aoi, Args... args) : 
      SingleSlaterBase(c,aoi.memManager(),args...), WaveFunctionBase(c,aoi.memManager(),args...),
      QuantumBase(c,aoi.memManager(),args...), WaveFunction<MatsT,IntsT>(c,aoi,args...)
      //, coreType(NON_RELATIVISTIC), orthoType(LOWDIN) 
    { 
      // Allocate SingleSlater Object
      alloc(); 

      // Determine Real/Complex part of method string
      if(std::is_same<MatsT,double>::value) {
        refLongName_  = "Real ";
        refShortName_ = "R-";
      } else {
        refLongName_  = "Complex ";
        refShortName_ = "C-";
      
      }

      // Default to NRH
      setCoreH(coreType);

    }; // SingleSlater constructor

    // See include/singleslater/impl.hpp for documentation 
    // on the following constructors

    // Different type
    template <typename MatsU> 
      SingleSlater(const SingleSlater<MatsU,IntsT> &, int dummy = 0);
    template <typename MatsU> 
      SingleSlater(SingleSlater<MatsU,IntsT> &&     , int dummy = 0);

    // Same type
    SingleSlater(const SingleSlater<MatsT,IntsT> &);
    SingleSlater(SingleSlater<MatsT,IntsT> &&);     

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


    // Declarations from QuantumBase 
    // (see include/singleslater/quantum.hpp for docs)
    void formDensity();
    void computeEnergy();
    void computeMultipole(EMPerturbation &);
    void computeSpin();

    // Compute various core Hamitlonian
    void formCoreH(EMPerturbation&,CORE_HAMILTONIAN_TYPE); // Compute the CH
    inline void formCoreH(EMPerturbation &emPert) { formCoreH(emPert,coreType); }
//  void updateCoreH(EMPerturbation &);
    void computeNRCH(EMPerturbation&,oper_t_coll&); // Non-relativistic CH
    void computeX2CCH(EMPerturbation&,oper_t_coll&); // X2C CH (aointegrals_rel.cxx)
    void compute4CCH(std::vector<libint2::Shell>&, dcomplex*); // 4C CH
    void computeOrtho();  // Evaluate orthonormalization transformations

    void addMagPert(EMPerturbation&, oper_t_coll&);

    // Method specific properties
    void populationAnalysis();
    void methodSpecificProperties() {
      populationAnalysis();
    }





    // Form a fock matrix (see include/singleslater/fock.hpp for docs)
    virtual void formFock(EMPerturbation &, bool increment = false, double xHFX = 1.);
    void formGD(EMPerturbation &, bool increment = false, double xHFX = 1.);

    // Form initial guess orbitals
    // see include/singleslater/guess.hpp for docs)
    void formGuess();
    void CoreGuess();
    void SADGuess();
    void RandomGuess();
    void ReadGuessMO();
    void ReadGuess1PDM();
    

    // Transformations to and from the orthonormal basis
    // see include/singleslater/ortho.hpp for docs
      
    template <typename TT> void Ortho1Trans(TT* A, TT* TransA); 
    template <typename TT> void Ortho2Trans(TT* A, TT* TransA); 
    template <typename TT> void Ortho1TransT(TT* A, TT* TransA);
    template <typename TT> void Ortho2TransT(TT* A, TT* TransA);


    // SCF procedural functions (see include/singleslater/scf.hpp for docs)

    // Transformation functions to and from the orthonormal basis
    void ao2orthoFock();
    void ortho2aoDen();
    void ortho2aoMOs();

    // Evaluate convergence
    bool evalConver(EMPerturbation &);

    // Obtain new orbitals
    void getNewOrbitals(EMPerturbation &, bool frmFock = true);
    void ConventionalSCF(bool modF);
    void NewtonRaphsonSCF();
    virtual MatsT* getNRCoeffs() = 0;

    // Misc procedural
    void diagOrthoFock();
    void FDCommutator(oper_t_coll &);
    virtual void saveCurrentState();
    virtual void formDelta();
    void orthoAOMO();
    void SCFInit();
    void SCFFin();

    // Stability and reopt
    virtual std::pair<double,MatsT*> getStab() = 0;
    bool checkStability();

    // Print functions
    void printFock(std::ostream& )    ;
    void print1PDMOrtho(std::ostream&);
    void printGD(std::ostream&)       ;
    void printJ(std::ostream&)        ;
    void printK(std::ostream&)        ;
    void printMiscProperties(std::ostream&);
    void printMOInfo(std::ostream&); 
    virtual void printFockTimings(std::ostream&);

    // SCF extrapolation functions (see include/singleslater/extrap.hpp for docs)
    void allocExtrapStorage();
    void deallocExtrapStorage();
    void modifyFock();
    void fockDamping();
    void scfDIIS(size_t);


	void swaporbit(bool, bool, std::pair<int,int>, std::pair<int,int>);
	void readpdm();
	void savenpy(std::string);
    // MO Transformations
    void MOFOCK();

  }; // class SingleSlater

}; // namespace ChronusQ


// Include headers for specializations of SingleSlater
#include <singleslater/hartreefock.hpp> // HF specialization
#include <singleslater/kohnsham.hpp>    // KS specialization

#endif
