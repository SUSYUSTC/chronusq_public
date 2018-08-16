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
#ifndef __INCLUDED_AOINTEGRALS_HPP__
#define __INCLUDED_AOINTEGRALS_HPP__

#include <chronusq_sys.hpp>
#include <molecule.hpp>
#include <basisset/basisset_def.hpp>
#include <memmanager.hpp>
#include <libint2/engine.h>
#include <aointegrals/inhouseaointegral.hpp>
#include <util/files.hpp>
#include <fields.hpp>

namespace ChronusQ {

  /**
   *   \brief The AOIntegrals class. A class to handle the evaluation and
   *   storage of integrals in the atomic orbital (AO) gaussian--type
   *   orbital (GTO) basis. General functionality extends to both the 
   *   contracted (CGTO) and uncontracted (primitive, PGTO) bases.
   *
   *   Also handles the contraction of 2-body (3,4 index) integrals
   *   with 1-body (2 index) operators.
   *
   *   \warning Currently only functional for real GTOs.
   */

  enum TWOBODY_CONTRACTION_TYPE {
    COULOMB, ///< (mn | kl) X(lk)
    EXCHANGE,///< (mn | kl) X(nk)
    PAIR     ///< (mn | kl) X(nl)
  }; ///< 2-Body Tensor Contraction Specification


  enum CORE_HAMILTONIAN_TYPE {
    NON_RELATIVISTIC,
    RELATIVISTIC_X2C_SPIN_FREE,
    RELATIVISTIC_X2C_1E,
    RELATIVISTIC_X2C_2E,
    RELATIVISTIC_4C
  };

  struct OneETerms {
    bool finiteWidthNuc;
    bool coreH; //overlap, kinetic, potential
    bool relativistic; //spin-orbit, scalar relativity
  };

  /**
   *  The TwoBodyContraction struct. Stores information
   *  pertinant for a two body operator contraction with
   *  a one body (2 index) operator. z.B. The density matrix.
   */ 
  template <typename T>
  struct TwoBodyContraction {

    bool eval;

    T*  X;  ///< 1-Body (2 index) operator to contraction
    T*  AX; ///< 1-Body (2 index) storage for the contraction

    bool HER; ///< Whether or not X is hermetian
    
    TWOBODY_CONTRACTION_TYPE contType;

  }; // struct TwoBodyContraction

  enum CONTRACTION_ALGORITHM {
    DIRECT,
    INCORE,
    DENFIT
  }; ///< 2-e Integral Contraction Algorithm

  enum ORTHO_TYPE {
    LOWDIN,
    CHOLESKY
  }; ///< Orthonormalization Scheme





  /**
   *  \brief Abstract Base class for AOIntegrals
   *
   *  Stores type independent members and interfaces for templated the
   *  AOIntegrals class
   *
   */
  struct AOIntegralsBase { 

    size_t nTT_;  ///< Reduced number of basis functions \f$ N_B(N_B+1)/2 \f$
    size_t nSQ_;  ///< Squared basis functions \f$ N_B^2\f$
    size_t npTT_; ///< Reduced number of primitive functions \f$ N_P(N_P+1)/2 \f$
    size_t npSQ_; ///< Squared primitive functions \f$ N_P^2\f$

    CQMemManager &memManager_; ///< CQMemManager to allocate matricies
    Molecule     &molecule_;   ///< Molecule object for nuclear potential
    BasisSet     &basisSet_;   ///< BasisSet for the GTO basis defintion

    // Control Variables
    CONTRACTION_ALGORITHM contrAlg = DIRECT;///< Alg for 2-body contraction

    double threshSchwartz = 1e-12; ///< Schwartz screening threshold

    
    SafeFile savFile; ///< Hard storage of integrals

    // Default copy and move ctors
    AOIntegralsBase( const AOIntegralsBase & ) = default;
    AOIntegralsBase( AOIntegralsBase && )      = default;

    // Remove default ctor
    AOIntegralsBase() = delete;

    /**
     * \brief Constructor.
     *
     *  \param [in] memManager Memory manager for matrix allocation
     *  \param [in] mol        Molecule object for molecular specification
     *  \param [in] basis      The GTO basis for integral evaluation
     */
    AOIntegralsBase(CQMemManager &mem, Molecule &mol, BasisSet &basis) :
      memManager_(mem), molecule_(mol), basisSet_(basis) {

      nTT_  = basis.nBasis * ( basis.nBasis + 1 ) / 2;
      nSQ_  = basis.nBasis * basis.nBasis;
      npTT_ = basis.nPrimitive * ( basis.nPrimitive + 1 ) / 2;
      npSQ_ = basis.nPrimitive * basis.nPrimitive;

    }

    // Getters
    CQMemManager& memManager() { return memManager_; }
    BasisSet&     basisSet()   { return basisSet_;   }
    Molecule&     molecule()   { return molecule_;   }

    // Interfaces
    virtual void computeAOOneE(EMPerturbation&,OneETerms&) = 0; 
    virtual void computeERI(EMPerturbation&) = 0;    

    // Print (see src/aointegrals/print.cxx for docs)
    template <typename G> 
      friend std::ostream & operator<<(std::ostream &, const AOIntegralsBase& );

  };







  /**
   *  \brief Templated class to handle the evaluation and storage of 
   *  integral matrices representing quantum mechanical operators in
   *  a finite basis set.
   *
   *  Templated over storage type (IntsT) to allow for a seamless
   *  interface to both GTO and GIAO basis sets.
   */
  template <typename IntsT>
  class AOIntegrals : public AOIntegralsBase {

  private:

    typedef IntsT* oper_t; ///< Storage of an operator
    typedef std::vector<oper_t> oper_t_coll; ///< A collection of operators

    // General wrapper for 1-e integrals
    // See src/aointegrals/aointegrals_onee_drivers.cxx for documentation

    template <size_t NOPER, bool SYMM, typename F>
    std::vector<IntsT*> OneEDriverLocal(const F&,std::vector<libint2::Shell>&);
    template <size_t NOPER, bool SYMM, typename F>
    std::vector<IntsT*> OneEDriverLocalGTO(const F&,std::vector<libint2::Shell>&);
    template <size_t NOPER, bool SYMM, typename F>
    std::vector<IntsT*> OneEDriverLocalGIAO(const F&,std::vector<libint2::Shell>&);

    public:
    // 1-e builders for in-house integral code
    std::vector<IntsT*> OneEDriverLibint(libint2::Operator,std::vector<libint2::Shell>&);

    double* schwartz = nullptr; ///< Schwartz bounds for the ERIs

    oper_t overlap   = nullptr;   ///< Overlap matrix 
    oper_t kinetic   = nullptr;   ///< Kinetic matrix 
    oper_t potential = nullptr; ///< Nuclear potential matrix 

    oper_t_coll lenElecDipole;     ///< Electric Dipole matrix     (length)
    oper_t_coll lenElecQuadrupole; ///< Electric Quadrupole matrix (length)
    oper_t_coll lenElecOctupole;   ///< Electric Octuupole matrix  (length)

    oper_t_coll velElecDipole;     ///< Electric Dipole matrix     (velocity)
    oper_t_coll velElecQuadrupole; ///< Electric Quadrupole matrix (velocity)
    oper_t_coll velElecOctupole;   ///< Electric Octuupole matrix  (velocity)

    oper_t_coll magDipole;     ///< Electric Dipole matrix     (length)
    oper_t_coll magQuadrupole; ///< Electric Quadrupole matrix (length)

    // Relativistic integrals
    oper_t       PVdotP = nullptr;
    oper_t_coll  PVcrossP;  

    // 2-e Storage
    oper_t ERI = nullptr;    ///< Electron-Electron repulsion integrals (4 index) 

    // Constructors
    
    // Disable default constructor
    AOIntegrals() = delete;

    /**
     *  AOIntegrals Constructor. Constructs an AOIntegrals object.
     *
     *  \param [in] memManager Memory manager for matrix allocation
     *  \param [in] mol        Molecule object for molecular specification
     *  \param [in] basis      The GTO basis for integral evaluation
     */ 
    AOIntegrals(CQMemManager &memManager, Molecule &mol, BasisSet &basis) :
      AOIntegralsBase(memManager,mol,basis){ }

    // See src/aointegrals/aointegrals.cxx for documentation 
    // onf the following constructors

    AOIntegrals(const AOIntegrals &); // Copy constructor
    AOIntegrals(AOIntegrals &&);      // Move constructor

    /**
     *  Copy assignment.
     *
     *  Assign one AOIntegrals object to another by a copy
     */ 
    AOIntegrals& operator=(const AOIntegrals &) = default;
    
    /**
     *  Move assignment.
     *
     *  Assign one AOIntegrals object to another by a move
     */ 
    AOIntegrals& operator=(AOIntegrals &&) = default;

    /**
     *  Destructor.
     *
     *  Destructs an AOIntegrals object
     */ 
    ~AOIntegrals() { dealloc(); }
    

    // Member functions


    // Memory

    // Deallocation (see src/aointegrals/aointegrals.cxx for docs)
    void dealloc();

    // Integral evaluation
    // (see src/aointegrals/aointegrals_onee/twoe_drivers.cxx for docs)

    void computeAOOneE(EMPerturbation&,OneETerms&);     // Evaluate the 1-e ints (general)
    void computeAOOneEGTO(OneETerms&);                  // Evaluate the 1-e ints in the CGTO basis
    void computeAOOneEGIAO(EMPerturbation&,OneETerms&); // Evaluate the 1-e ints in the GIAO basis

    void computeERI(EMPerturbation&);                   // Evaluate ERIs (general)
    void computeERIGTO();                               // Evaluate ERIs in the CGTO basis
    void computeERIGIAO(EMPerturbation&);               // Evaluate ERIs in the GIAO basis
    void computeSchwartz();                             // Evaluate schwartz bounds (currently implemented for CGTOs 


    // Integral contraction

    /**
     *  Contract the two body potential with one body (2 index) operators.
     *
     *  Smartly determines whether to do the contraction directly, incore
     *  or using density fitting depending on context
     *
     *  \param [in/ont] contList List of one body operators for contraction.
     */ 
    template <typename TT>
    void twoBodyContract(
        MPI_Comm comm,
        const bool screen,
        std::vector<TwoBodyContraction<TT>> &contList, 
        EMPerturbation &pert) {

      if( contrAlg == INCORE ) twoBodyContractIncore(comm,contList);
      else if( contrAlg == DIRECT ) 
        twoBodyContractDirect(comm,screen,contList,pert);

    };
    template <typename TT>
    void twoBodyContract(
        MPI_Comm comm,
        const bool screen,
        std::vector<TwoBodyContraction<TT>> &contList) { 

      EMPerturbation pert;
      twoBodyContract(comm,screen,contList,pert);


    };
    
    template <typename TT>
    inline void twoBodyContract(
        MPI_Comm comm, 
        std::vector<TwoBodyContraction<TT>> &contList, 
        EMPerturbation &pert) {

      twoBodyContract(comm,true,contList,pert);

    }

    template <typename TT>
    inline void twoBodyContract(
        MPI_Comm comm, 
        std::vector<TwoBodyContraction<TT>> &contList) { 

      twoBodyContract(comm,true,contList);

    }

    // INCORE contraction routines
    // Perform the two body contraction incore (using the rank-4 ERI tensor)
    // see include/aointegrals/contract/incore.hpp for docs.
      
    template <typename TT>
    void twoBodyContractIncore(MPI_Comm, std::vector<TwoBodyContraction<TT>>&);

    template <typename TT>
    void JContractIncore(MPI_Comm, TwoBodyContraction<TT> &);

    template <typename TT>
    void KContractIncore(MPI_Comm, TwoBodyContraction<TT> &);



    // DIRECT contraction routines
    // Perform the two body contraction directly
    // see include/aointegrals/contract/direct.hpp for docs.
      
    template <typename TT>
    void twoBodyContractDirect(
        MPI_Comm,
        const bool,
        std::vector<TwoBodyContraction<TT>>&, 
        EMPerturbation &pert);

    template <typename TT>
    void directScaffold(
        MPI_Comm,
        const bool,
        std::vector<TwoBodyContraction<TT>>&);

    template <typename TT>
    void directScaffold(
        MPI_Comm,
        const bool,
        std::vector<TwoBodyContraction<TT>>&, 
        EMPerturbation &pert);

    template <typename TT>
    void JContractDirect(
        MPI_Comm,
        const bool,
        TwoBodyContraction<TT> &);

    template <typename TT>
    void KContractDirect(
        MPI_Comm,
        const bool,
        TwoBodyContraction<TT> &);


    
  }; // class AOIntegrals

  extern std::array<std::array<double,25>,3201> FmTTable;

  void generateFmTTable();  


}; // namespace ChronusQ

#endif
