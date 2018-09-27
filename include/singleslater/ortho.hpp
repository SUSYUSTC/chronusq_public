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
#ifndef __INCLUDED_SINGLESLATER_ORTHO_HPP__
#define __INCLUDED_SINGLESLATER_ORTHO_HPP__

#include <singleslater.hpp>
#include <cqlinalg/blas3.hpp>


namespace ChronusQ {

  /**
   *  \brief Performs the transformation \f$ A' = O_1 A O_1^T\f$ for
   *  a general matrix \f$A\f$ 
   *
   *  \f$ O_1 \f$ is the orthonormalization matrix stored in 
   *  AOIntegrals::ortho1 (see AOIntegrals::ORTHO_TYPE and
   *  AOIntegrals::computeOrtho for details).
   *
   */ 
  template <typename MatsT, typename IntsT>
  template <typename TT> 
  void SingleSlater<MatsT,IntsT>::Ortho1Trans(TT* A, TT* TransA) {

    size_t NB = this->aoints.basisSet().nBasis;
    size_t nSQ = NB*NB;

    // Make sure that the incoming matricies are of the right size
    assert(memManager.template getSize(A) == nSQ);
    assert(memManager.template getSize(TransA) == nSQ);

    // Allocate scratch space
    TT* SCR = memManager.template malloc<TT>(nSQ);

    // If necessary, create a complex copy of ortho1
    TT* O1 = reinterpret_cast<TT*>(ortho1);
    if(std::is_same<TT,dcomplex>::value) {
      O1 = memManager.template malloc<TT>(nSQ);
      std::copy_n(ortho1,nSQ,O1);
    }

    // Perform transformation
    Gemm('N', 'N', NB, NB, NB, TT(1.), O1, NB, A, NB, TT(0.), SCR, NB);
    Gemm('N', 'C', NB, NB, NB, TT(1.), SCR, NB, O1, NB, TT(0.), TransA, NB);

    // Free up scratch space
    memManager.free(SCR);
    if(std::is_same<TT,dcomplex>::value) memManager.free(O1);
    
  }; // AOIntegrals::Ortho1Trans

  /**
   *  \brief Performs the transformation \f$ A' = O_2 A O_2^T\f$ for
   *  a general matrix \f$A\f$ 
   *
   *  \f$ O_2 \f$ is the orthonormalization matrix stored in 
   *  AOIntegrals::ortho1 (see AOIntegrals::ORTHO_TYPE and
   *  AOIntegrals::computeOrtho for details).
   *
   */ 
  template <typename MatsT, typename IntsT>
  template <typename TT> 
  void SingleSlater<MatsT,IntsT>::Ortho2Trans(TT* A, TT* TransA) {

    size_t NB = this->aoints.basisSet().nBasis;
    size_t nSQ = NB*NB;
    // Make sure that the incoming matricies are of the right size
    assert(memManager.template getSize(A) == nSQ);
    assert(memManager.template getSize(TransA) == nSQ);

    // Allocate scratch space
    TT* SCR = memManager.template malloc<TT>(nSQ);

    // If necessary, create a complex copy of ortho2
    TT* O2 = reinterpret_cast<TT*>(ortho2);
    if(std::is_same<TT,dcomplex>::value) {
      O2 = memManager.template malloc<TT>(nSQ);
      std::copy_n(ortho2,nSQ,O2);
    }

    // Perform transformation
    Gemm('N', 'N', NB, NB, NB, TT(1.), O2, NB, A, NB, TT(0.), SCR, NB);
    Gemm('N', 'C', NB, NB, NB, TT(1.), SCR, NB, O2, NB, TT(0.), TransA, NB);

    // Free up scratch space
    memManager.free(SCR);
    if(std::is_same<TT,dcomplex>::value) memManager.free(O2);
    
  }; // AOIntegrals::Ortho2Trans

  /**
   *  \brief Performs the transformation \f$ A' = O_1^T A O_1\f$ for
   *  a general matrix \f$A\f$ 
   *
   *  \f$ O_1 \f$ is the orthonormalization matrix stored in 
   *  AOIntegrals::ortho1 (see AOIntegrals::ORTHO_TYPE and
   *  AOIntegrals::computeOrtho for details).
   *
   */ 
  template <typename MatsT, typename IntsT> 
  template <typename TT>
  void SingleSlater<MatsT,IntsT>::Ortho1TransT(TT* A, TT* TransA) {

    size_t NB = this->aoints.basisSet().nBasis;
    size_t nSQ = NB*NB;

    // Make sure that the incoming matricies are of the right size
    assert(memManager.template getSize(A) == nSQ);
    assert(memManager.template getSize(TransA) == nSQ);

    // Allocate scratch space
    TT* SCR = memManager.template malloc<TT>(nSQ);

    // If necessary, create a complex copy of ortho1
    TT* O1 = reinterpret_cast<TT*>(ortho1);
    if(std::is_same<TT,dcomplex>::value) {
      O1 = memManager.template malloc<TT>(nSQ);
      std::copy_n(ortho1,nSQ,O1);
    }

    // Perform transformation
    Gemm('C', 'N', NB, NB, NB, TT(1.), O1, NB, A, NB, TT(0.), SCR, NB);
    Gemm('N', 'N', NB, NB, NB, TT(1.), SCR, NB, O1, NB, TT(0.), TransA,
      NB);

    // Free up scratch space
    memManager.free(SCR);
    if(std::is_same<TT,dcomplex>::value) memManager.free(O1);
    
  }; // AOIntegrals::Ortho1TransT

  /**
   *  \brief Performs the transformation \f$ A' = O_2^T A O_2\f$ for
   *  a general matrix \f$A\f$ 
   *
   *  \f$ O_2 \f$ is the orthonormalization matrix stored in 
   *  AOIntegrals::ortho1 (see AOIntegrals::ORTHO_TYPE and
   *  AOIntegrals::computeOrtho for details).
   *
   */ 
  //Added by Jiace
  template <typename MatsT, typename IntsT>
  template <typename TT> 
  void SingleSlater<MatsT,IntsT>::Ortho2TransT(TT* A, TT* TransA) {

    size_t NB = this->aoints.basisSet().nBasis;
    size_t nSQ = NB*NB;
    // Make sure that the incoming matricies are of the right size
    assert(memManager.template getSize(A) == nSQ);
    assert(memManager.template getSize(TransA) == nSQ);

    // Allocate scratch space
    TT* SCR = memManager.template malloc<TT>(nSQ);

    // If necessary, create a complex copy of ortho2
    TT* O2 = reinterpret_cast<TT*>(ortho2);
    if(std::is_same<TT,dcomplex>::value) {
      O2 = memManager.template malloc<TT>(nSQ);
      std::copy_n(ortho2,nSQ,O2);
    }

    // Perform transformation
    Gemm('C', 'N', NB, NB, NB, TT(1.), O2, NB, A, NB, TT(0.), SCR, NB);
    Gemm('N', 'N', NB, NB, NB, TT(1.), SCR, NB, O2, NB, TT(0.), TransA, NB);

    // Free up scratch space
    memManager.free(SCR);
    if(std::is_same<TT,dcomplex>::value) memManager.free(O2);
    
  }; // AOIntegrals::Ortho2Trans

}; // namespace ChronusQ

#endif
