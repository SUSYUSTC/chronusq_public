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

#ifndef __INCLUDED_REALTIME_ENUMS_HPP__
#define __INCLUDED_REALTIME_ENUMS_HPP__

namespace ChronusQ {

  enum IntegrationAlgorithm {
    MMUT,
    ExpMagnus2
  };

  enum PropagationStep {
    ForwardEuler,
    ModifiedMidpoint,
    ExplicitMagnus2,
	ExpotentialMM
  };

  enum PropagatorAlgorithm {
    Diagonalization,
    TaylorExpansion,
    ChebyshevExpansion
  };

  enum FieldEnvelopeTyp {
    Constant,
    LinRamp,
    Gaussian,
    Step
  };

};

#endif
