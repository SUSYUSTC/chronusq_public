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
#ifndef __INCLUDED_REALTIME_HPP__
#define __INCLUDED_REALTIME_HPP__

#include <chronusq_sys.hpp>
#include <cerr.hpp>
#include <memmanager.hpp>
#include <singleslater.hpp>


// RT Headers
#include <realtime/enums.hpp>
#include <realtime/fields.hpp>



namespace ChronusQ {

  /**
   *  \brief A struct to store information pertinent to the time
   *  propagation procedure.
   */ 
  struct IntegrationScheme {

    IntegrationAlgorithm intAlg  = MMUT;         ///< Integration Algorithm
    PropagationStep      rstStep = ForwardEuler; ///< Restart Step
    PropagatorAlgorithm  prpAlg  = Diagonalization; ///< exp(-iF) Algorithm

    double tMax    = 0.1;  ///< Max simulation time in AU
    double deltaT  = 0.01; ///< Time-step in AU

    size_t iRstrt  = 50; ///< Restart every N steps

  }; // struct IntegrationScheme

  /**
   *  \brief A struct to store information pertinent to the current
   *  state of the time propagation
   */ 
  struct IntegrationProgress {

    double  xTime = 0.; ///< Current time point
    size_t  iStep = 0;  ///< Step index of current time point
    double  stepSize;   ///< Current step size

    PropagationStep curStep;  ///< Current integration step

  };

  /**
   *  \brief A struct to store the property data obtained throughout the
   *  RealTime simulation
   */ 
  struct IntegrationData {

    std::vector<double> Time;
    std::vector<double> Energy;
    std::vector<std::array<double,3>> ElecDipole;

    // Field
    std::vector<std::array<double,3>> ElecDipoleField;
  };




  struct RealTimeBase {

    SafeFile savFile; ///< Data File

    IntegrationScheme intScheme;   ///< Integration scheme (MMUT, etc)
    TDEMPerturbation  pert;        ///< TD field perturbation

    IntegrationProgress curState;  ///< Current state of the time propagation
    IntegrationData     data;      ///< Data collection

	size_t scaling_first=5;
	bool is_swap_alpha=false;
	bool is_swap_beta=false;
	std::pair<int,int> swap_alpha;
	std::pair<int,int> swap_beta;
	size_t start_time=0;

    RealTimeBase()                     = delete;
    RealTimeBase(const RealTimeBase &) = delete;
    RealTimeBase(RealTimeBase &&)      = delete;


    RealTimeBase( CQMemManager &memManager): memManager_(memManager){ }


    // RealTimeBase procedural functions
    virtual void doPropagation()         = 0;

    // Progress functions
    void printRTHeader();
    void printRTStep();
    void appendStepRecord();


    /**
     *  \brief Adds a field to the time-dependent electromagnetic
     *  perturbation.
     *
     *  Calls TDEMPerturbation::addField. See include/realtime/fields.hpp
     *  for proper documentation.
     */ 
    template <typename... Args>
    inline void addField(Args... args){ pert.addField(args...); }

  protected:

    CQMemManager     &memManager_; ///< Memory manager

  public:
	bool restart=false;

  };


  template <template <typename, typename> class _SSTyp, typename IntsT>
  class RealTime : public RealTimeBase {
   
    typedef dcomplex*                 oper_t;
    typedef std::vector<oper_t>       oper_t_coll;

    SingleSlaterBase         *reference_ = nullptr;  ///< Initial conditions
    _SSTyp<dcomplex,IntsT>    propagator_; ///< Object for time propagation

    oper_t_coll DOSav;
    oper_t_coll UH;
    
  public:


    // Constructors

    // Disable default, copy and move constructors
    RealTime()                 = delete;
    RealTime(const RealTime &) = delete;
    RealTime(RealTime &&)      = delete;


    /**
     *  \brief RealTime Constructor.
     *
     *  Stores references to a "reference" SingleSlater object and
     *  CQMemManager and makes a copy of the reference into a complex
     *  SingleSlater object for the propagation.
     */ 
    template <typename RefMatsT>
    RealTime(_SSTyp<RefMatsT,IntsT> &reference) : 
      RealTimeBase(reference.memManager),
      reference_(&reference), propagator_(reference) { 

      alloc(); 

    }; // RealTime constructor
  
    ~RealTime(){ dealloc(); }


    // RealTime procedural functions
    void doPropagation(); // From RealTimeBase
    void formPropagator();
    void formFock(bool,double t);
    void propagateWFN();

    // Progress functions
    void printRTHeader();
    void printRTStep();
    void appendStepRecord();


    // Memory functions
    void alloc();
    void dealloc();

	//modify
	void DoIteration(EMPerturbation&, bool);
	void DoStart(EMPerturbation&);
	void DoMM(EMPerturbation&);

  }; // class RealTime
  

}; // namespace ChronusQ

#endif
