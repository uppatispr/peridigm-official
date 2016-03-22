/*! \file Peridigm_Compute_CrackTip_Location.hpp */

#ifdef COMPUTE_CLASS

ComputeClass(CrackTip_Location,Compute_CrackTip_Location)

#else

#ifndef PERIDIGM_COMPUTE_CRACKTIP_LOCATION_HPP
#define PERIDIGM_COMPUTE_CRACKTIP_LOCATION_HPP

#include "Peridigm_Compute.hpp"

namespace PeridigmNS {

  //! Base class for determing the crack tip location for 1-D crack growth
  class Compute_CrackTip_Location : public PeridigmNS::Compute {

  public:
    //! Standard constructor.
    Compute_CrackTip_Location( Teuchos::RCP<const Teuchos::ParameterList> params,
                             Teuchos::RCP<const Epetra_Comm> epetraComm_,
                             Teuchos::RCP<const Teuchos::ParameterList> computeClassGlobalData_);

    //! Destructor.
    virtual  ~Compute_CrackTip_Location();

    //! Returns a vector of field IDs corresponding to the variables associated with the compute class.
    virtual std::vector<int> FieldIds() const { return m_fieldIds; }

    //! Perform computation.
    virtual int compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks  ) const;

    //! Compute the crack tip location. 
    int computeCrackTip( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks, bool storeLocal ) const;

  private:

    // field ids for all relevant data
    std::vector<int> m_fieldIds;
    int m_crackTipFieldId;
    int m_damageFieldId;
    int m_modelCoordinatesFieldId;
    int m_timeStepCountFieldId;		// INPRINCE EDIT: timeStepCount -> which time step is the simulation on? 
    int m_horizonMultiple;		// INPRINCE EDIT: horizon multiple as an upper limit for the change in crack tip location in one time step
    int m_initialCTX;			// INPRINCE EDIT: initial position of crack tip location (x-coordinate)
    int m_initialCTY;			// INPRINCE EDIT: initial position of crack tip location (y-coordinate)
    int m_initialCTZ;			// INPRINCE EDIT: initial position of crack tip location (z-coordinate)

    // Damage criteria
    double m_damageLowerLimit;		// (lower damage limit) Some percentage damage above which damage is likely caused by fracture propagation
    double m_damageUpperLimit;		// (upper damage limit) Damage limit, above which the nodes are disconnected from body due to instability

    // Direction of major crack growth
    int m_growthDirection;		// Growth Direction (1 - xP, 2 - yP, 3 - zP, 4 - xN, 5 - yN, 6 - zN)

    // Label for output variable
    std::string m_outputLabel;
  };
}

#endif // PERIDIGM_COMPUTE_KINETIC_ENERGY_HPP
#endif // COMPUTE_CLASS
