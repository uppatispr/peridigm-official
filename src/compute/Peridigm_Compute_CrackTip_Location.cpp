/*! \file Peridigm_Compute_CrackTip_Location.cpp */

#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <ctgmath>

#include "Peridigm_Compute_CrackTip_Location.hpp"
#include "Peridigm_Field.hpp"
#include "Peridigm_HorizonManager.hpp"

//! Standard constructor.
PeridigmNS::Compute_CrackTip_Location::Compute_CrackTip_Location(Teuchos::RCP<const Teuchos::ParameterList> params,
                                                           Teuchos::RCP<const Epetra_Comm> epetraComm_,
                                                           Teuchos::RCP<const Teuchos::ParameterList> computeClassGlobalData_)
  : Compute(params, epetraComm_, computeClassGlobalData_), m_crackTipFieldId(-1), m_damageFieldId(-1), m_modelCoordinatesFieldId(-1),
    m_timeStepCountFieldId(-1)									// INPRINCE EDIT: initializing timeStepCountFieldId
{
  if (params->isParameter("Damage Lower Limit"))
	m_damageLowerLimit = params->get<double>("Damage Lower Limit");					
  else
	m_damageLowerLimit = 0.25;

  if (params->isParameter("Damage Upper Limit"))
	m_damageUpperLimit = params->get<double>("Damage Upper Limit");
  else
	m_damageUpperLimit = 0.8;

  m_growthDirection = params->get<int>("Crack Growth Direction");

  if ((m_growthDirection != 1) && (m_growthDirection != 2) && (m_growthDirection != 3) && (m_growthDirection != 4) && (m_growthDirection != 5) && (m_growthDirection != 6))
	TEUCHOS_TEST_FOR_EXCEPT_MSG(true, "**** Error: CrackTip_Location compute class can only accept crack growth directions from the set {1, 2, 3, 4, 5, 6} corresponding to the directions {xP, yP, zP, xN, yN, zN}.\n");

  m_horizonMultiple = params->get<double>("Horizon Multiple");			// INPRINCE EDIT: multiple of horizon as an upper limit of the change in crack tip location in one time step

  m_outputLabel = params->get<std::string>("Output Label");

  m_initialCTX = params->get<double>("Initial X");				// INPRINCE EDIT: initiat known x-coordinate location of the crack tip (due to precrack)
  m_initialCTY = params->get<double>("Initial Y");				// INPRINCE EDIT: initiat known y-coordinate location of the crack tip (due to precrack)
  m_initialCTZ = params->get<double>("Initial Z");				// INPRINCE EDIT: initiat known z-coordinate location of the crack tip (due to precrack)

  FieldManager& fieldManager = FieldManager::self();
  m_crackTipFieldId = fieldManager.getFieldId(PeridigmField::GLOBAL, PeridigmField::VECTOR, PeridigmField::TWO_STEP, "CrackTip_Location");
  m_damageFieldId = fieldManager.getFieldId(PeridigmField::ELEMENT, PeridigmField::SCALAR, PeridigmField::TWO_STEP, "Damage");
  m_modelCoordinatesFieldId = fieldManager.getFieldId(PeridigmField::NODE, PeridigmField::VECTOR, PeridigmField::CONSTANT, "Model_Coordinates");
  m_timeStepCountFieldId = fieldManager.getFieldId(PeridigmField::GLOBAL, PeridigmField::SCALAR, PeridigmField::CONSTANT, "TimeStep_Count");	// INPRINCE EDIT: get actual vlaue of timeStepCountFieldId
  m_fieldIds.push_back(m_crackTipFieldId);
  m_fieldIds.push_back(m_damageFieldId);
  m_fieldIds.push_back(m_modelCoordinatesFieldId);
  m_fieldIds.push_back(m_timeStepCountFieldId);		// INPRINCE EDIT: add timeStepCountFieldId to the list of fieldIds used in this compute class
}

//! Destructor.
PeridigmNS::Compute_CrackTip_Location::~Compute_CrackTip_Location(){}

//! Perform computation. 
int PeridigmNS::Compute_CrackTip_Location::compute( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks ) const
{
  	bool storeLocal = true;
  	int result = computeCrackTip(blocks, storeLocal);
  	return result;
}

//! Computation function definition. 
int PeridigmNS::Compute_CrackTip_Location::computeCrackTip( Teuchos::RCP< std::vector<PeridigmNS::Block> > blocks, bool storeLocal ) const
{ 
	int retval;
	
	// STEP 1: COLLECT DAMAGED NODES AND CORRESPONDING COORIDNATES AS CANDIDATES FOR CRACK TIP NODES AND COORDINATE LOCATIONS
 
	// INPRINCE EDIT: Determine the current time step number
	Teuchos::RCP<Epetra_Vector> timeStepCount_data = blocks->begin()->getData(m_timeStepCountFieldId, PeridigmField::STEP_NONE);
	int stepCount = (*timeStepCount_data)[0];
	std::cout << std::endl << "Time step number: " << stepCount << std::endl;		// INPRINCE EDIT

	// Define variables and assign them based on which direction is relevant for tracking crack tip and whether it's positive or negative
	int xyzDir;
	bool dirPos;

	if (m_growthDirection == 1 || m_growthDirection == 4)
		xyzDir = 0;
	else if (m_growthDirection == 2 || m_growthDirection == 5)
		xyzDir = 1;
	else if (m_growthDirection == 3 || m_growthDirection == 6)
		xyzDir = 2;

	if (m_growthDirection == 1 || m_growthDirection == 2 || m_growthDirection == 3)
		dirPos = true;
	else
		dirPos = false;

	// Declare necessary entities 
	PeridigmNS::HorizonManager& horizonManager = PeridigmNS::HorizonManager::self();		// INPRINCE EDIT: declare the horizon manager
	std::string blockName;
	double blockHorizon;
	double CTn[3];											// INPRINCE EDIT: previous time step CT_X, CT_Y, CT_Z values declared
	Teuchos::RCP<Epetra_Vector> prevCT, damage, xyz_coord;
	std::vector<int> damaged_nodeIds;
	std::vector<double> damaged_coord, damaged_x, damaged_y, damaged_z;
	std::vector<Block>::iterator blockIt;

	// INPRINCE EDIT: If it is not the initial time (step 0) or the first simulation time step (step 1), determine the previous time step crack tip coordinates
	if (stepCount > 1) {
		prevCT = blocks->begin()->getData(m_crackTipFieldId, PeridigmField::STEP_N);
		CTn[0] = (*prevCT)[0];
		CTn[1] = (*prevCT)[1];
		CTn[2] = (*prevCT)[2];
		std::cout << "Previous CT location (at time step " << stepCount-1 << "): (" << CTn[0] << ", " << CTn[1] << ", " << CTn[2] << ")" << std::endl;	// INPRINCE EDIT
	}

	// INPRINCE EDIT: Iterate over the blocks
	for(blockIt = blocks->begin() ; blockIt != blocks->end() ; blockIt++)
	{
		// INPRINCE EDIT: Determine horizon for block
		blockName = blockIt->getName();
		blockHorizon = horizonManager.getBlockConstantHorizonValue(blockName);	

		Teuchos::RCP<NeighborhoodData> neighborhoodData = blockIt->getNeighborhoodData();
		const int numOwnedPoints = neighborhoodData->NumOwnedPoints();
		const int* ownedIDs = neighborhoodData->OwnedIDs();

		damage = blockIt->getData(m_damageFieldId, PeridigmField::STEP_NP1);
		xyz_coord = blockIt->getData(m_modelCoordinatesFieldId, PeridigmField::STEP_NONE);
		
		// Collect values
		double *damage_values = damage->Values();
		double *xyz_coord_vals = xyz_coord->Values();

		// on a per-element basis, check the damage at each element/node to see if it falls in a certain range
		// If damage criteria is satisfied and is close enough to previous crack tip coordinates, record the nodeID of the damaged node and its coordinates
		int numElements = numOwnedPoints;
		double damage_atNode;
		double x_coord, y_coord, z_coord, relevantCoord;

		for (int i=0;i<numElements;i++) 
		{
			int ID = ownedIDs[i];
			damage_atNode = damage_values[ID];
			bool isAppropriateNode = false;
		
			if (damage_atNode >= m_damageLowerLimit && damage_atNode <= m_damageUpperLimit) {
				x_coord = xyz_coord_vals[3*ID];
				y_coord = xyz_coord_vals[3*ID+1];
				z_coord = xyz_coord_vals[3*ID+2];
				relevantCoord = x_coord + (y_coord-x_coord)*(xyzDir==1) + (z_coord-x_coord)*(xyzDir==2);

				if (stepCount <= 1) 
					isAppropriateNode = true;
				else if (sqrt((x_coord-CTn[0])*(x_coord-CTn[0])+(y_coord-CTn[1])*(y_coord-CTn[1])+(z_coord-CTn[2])*(z_coord-CTn[2])) < m_horizonMultiple*blockHorizon)
					isAppropriateNode = true;
			}
			
			if (isAppropriateNode) {
				damaged_nodeIds.push_back(ID);
				damaged_x.push_back(x_coord);
				damaged_y.push_back(y_coord);
				damaged_z.push_back(z_coord);
				damaged_coord.push_back(relevantCoord);
			}
		}
	}

	// INPRINCE EDIT: PRINT OUT DAMAGED COORD LIST
	/*std::vector<double>::iterator coordIt;
	std::cout << std::endl << "Damaged x-coordinates of nodes: ";
	for (coordIt = damaged_coord.begin(); coordIt != damaged_coord.end(); coordIt++)
		std::cout << *coordIt << ", ";	
	std::cout << std:: endl;
	*/

	// Sanity check: make sure the number of coordinate positions collected matches the number of damaged nodes
	if (damaged_nodeIds.size() != damaged_coord.size())
	 {
		retval = 1;
		return (retval);
	 }	

	// STEP 2: FIND LOCAL MAXIMA OR MINIMA FOR CRACK TIP CANDIDATES ON EACH PROCESSOR (IF ANY) 
	
	// When there are no damaged nodes, we just use the initial crack tip location as the crack has not propagated
	int damage_indicator = 0;								// Damage boolean (defined as integer for EpetraComm comparison)
	if (damaged_coord.size() != 0)
		damage_indicator = 1;

	// Once all the appropriate (initial crack related) damaged node coordinates (in appropriate provided crack direction) are collected from
	// all blocks in a single processor, find the largest or smallest based on the growth direction.
	double extrema;

	if (damage_indicator == 1) 
	 {
		std::vector<double>::iterator extrema_index;
		
		if (dirPos)
	 	 {
			extrema_index = std::max_element(damaged_coord.begin(), damaged_coord.end());	
			extrema = *extrema_index;							// local max x/y/z of damaged nodes
	 	 }
		else
	 	 {
			extrema_index = std::min_element(damaged_coord.begin(), damaged_coord.end());
			extrema = *extrema_index;							// local min x/y/z of damaged nodes
	 	 }
	 }
	else
	 {
		// If damage doesn't exist on this processor, assign extrema value such that comparison across processors isn't affected by this value
		if (dirPos)
			extrema = -std::numeric_limits<double>::max();
		else
			extrema = std::numeric_limits<double>::max();
	 }

	// STEP 3: FIND CRACK TIP COORDINATES FROM ALL PROCESSOR COMPARISON AND STORE AS A GLOBAL QUANITITY FOR CURRENT TIME STEP
	bool noCrack = true;			// If boolean quantity is true, there are no damaged nodes anywhere and crack tip is just the initial location

	// Compare the damage_indicators on all the processors and determine if there is a crack under any processor
	int local_indicator, global_indicator;
	local_indicator = damage_indicator;

	epetraComm()->SumAll(&local_indicator, &global_indicator, 1);

	if (global_indicator != 0) 
		noCrack = false;

	// Determind crack tip coordinates
	double CT[3]; 										// Crack tip x, y and z coordinate declarations
	int globalExtrema_PID = std::numeric_limits<int>::max();
	bool globalExtremaNode_on_processor = false;						// Is the damaged node with largest coordinate position
												// in the crack growth direction on this processor? 

	if (!noCrack)										// If there is crack damage on at least on processor
	 {
		double local_extrema, global_extrema;
		local_extrema = extrema;

		if (m_growthDirection == 1 || m_growthDirection == 2 || m_growthDirection == 3)	
			epetraComm()->MaxAll(&local_extrema, &global_extrema, 1);
		else
			epetraComm()->MinAll(&local_extrema, &global_extrema, 1);


		// Determine if extrema is on this processor and find its x,y,z-coordinates
		std::vector<double>::iterator extremaNode_finder;
		std::vector<double>::iterator extremaNode_coord;
		int dist_to_extrema;
	
		if (global_extrema == extrema)			// If global_extrema matches local_extrema
	 	 {
			globalExtremaNode_on_processor = true;
			globalExtrema_PID = epetraComm()->MyPID();
			extremaNode_finder = find(damaged_coord.begin(), damaged_coord.end(), extrema);
			dist_to_extrema = std::distance(damaged_coord.begin(), extremaNode_finder);

			extremaNode_coord = damaged_x.begin() + dist_to_extrema;
			CT[0] = *extremaNode_coord;								// X-coordinate of the crack tip
		
			extremaNode_coord = damaged_y.begin() + dist_to_extrema;
			CT[1] = *extremaNode_coord;								// Y-coordinate of the crack tip

			extremaNode_coord = damaged_z.begin() + dist_to_extrema;
			CT[2] = *extremaNode_coord;								// Z-coordinate of the crack tip
	 	 }
	 }
	else
	 {
		globalExtrema_PID = epetraComm()->MyPID();
		CT[0] = m_initialCTX;
		CT[1] = m_initialCTY;
		CT[2] = m_initialCTZ;
	 }

	// Synchronize the processors
	epetraComm()->Barrier();	

	// Determind the PID of the processor that contains the global_extrema
	int localPID, minMatchingPID;		
	localPID = globalExtrema_PID;
	epetraComm()->MinAll(&localPID, &minMatchingPID, 1);

	// Broadcast this value to all processors
	epetraComm()->Broadcast(&CT[0], 3, minMatchingPID);		

	// INPRINCE EDIT: Change CT values at timestep 2 for diagnostic purposes:
	if (stepCount == 2) {
		CT[0] = 1;
		CT[1] = 2;
		CT[2] = 3;
	}

	// INPRINCE EDIT: see what CT values are being assigned
	std::cout << std::endl << "Current processor CT location at time step " << stepCount << ": (" << CT[0] << ", " << CT[1] << ", " << CT[2] << ")" << std::endl;

	// Store this coordinate information (store in first block, since block globals are static)
	Teuchos::RCP<Epetra_Vector> data = blocks->begin()->getData(m_crackTipFieldId, PeridigmField::STEP_NP1);	
	(*data)[0] = CT[0];
	(*data)[1] = CT[1];
	(*data)[2] = CT[2];

	return(0);

}
