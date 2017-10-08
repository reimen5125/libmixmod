/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/BinaryParameter.h  description
    copyright            : (C) MIXMOD Team - 2001-2016
    email                : contact@mixmod.org
 ***************************************************************************/

/***************************************************************************
    This file is part of MIXMOD
    
    MIXMOD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MIXMOD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MIXMOD.  If not, see <http://www.gnu.org/licenses/>.

    All informations available on : http://www.mixmod.org                                                                                               
***************************************************************************/
#ifndef XEMBinaryParameter_H
#define XEMBinaryParameter_H

#include "mixmod/Utilities/Util.h"
#include "mixmod/Kernel/Parameter/Parameter.h"

namespace XEM {

class GaussianParameter;
class BinaryParameter;

/**
  @brief Base class for XEMBinaryParameter(s)
  @author F. Langrognet
 */

class BinaryParameter : public Parameter {

public:

	//----------------------------
	// constructors / desctructors
	// ---------------------------

	/// Default constructor
	BinaryParameter();

	/// Constructor
	// called by XEMModel (via XEMBinary...Parameter)
	BinaryParameter(Model * iModel, ModelType * iModelType, int * tabNbModality);

	// constructor
	// called if USER initialisation
	BinaryParameter(int iNbCluster, int iPbDimension, 
			ModelType * iModelType, int * tabNbModality);

	/// Constructor
	BinaryParameter(const BinaryParameter * iParameter);

	/// Destructor
	virtual ~BinaryParameter();

	/// Comparison operator
	virtual bool operator ==(const BinaryParameter & param) const;


	/// reset to default values
	virtual void reset();

	/// create the same parameter than this but after updating because without xi0
	Parameter * createParameter(Model * iModel, int i0, int ki0);


	//----------
	// selectors
	//----------

	/// get TabCenter
	int ** getTabCenter() const;

	/// get _tabNbModality
	int * getTabNbModality() const;

	/// get total number of modality
	int getTotalNbModality() const;

	//----------------
	// compute methods
	//----------------

	void getAllPdf(float ** tabFik, float * tabProportion) const;

	/** @brief Compute probability density
	@param iSample  Probability for sample iSample
	@param kCluster Probability in class kCluster
	 */

	virtual float getPdf(int iSample, int kCluster) const = 0;

	/** @brief Compute log probability density
	@param iSample  Probability for sample iSample
	@param kCluster Probability in class kCluster
	 */
	virtual float getLogPdf(int iSample, int kCluster) const = 0;

	/** Compute normal probability density function
		 for x vector and kCluster th cluster
	 */
	//float getPdf(RowVector x, int kCluster);
	virtual float getPdf(Sample * x, int kCluster) const = 0;

	/// getlogLikelihoodOne (one cluster)
	virtual float getLogLikelihoodOne() const = 0;

	/// compute Tik for xi (i=0 -> _nbSample-1) when underflow
	virtual void computeTikUnderflow(int i, float ** tabTik);

	/// Compute table of centers of the samples for each cluster
	void computeTabCenter();

	/// Compute scatter(s) 
	virtual void computeScatter() = 0;

	/// Compute random scatter(s)
	virtual void computeRandomScatter() = 0;

	///recopy scatter from param (used for init  : USER)
	virtual void recopyScatter(Parameter * iParam) = 0;


	//---------------
	// initialization
	//---------------

	/// init user
	void initUSER(Parameter * iParam);

	/// updateForInitRANDOMorUSER_PARTITION
	void updateForInitRANDOMorUSER_PARTITION(
			Sample ** tabSampleForInit, bool * tabClusterToInitialize);

	/// initialize attributes before an InitRandom 
	void initForInitRANDOM();

	/// initialize attributes for init USER_PARTITION
	/// outputs :
	/// -  nbInitializedCluster
	/// - tabNotInitializedCluster (array of size _nbCluster)
	void initForInitUSER_PARTITION(int & nbInitializedCluster, 
			bool * tabNotInitializedCluster, Partition * initPartition);

	/// computeTabCenterInitUSER_PARTITIONoutputs :
	/// -  nbInitializedCluster
	/// - tabNotInitializedCluster (array of size _nbCluster)
	void computeTabCenterInitUSER_PARTITION(int & nbInitializedCluster, 
			bool * tabNotInitializedCluster, Partition * initPartition);


	//-----------
	// Algorithms
	//-----------

	/// Maximum a posteriori step method
	void MAPStep();

	/// Expectation step method
	//void EStep();

	/// Maximization step method
	void MStep();


	//---------------
	// input / output
	//---------------

	// edit (for debug)
	void edit();

	/// editScatter (for debug)
	virtual void editScatter(int k) = 0;

	/// Edit
	void edit(std::ofstream & oFile, bool text = false);

	/// editScatter 
	virtual void editScatter(std::ofstream & oFile, int k, bool text = false) = 0;

	// Read Parameters in input file
	void input(std::ifstream & fi);

	// Read Parameters in input containers
	void input(
			float * proportions, 
			float ** centers, 
			float *** scatters);

	// Read Scatter in input file
	virtual void inputScatter(std::ifstream & fi, int k) = 0;
	virtual void inputScatter(float *** scatters) = 0;

	/// recopie sans faire construction / destruction
	// utilise par SMALL_EM, CEM_INIT, SEM ...
	void recopy(Parameter * otherParameter);

	///create Scatter from "Binary Parameter Ekjh"
	virtual void createScatter(float *** scatter) = 0;
	virtual float *** scatterToArray() const = 0;

	void updateForCV(Model * originalModel, CVBlock & CVBlock);

protected:

	///compute TabCenter if there is only One Cluster
	// _tabCenter will not be changed
	void getTabCenterIfOneCluster(int * tabCenter, float * tabNbSampleInMajorModality, 
			float ** tabNbSamplePerModality = NULL) const;

	/// Table of centers vector of each cluster
	int ** _tabCenter;

	/// Table of modality
	int * _tabNbModality;

	/// Total number of modality
	int _totalNbModality;
};

/// compute Pdf in case nbCluster=1 (Scatter is a scalar)
float computePdfOneCluster(Sample * x, int * Center, float Scatter, int * tabNbModality);

/// compute Pdf in case nbCluster=1 (Scatter is a array of float, depends on variables)
float computePdfOneCluster(Sample * x, int * Center, 
		float * Scatter, int * tabNbModality);

/// compute Pdf in case nbCluster=1 
// (Scatter is a array of float*float, depends on variables and modalities)
float computePdfOneCluster(Sample * x, int * Center, 
		float ** Scatter, int * tabNbModality);

//---------------
// inline methods
//---------------

inline int ** BinaryParameter::getTabCenter() const {
	return _tabCenter;
}

inline int * BinaryParameter::getTabNbModality() const {
	return _tabNbModality;
}

inline int BinaryParameter::getTotalNbModality() const {
	return _totalNbModality;
}

}

#endif
