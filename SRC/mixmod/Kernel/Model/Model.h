/***************************************************************************
                             SRC/mixmod/Kernel/Model/Model.h  description
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
#ifndef XEMMODEL_H
#define XEMMODEL_H

#include "mixmod/Utilities/Util.h"
#include "mixmod/Utilities/Error.h"
#include "mixmod/Kernel/IO/Data.h"
#include "mixmod/Kernel/Parameter/Parameter.h"

namespace XEM {

/**
	  @brief Base class for Model(s)
	  @author F Langrognet
 */

// pre-declaration
class Parameter;
class GaussianData;
class BinaryData;
class ClusteringStrategyInit;
class Partition;
class Sample;
class ModelType;
class LabelDescription;

class Model {

public:

	/// Default constructor
	Model();

	//clone the model
	virtual Model * clone();
	/// Constructor
	Model(Model * iModel);

	/// Constructor
	Model(ModelType * modelType, int nbCluster, Data *& data, Partition * knownPartition);

	/// Destructor
	virtual ~Model();

	void updateForCV(Model * originalModel, CVBlock & CVBlock);
	
	
	//-------
	// select
	//-------

	/** @brief Selector
		@return The parameters
	 */
	Parameter * getParameter();
	GaussianParameter * getGaussianParameter();
	BinaryParameter * getBinaryParameter();

	/** @brief Selector
		@return The current number for cluster
	 */
	int getNbCluster();

	/** @brief Selector
		@return The current data
	 */
	Data * getData();

	/**
	 * @brief Return Gaussian data
	 * @return
	 */
	GaussianData* getGaussianData();
	/**
	 * @brief Return Binary data
	 * @return
	 */
	BinaryData* getBinaryData();

	/** @brief Selector
		@return The type of Error
	 */
	Exception& getErrorType() const;

	/** @brief Selector
	 @return The type of the model
	 */
	ModelType * const & getModelType() const;

	/** @brief Selector
		@return The number of samples
	 */
	int getNbSample();

	/** @brief Selector
		@return Table of Fik of each cluster : probabilitites: _fik = pk * f(xi,muk,Sk)
	 */
	float ** getTabFik();

	/// return _tabSumF
	float * getTabSumF();

	/** @brief Selector
	    @return Table of Tik of each cluster :
	            conditional probabilities that xi arises from the k-th 
	            mixture component, 0 <= tik[i]k0] <= 1
	 */
	float ** getTabTik();

	/** @brief Selector
		@return Table of Zik zik[i][k0] = 1 if xi arises from the k0-th mixture component, 0 else
	 */
	int ** getTabZikKnown();

	float ** getTabCik();

	/// getTabZikKnown
	bool * getTabZiKnown();

	/** @brief Selector
		@return Table of number of elements in each cluster
	 */
	float * getTabNk();

	bool getDeleteData();

	
	//---------
	// compute
	//--------

	/// compute _fik
	void computeFik();

	/// Compute the number of points in each class
	void computeNk();

	/** @brief Compute the log-likelihood
		@return The log-likelihood
	 */
	float getLogLikelihood(bool fikMustBeComputed);


	/** @brief Compute the log-likelihood with one cluster
		@return The log-likelihood
	 */
	float getLogLikelihoodOne();

	/** @brief Compute the entropy
		@return The entropy
	 */
	float getEntropy();

	/** @brief Compute the entropy matrix (for massiccc's visualization)
		@return The entropy matrix
	 */
	vector< vector<float> > getEntropyMatrix();

	/** @brief Compute the completed log-likelihood
		@return The completed log-likelihood
	 */
	float getCompletedLogLikelihood();

	/** get completed LL (if CEM) or LL (elseif)*/
	float getCompletedLogLikelihoodOrLogLikelihood();

	/// return the number of free parameters
	int getFreeParameter();

	/** @brief Selector
		@return Log of the weight total
	 */
	float getLogN();

	/// getLabel and partition
	/// label=1...nbSample
	void getLabelAndPartitionByMAPOrKnownPartition(int * label, int ** partition);

	/// get label of the ith individual (i=0 .... nbSample-1) by MAP (or known label)
	/// return value in [0 nbCluster-1]
	int getLabelByMAPOrKnownPartition(int i);

	/// get knownLabel of the ith individual (i=0 .... nbSample-1)
	/// return value in [0 nbCluster-1]
	/// throw an error if the label is unknown
	int getKnownLabel(int i);

	/// getPostProba
	float ** getPostProba();


	//--------
	// compute
	//--------

	/** @brief Compute the label of the i0-th point of the sample
		@return The label of i0 (i0=0 -> _nBSample -1)
	 */
	int computeLabel(int i0);

	/** @brief Compute the label of new point x
		@return The label of x
	 */
	int computeLabel(Sample * x);


	//------
	// algo
	//------

	/// Maximum a posteriori step method
	void MAPstep();

	/// Expectation step method
	void Estep();

	/// Maximization step method
	void Mstep();

	/// Stochastic classification step method
	void Sstep();

	/// Classification step method
	void Cstep();


	//-----
	// init
	//-----

	/// Random center initialization of the parameters of the model
	void initRANDOM(int nbTry);

	/// random step for init RANDOM or USER_PARTITION
	void randomForInitRANDOMorUSER_PARTITION(
			bool * tabIndividualCanBeUsedForInitRandom, bool * tabClusterToInitialize);

	/// User initialization of the parameters of the model
	void initUSER(Parameter * initParameter);

	/// User partition initialization of the parameters of the model
	void initUSER_PARTITION(Partition * initPartition, int nbTryInInit = defaultNbTryInInit);

	// set name of the algorithm
	void setParameter(Parameter * parameter);

	// set name of the algorithm
	void setAlgoName(AlgoName algoName);

	AlgoName getAlgoName();

	// set an error for the model
	void setError(Exception& errorType);

	/// Fix label Known
	void FixKnownPartition(Partition *& y);

	// edit debug information
	void editDebugInformation();
	void editFik();
	void editCik();
	void editTik();
	void editNk();

protected:

	/// type of the model
	ModelType * _modelType;

	/// Number of clusters
	int _nbCluster;

	/// Number of samples
	int _nbSample;

	/// Current data
	Data * _data;
	bool _deleteData;

	/// parameter of model
	Parameter * _parameter;

	/// Probabilities: _fik = pk * f(xi,muk,Sk)
	/// dim : _nbSample * _nbCluster
	float ** _tabFik;

	/// table of sum of _tabFik for all k (dim : _nbSample)
	float * _tabSumF;

	/// Conditional probabilities that x(i) arises from the k-th mixture component, 
	/// 0 <= tik[i][k0] <= 1. dim : _nbSample * _nbCluster
	float ** _tabTik;

	/// zikKnown : _tabZikKonwn[i][k0] = 1 if xi arises from the k0-th mixture component, 0 else.
	/// dim : _nbSample * _nbCluster
	int ** _tabZikKnown;

	/** classification array for individual i and class k
	  // if zikKnown 
	  //		cik = zikKnown
	  //	else :
	  //		cik = tik if EM
	  //		cik = zik by MAP rule if CEM or MAP
	  //		cik = 'random' if SEM
	 */
	float ** _tabCik;


	/// is the label zik known (fixed)
	bool * _tabZiKnown;


	/// Number of points in each class
	float * _tabNk;

	// name of the algorithm
	AlgoName _algoName;

	// Error handler
	Error _error;
};

//--------------
//inline methods
//--------------

inline bool * Model::getTabZiKnown() {
	return _tabZiKnown;
}

inline int ** Model::getTabZikKnown() {
	return _tabZikKnown;
}

inline float ** Model::getTabCik() {
	return _tabCik;
}

inline float ** Model::getTabTik() {
	return _tabTik;
}

inline float ** Model::getTabFik() {
	return _tabFik;
}

inline float * Model::getTabSumF() {
	return _tabSumF;
}

inline float * Model::getTabNk() {
	return _tabNk;
}

inline int Model::getNbCluster() {
	return _nbCluster;
}

inline Data * Model::getData() {
	return _data;
}

inline GaussianData* Model::getGaussianData() {
	return _data->getGaussianData();
}

inline BinaryData* Model::getBinaryData() {
	return _data->getBinaryData();
}

inline Parameter * Model::getParameter() {
	return _parameter;
}

inline GaussianParameter * Model::getGaussianParameter() {
	return _parameter->getGaussianParameter();
}

inline BinaryParameter * Model::getBinaryParameter() {
	return _parameter->getBinaryParameter();
}

inline int Model::getNbSample() {
	return _nbSample;
}

inline float ** Model::getPostProba() {
	return _tabTik;
}

inline ModelType * const & Model::getModelType() const {
	return _modelType;
}

inline Exception& Model::getErrorType() const {
	return _error.getError();
}

}

#endif
