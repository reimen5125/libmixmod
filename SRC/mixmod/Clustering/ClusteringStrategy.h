/***************************************************************************
                             SRC/mixmod/Clustering/ClusteringStrategy.h  description
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
#ifndef XEMClusteringStrategy_H
#define XEMClusteringStrategy_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

// pre-declaration
class ClusteringStrategyInit;
class ModelType;
class Parameter;
class Data;
class Model;
class Partition;

/**
 \class XEMClusteringOutput
 Main class for Mixmod Strategy
 @author F. Langrognet
		@date 2012
		@brief XEMClusteringStrategy class
 */
class ClusteringStrategy {

public:

	/// Default constructor
	ClusteringStrategy();

	/// Constructor
	ClusteringStrategy(const ClusteringStrategy & strategy);

	/// Destructor
	~ClusteringStrategy();

	///--------------------------------
	/// Strategy initialisation methods
	///--------------------------------

	///getStrategyInit
	const ClusteringStrategyInit * getStrategyInit() const;

	/// setStrategyInit
	void setStrategyInit(ClusteringStrategyInit * iStrategyInit);

	/// setStrategyInit
	void setStrategyInit(StrategyInitName  initName, Data *& data, int nbNbCluster,
			int * tabNbCluster, ModelType * modelType);

	/// setStrategyInitName
	void setStrategyInitName(StrategyInitName initName);

	/// setInitParam
	void setInitParam(std::string & paramFileName, int position);

	/// setInitParam
	void setTabInitParameter(Parameter ** tabInitParameter, int nbInitParameter);

	/// setInitPartition
	void setInitPartition(std::string & partitionFileName, int position);

	/// setInitPartition
	void setInitPartition(Partition * part, int position);

	/// setTabPartition
	void setTabPartition(Partition ** tabPartition, int nbPartition);

	/// getNbTryInInit
	const int getNbTryInInit() const;

	/// setNbTryInInit
	void setNbTryInInit(int nbTry);

	/// getNbIterationInInit
	const int getNbIterationInInit() const;

	/// set NbIterationInInit
	void setNbIterationInInit(int nbIteration);

	/// getEpsilonInInit
	const float getEpsilonInInit() const;

	/// setEpsilonInInit
	void setEpsilonInInit(float epsilon);

	/// getStopNameInInit
	const AlgoStopName getStopNameInInit() const;

	/// setStopNameInInit
	void setStopNameInInit(AlgoStopName stopName);

	///--------------------------------
	/// Algo methods
	///--------------------------------

	///getAlgo[i]
	const Algo * getAlgo(int index) const;

	/// setAlgo
	void setAlgo(AlgoName algoName, int position);

	/// addAlgo (and the end of the list)
	void addAlgo(AlgoName algoName);

	/// removeAlgo
	void removeAlgo(unsigned int  position);

	/// getTabAlgo
	std::vector<Algo*> const & getTabAlgo() const;
	std::vector<Algo*> & getTabAlgo();

	/// insertAlgo
	void insertAlgo(AlgoName algoName, int position);

	/// setAlgoStopRuleTypeValue
	void setAlgoStopRule(AlgoStopName stopName, int position);
	void setAlgoIteration(  int position, int nbIterationValue);
	void setAlgoEpsilon( int position, float epsilonValue);

	///nbTry
	const int getNbTry()const;
	void setNbTry(int nbTry);

	const int getNbAlgo() const;

	/// Input strategy (FLAT FORMAT)
	// TODO XEMInput : a enlever
	void input_FLAT_FORMAT(std::ifstream & fi, Data *& data, int nbNbCluster,
			int * tabNbCluster, ModelType * modelType);

	/// Run method
	void run(Model *& model) const;

	bool verify();

	//void edit(std::ofstream & oFile);
	void edit(std::ostream & out);

	friend std::ostream & operator << (std::ostream & fo, ClusteringStrategy & strategy);

private:

	/// Number of try in the strategy
	int _nbTry;

	/// strategyInit
	ClusteringStrategyInit * _strategyInit;

	/// Number of algorithm in the strategy
	int _nbAlgo;

	/// Table of algorithm
	std::vector<Algo*> _tabAlgo;    // aggregate

  /// One try (init+algo(s))
	void oneTry(Model *& model, bool doThrow=false) const;
};

inline const ClusteringStrategyInit * ClusteringStrategy::getStrategyInit() const {
	return _strategyInit;
}

inline const Algo * ClusteringStrategy::getAlgo(int index) const {
	return _tabAlgo[index];
}

inline std::vector<Algo*> const & ClusteringStrategy::getTabAlgo() const {
	return _tabAlgo;
}

inline std::vector<Algo*> & ClusteringStrategy::getTabAlgo() {
	return _tabAlgo;
}

inline const int ClusteringStrategy::getNbAlgo() const {
	return _nbAlgo;
}

inline const int ClusteringStrategy::getNbTry()const {
	return _nbTry;
}

}

#endif
