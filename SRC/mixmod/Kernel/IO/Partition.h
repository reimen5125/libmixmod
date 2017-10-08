/***************************************************************************
                             SRC/mixmod/Kernel/IO/Partition.h  description
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
#ifndef XEMPARTITION_H
#define XEMPARTITION_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

// pre-declaration
class Label;

/** @brief Base class for Partition(s)
@author F Langrognet
 */

class Partition {

public:

	/// Default constructor
	Partition();

	/// Constructor
	Partition(Partition * iPartition);

	/// Constructor 
	Partition(int nbSample, int nbCluster, const NumericPartitionFile & partitionFile);

	/// constructor from a XEMLabel
	Partition(const Label * label, int nbCluster);

	/// Constructor
	// used in DCV context
	Partition(Partition * originalPartition, CVBlock & block);

	/// Destructor
	virtual ~Partition();

	/** @brief Set the dimension of the Partition
	@param nbSample The number of samples
	@param nbCluster The number of clusters
	 */
	void setDimension(int nbSample, int nbCluster);

	int getGroupNumber(int idxSample);

	/** verify if partition is complete 
	- each line has one (and only one) '1'
	- each cluster appears at least one time
	 */
	bool isComplete();
	int** getTabValue() const;
	int* getTabValueI(int index) const;

	///get Format
	const NumericPartitionFile & getPartitionFile() const;
    void setPartitionFile(std::string f, TypePartition::TypePartition type);
	/// get Number of samples
	int getNbSample()const;

	/// get Number of clusters
	int getNbCluster() const;

	bool operator==(Partition & otherPartition);

	/** @brief Friend method
	@return Operator >> overloaded to read Partition from input files
	 */
	friend std::ifstream & operator >>(std::ifstream & fi, Partition & partition);

	/** @brief Friend method
	@return Operator << overloaded to write Partition 
	 */
	friend std::ostream & operator <<(std::ostream & fo, const Partition & partition);

//TODO
//private:
	
	/// Number of samples
	int _nbSample;

	/// Number of clusters
	int _nbCluster;

	int ** _tabValue; // aggregate

	///name of partition
	NumericPartitionFile _partitionFile;

	bool _deleteValues;
};

inline int ** Partition::getTabValue() const {
	return _tabValue;
}

inline int * Partition::getTabValueI(int index) const {
	return _tabValue[index];
}

inline const NumericPartitionFile & Partition::getPartitionFile() const {
	return _partitionFile;
}

inline int Partition::getNbSample()const {
	return _nbSample;
}

inline int Partition::getNbCluster() const {
	return _nbCluster;
}

}

#endif
