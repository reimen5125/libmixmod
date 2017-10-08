/***************************************************************************
                             SRC/mixmod/Kernel/Model/ModelType.h  description
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
#ifndef XEMMODELTYPE_H
#define XEMMODELTYPE_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

/**
  @brief Base class for ModelType(s)
  @author F Langrognet
 */

class ModelType {

public:

	/// Default constructor
	ModelType();

	// constructor
	ModelType(ModelName name, int nbSubDimensionFree = 0);

	// copy constructor
	ModelType(const ModelType & iModelType);

	/// Destructor
	~ModelType();

	/// Comparison operator
	bool operator ==(const ModelType & modelType) const;

	/// Input model type
	void input(std::ifstream & fi, int nbCluster);

	/// name of the model
	ModelName _nameModel;

	ModelType* clone();

	//// list of number of subDimensionEqual
	//int _nbSubDimensionEqual;
	//// list of number of subDimensionFree
	//int _nbSubDimensionFree;

	/// list of subDimensionEqual
	int _subDimensionEqual;

	/// _nbSubDimensionFree : size of array _tabSubDimensionFree
	int _nbSubDimensionFree;

	/// array of subDimensionFree
	int * _tabSubDimensionFree;

	/// getModelName
	const ModelName & getModelName() const;

	/// getSubDimensionEqual
	const int & getSubDimensionEqual() const;

	/// getTabSubDimensionFree
	const int * getTabSubDimensionFree() const;

	///getTabSubDimensionFreeI
	const int & getTabSubDimensionFreeI(int index) const;

	/// setSubDimensionFree
	void setTabSubDimensionFree(int iTabSubDimensionFree, int position);

	/// setSubDimensionEqual
	void setSubDimensionEqual(int iSubDimensionEqual);

	/// <<
	friend std::ostream & operator<<(std::ostream & fo, ModelType & modelType);

	// print out model type
	void print(std::ostream & flux) const;
	// print out model type short cut
	void printShortcut(std::ostream & flux) const;
	/// editModelType
	void edit(std::ostream & oFile);
};

inline const ModelName & ModelType::getModelName() const {
	return _nameModel;
}

inline const int & ModelType::getSubDimensionEqual() const {
	return _subDimensionEqual;
}

inline const int * ModelType::getTabSubDimensionFree() const {
	return _tabSubDimensionFree;
}

inline const int & ModelType::getTabSubDimensionFreeI(int index) const {
	return _tabSubDimensionFree[index];
}


}

#endif
