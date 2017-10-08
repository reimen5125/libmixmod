/***************************************************************************
                             SRC/mixmod/Kernel/IO/BinarySample.h  description
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
#ifndef XEMBINARYSample_H
#define XEMBINARYSample_H

#include "mixmod/Kernel/IO/Sample.h"

namespace XEM {

/**
  @brief Base class for Sample
  @author F Langrognet 
 */

class BinarySample : public Sample {

public:

	/// Constructor
	BinarySample();

	/// Constructor
	BinarySample(int pbDimension);

	/// Constructor
	BinarySample(BinarySample * iSample);

	/// Constructor
	BinarySample(int pbDimension, int * tabValue);

	/// Destructor
	virtual ~BinarySample();

	/// Set value vector of sample
	void setDataTabValue(int * tabValue);

	/// Set one value of sample
	void setDataValue(int idxDim, int iValue);

	/// get value vector of sample
	int * getTabValue() const;

	/// get one value of sample
	int getDataValue(int idxDim) const;


protected:

	/// Vector of sample value
	int * _value;
};

//---------------
// inline methods
//---------------

inline int * BinarySample::getTabValue() const {
	return _value;
}

inline int BinarySample::getDataValue(int idxDim) const {
	return _value[idxDim];
}

}

#endif
