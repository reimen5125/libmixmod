/***************************************************************************
                             SRC/mixmod/Kernel/IO/GaussianSample.h  description
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
#ifndef XEMGAUSSIANSample_H
#define XEMGAUSSIANSample_H

#include "mixmod/Kernel/IO/Sample.h"

namespace XEM {

/**
  @brief Base class for Sample
  @author F Langrognet
 */

class GaussianSample : public Sample {

public:

	/// Constructor
	GaussianSample();

	/// Constructor
	GaussianSample(int pbDimension);

	/// Constructor
	GaussianSample(GaussianSample * iSample);

	/// Constructor
	GaussianSample(int pbDimension, float * tabValue);

	/// Destructor
	virtual ~GaussianSample();

	/// Set value vector of sample
	void setDataTabValue(float * tabValue);

	/// Set one value of sample
	void setDataValue(int idxDim, float value);

	/// get value vector of sample
	float * getTabValue() const;

	/// get one value of sample
	float getDataValue(int idxDim) const;

protected:

	/// Vector of sample value
	float * _value;
};

//---------------
// inline methods
//---------------

inline float * GaussianSample::getTabValue() const {
	return _value;
}

inline float GaussianSample::getDataValue(int idxDim) const {
	return _value[idxDim];
}

}

#endif
