/***************************************************************************
                             SRC/mixmod/Kernel/IO/GaussianData.h  description
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
#ifndef XEMGAUSSIANDATA_H
#define XEMGAUSSIANDATA_H

#include "mixmod/Kernel/IO/Data.h"

namespace XEM {

/**
  @brief Base class for Gaussian Data
  @author F Langrognet
 */

class GaussianData : public Data {

public:

	/// Default constructor
	GaussianData();

	/// Constructor
	GaussianData(const GaussianData & iData);

	/// Constructor
	GaussianData(int nbSample, int pbDimension, const std::string & dataFileName);

	/// Constructor  (without fill matrix)
	GaussianData(int nbSample, int pbDimension);

	/// Constructor  (with matrix)
	GaussianData(int nbSample, int pbDimension, float ** matrix);

	/// Constructor
	GaussianData(int nbSample, int pbDimension, float weightTotal, Sample **& matrix, float * weight);

	/// Constructor
	// used in DCV context
	GaussianData(int nbSample, int pbDimension, Data * originalData, CVBlock & block);

	/// Destructor
	virtual ~GaussianData();

	/** @brief Selector
		@return A copy of data
	 */
	virtual Data * clone() const;

	/**  @brief Copy
		 @return A copy data matrix
	 */
	virtual Sample ** cloneMatrix();

	//TODO a enlever XEMInput
	/** @brief Read data from gaussian data file
		@fi Gaussian Data file to read
	 */
	virtual void input(std::ifstream & fi);

	/** @brief Read data from XEMDataDescription
	 */
	virtual void input(const DataDescription & dataDescription);

	/** @brief Write gaussian data in output file
		@f0 Output file to write into
	 */
	virtual void output(std::ostream & fo);

	virtual bool verify()const;

	/// pointer to stored values
	float ** _yStore;

	float ** getYStore();

	float getInv2PiPow() const;

	float getHalfPbDimensionLog2Pi() const;

	float getPbDimensionLog2Pi() const;

	float * getTmpTabOfSizePbDimension() const;

protected:

	/// 1/ (2 * pi)^(d/2)
	float _Inv2PiPow;

	/// 0.5 * p * logf(2 * PI)
	float _halfPbDimensionLog2Pi;

	float _pbDimensionLog2Pi;

	// tableau de float longueur _pbDimension 
	// utilisé pour ne pas avoir à faire sans arret allocation / destruction
	// utilisé par XEMnorme, getLogLikeLihoodOne
	float * __tmpTabOfSizePbDimension;

	bool _deleteSamples;
};

inline float ** GaussianData::getYStore() {
	return _yStore;
}

inline float GaussianData::getInv2PiPow() const {
	return _Inv2PiPow;
}

inline float GaussianData::getHalfPbDimensionLog2Pi()const {
	return _halfPbDimensionLog2Pi;
}

inline float GaussianData::getPbDimensionLog2Pi()const {
	return _pbDimensionLog2Pi;
}

inline float* GaussianData::getTmpTabOfSizePbDimension()const {
	return __tmpTabOfSizePbDimension;
}

}

#endif
