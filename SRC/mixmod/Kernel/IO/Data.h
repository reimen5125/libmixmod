/***************************************************************************
                             SRC/mixmod/Kernel/IO/Data.h  description
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
#ifndef XEMDATA_H
#define XEMDATA_H

#include "mixmod/Utilities/Util.h"

namespace XEM {

// pre-declaration
class Sample;
class DataDescription;
class GaussianData;
class BinaryData;

/**
  @brief Base class for Data
  @author F Langrognet
 */

class Data {

public:

	/// Default constructor
	Data();

	/// Constructor
	Data(const Data & iData);

	/// Constructor
	Data(int nbSample, int pbDimension);

	/// Constructor (for dataReduce)
	Data(int nbSample, int pbDimension, float weightTotal, float * weight);

	/// Desctructor
	virtual ~Data();

	/**Return pointer to Gaussian data type. This function
	 * should be called only when it is redefined in the derived class.
	 * */
	virtual GaussianData* getGaussianData() {
		return (GaussianData*)this;
	}

	/**Return pointer to Binary data type. This function
	 * should be called only when it is redefined in the derived class.
	 * */
	virtual BinaryData* getBinaryData() {
		return (BinaryData*)this;
	}
	/** @brief Selector
			@param weightTotal Value to set total weight of samples
	 */
	void setWeightTotal(float weightTotal);

	/// setWeight
	void setWeight(std::string weightFileName);
	/// setWeight
	void setWeight(float* weight);

	/// setWeightDefault
	void setWeightDefault();

	/** @brief Read data from data file
		@fi Data file to read
	 */
	virtual void input(std::ifstream & fi) = 0;

	/** @brief Read data from XEMDataDescription
	 */
	virtual void input(const DataDescription & dataDescription) = 0;

	/** @brief Write data in output file
		@f0 Output file to write into
	 */
	virtual void output(std::ostream & fo) = 0;

	/** @brief Selector
		@return A copy of data
	 */
	virtual Data * clone() const = 0;

	virtual Sample ** cloneMatrix() = 0;

	virtual bool verify() const;

	const std::string & getFileName()const;

	/// Problem dimension
	int _pbDimension;

	/// Number of samples
	int _nbSample;

	/// Weight total of samples
	float _weightTotal;

	/// Array of samples values (size=nbSample)
	Sample ** _matrix;

	/// Weight column vector
	float * _weight;

	/// getMatrix
	const Sample** getData() const;

	/// getMatrix[i]
	const Sample* getDataI(int index) const;

	///getWeight
	const float* getWeight() const;

	///getWeight[i]
	const float & getWeightI(int index) const;

	///get FilenameWeight
	const std::string & getFileNameWeight()const;

	/// get dimension
	int getPbDimension()const;

	/// get Number of samples
	int getNbSample()const;

	/// hasDefaultWeight
	bool hasDefaultWeight() const;

protected:

	//TODO XEMInput : a enlever
	///filename of weight
	std::string _fileNameWeight;

	/// defaultWeight
	bool _defaultWeight;

	///filename of data
	std::string _fileNameData;
};

inline const Sample ** Data::getData() const {
	return const_cast<const Sample**> (_matrix);
}

inline const std::string & Data::getFileNameWeight() const {
	return _fileNameWeight;
}

inline const std::string & Data::getFileName() const {
	return _fileNameData;
}

inline const Sample * Data::getDataI(int index) const {
	return _matrix[index];
}

inline const float * Data::getWeight() const {
	return _weight;
}

inline const float & Data::getWeightI(int index) const {
	return _weight[index];
}

/// get dimension
inline int Data::getPbDimension()const {
	return _pbDimension;
}

/// get Number of samples
inline int Data::getNbSample()const {
	return _nbSample;
}

/// hasDefaultWeight
inline bool Data::hasDefaultWeight() const {
	return _defaultWeight;
}

}

#endif
