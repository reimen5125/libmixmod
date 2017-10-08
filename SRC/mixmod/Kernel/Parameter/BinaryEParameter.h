/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/BinaryEParameter.h  description
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
#ifndef XEMBINARYEPARAMETER_H
#define XEMBINARYEPARAMETER_H

#include "mixmod/Kernel/Parameter/BinaryParameter.h"

namespace XEM {

/**
@author F LANGROGNET
 */
class BinaryEParameter : public BinaryParameter {

public:
	
	/// Default constructor
	BinaryEParameter();

	/// Constructor
	// called by XEMModel
	BinaryEParameter(Model * iModel, ModelType * iModelType, int * tabNbModality);

	/// Constructor
	BinaryEParameter(const BinaryEParameter * iParameter);

	/// Destructor
	~BinaryEParameter();

	/// Comparison operator
	virtual bool operator ==(const BinaryEParameter & param) const;

	/// reset to default values
	virtual void reset();

	/// clone
	Parameter * clone() const;

	/// selector :  return scatter value
	float getScatter() const;

	/// getFreeParameter
	int getFreeParameter() const;

	float getPdf(int iSample, int kCluster) const;

	float getLogPdf(int iSample, int kCluster) const;

	/** Compute probability density function
		 for x vector and kCluster th cluster
	 */
	float getPdf(Sample * x, int kCluster) const;

	/// getlogLikelihoodOne (one cluster)
	float getLogLikelihoodOne() const;

	/// Compute scatter(s) 
	void computeScatter();

	/// Compute random scatter(s) 
	void computeRandomScatter();

	///recopy scatter from param (used for init  : USER)
	void recopyScatter(Parameter * iParam);

	///create Scatter from "Binary Parameter Ekjh"
	void createScatter(float *** scatter);

	/// editScatter (for debug)
	void editScatter(int k);

	/// editScatter 
	void editScatter(std::ofstream & oFile, int k, bool text = false);

	// Read Scatter in input file
	void inputScatter(std::ifstream & fi, int k);
	void inputScatter(float *** scatters);

	float *** scatterToArray() const;

private:
	
	/// scatter
	float _scatter;
};

inline float BinaryEParameter::getScatter() const {
	return _scatter;
}

}

#endif
