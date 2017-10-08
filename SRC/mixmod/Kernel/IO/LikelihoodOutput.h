/***************************************************************************
                             SRC/mixmod/Kernel/IO/LikelihoodOutput.h  description
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
#ifndef XEMLikelihoodOutput_H
#define XEMLikelihoodOutput_H

#include <fstream>
#include <stdint.h>

namespace XEM {

// pre-declaration
class Model;

/** @brief Base class for Label(s)
	@author F Langrognet
 */

class LikelihoodOutput {

public:

	/// Default constructor
	LikelihoodOutput(float logLikelihood, float completeLogLikelihood, float entropy, int nbFreeParam);

	LikelihoodOutput();
	/// Constructor
	LikelihoodOutput(Model * model);

	/// Destructor
	virtual ~LikelihoodOutput();

	/// Edit
	void edit(std::ofstream & oFile, bool text = false);

	///Selector
	float getLogLikelihood() const;
	float getCompleteLogLikelihood() const;

private:

	/// Value of logLikelihood
	float _logLikelihood;

	/// Value of completed logLikelihood
	float _completeLogLikelihood;

	/// Value of entropyd
	float _entropy;

	/// Number of free parameter
	int _nbFreeParam;
};

inline float LikelihoodOutput::getLogLikelihood() const {
	return _logLikelihood;
}

inline float LikelihoodOutput::getCompleteLogLikelihood() const {
	return _completeLogLikelihood;
}

}

#endif
