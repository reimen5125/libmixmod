/***************************************************************************
                             SRC/mixmod/Kernel/IO/Label.h  description
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
#ifndef XEMLabel_H
#define XEMLabel_H

/** @brief Base class for Label(s)
	@author F Langrognet & A Echenim
 */

#include <iostream>
#include <vector>
#include <stdint.h>

namespace XEM {

// pre-declaration
class Model;
class LabelDescription;

/**
 \class XEMLabel
 @author F. Langrognet
		@date 2010
		@brief XEMLabel class
 */
class Label {

public:

	/// Default constructor
	Label();

	/// Constructor
	Label(int nbSample);

	/// Constructor
	Label(Model * model);

	Label(const Label & iLabel);

	/// Destructor
	virtual ~Label();

	/// Comparison operator
	bool operator ==(const Label & label) const;

	/// edit labels
	void edit(std::ostream & stream) const;

	/// get label
	int * getTabLabel() const;

	/// get label
	std::vector<int> const & getLabel() const;

	/// set Label
	void setLabel(int * label, int nbSample);
	void setLabel(std::vector<int> label, int nbSample);

	/// Selector
	int getNbSample() const;

	// get Error Rate
	const float getErrorRate(std::vector<int> const & label) const;

	// get getClassificationTab
	int** getClassificationTab(std::vector<int> const & label, int nbCluster) const;

	///input stream
	void input(std::ifstream & flux, int nbCluster);
	void input(const LabelDescription & labelDescription);

private:

	/// Number of samples
	int _nbSample;

	std::vector<int> _label;

};

inline std::vector<int> const & Label::getLabel() const {
	return _label;
}

inline int Label::getNbSample() const {
	return _nbSample;
}

inline void Label::setLabel(std::vector<int> label, int nbSample) {
	_nbSample = nbSample;
	_label = label;
}

}

#endif
