/***************************************************************************
                             SRC/mixmod/Kernel/IO/Label.cpp  description
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
#include "mixmod/Kernel/IO/Label.h"
#include "mixmod/Kernel/IO/LabelDescription.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/BinaryModel.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/IO/IndividualColumnDescription.h"
#include <algorithm>

namespace XEM {

//------------
// Constructor
//------------
Label::Label() {
	_nbSample = 0;
}

//------------
// Constructor
//------------
Label::Label(int nbSample) {
	_nbSample = nbSample;
	_label.resize(_nbSample);
}

//------------
// Constructor
//------------
Label::Label(Model * model) {

	if (model == NULL) {
		THROW(OtherException, internalMixmodError);
	}

	// compute tabLabel
	//-----------------
	//int * tabLabel_p = NULL;
    std::unique_ptr<int[]>  tabLabel;
	int nbCluster = model->getNbCluster();
	ModelType * modelType = model->getModelType();
	bool binary = isBinary(modelType->_nameModel);

	if (!binary || (binary && !DATA_REDUCE)) {
		_nbSample = model->getNbSample();
		//int ** tabPartition = new int*[_nbSample];
        std::unique_ptr<int*[], TabDeleter<int>>  tabPartition(new int*[_nbSample], TabDeleter<int>(_nbSample));
		for (int i = 0; i < _nbSample; i++) {
			tabPartition[i] = new int[nbCluster];
		}
		//tabLabel = new int[_nbSample];
        tabLabel.reset(new int[_nbSample]); //provides exception-safe deletion
		model->getLabelAndPartitionByMAPOrKnownPartition(tabLabel.get(), tabPartition.get());
		//for (int i = 0; i < _nbSample; i++) {
		//	delete [] tabPartition[i];
		//}
		//delete [] tabPartition;
	}

	else {
		//binary case
		const vector<int> & correspondenceOriginDataToReduceData =
				dynamic_cast<BinaryModel*> (model)->getCorrespondenceOriginDataToReduceData();
		_nbSample = correspondenceOriginDataToReduceData.size();
		//tabLabel = new int[_nbSample];
        tabLabel.reset(new int[_nbSample]); //provides exception-safe deletion
		//label et partition on reduceData
		int nbSampleOfDataReduce = model->getNbSample();
		//int * tabLabelReduce = new int[nbSampleOfDataReduce];
		std::unique_ptr<int[]> tabLabelReduce(new int[nbSampleOfDataReduce]);        
		//int ** tabPartitionReduce = new int*[nbSampleOfDataReduce];
        std::unique_ptr<int*[], TabDeleter<int>>  tabPartitionReduce(new int*[nbSampleOfDataReduce], TabDeleter<int>(nbSampleOfDataReduce));
		for (int i = 0; i < nbSampleOfDataReduce; i++) {
			tabPartitionReduce[i] = new int[nbCluster];
		}
		model->getLabelAndPartitionByMAPOrKnownPartition(tabLabelReduce.get(), tabPartitionReduce.get());

		//  float ** tabPostProbaReduce = NULL;
		//  tabPostProbaReduce = copyTab(estimation->getModel()->getPostProba(),
		//  nbSampleOfDataReduce, nbCluster); // copy

		// convert labelReduce, partitionReduce, postProbaReduce to label, partition, postProba
		for (int i = 0; i < _nbSample; i++) {
			tabLabel[i] = tabLabelReduce[correspondenceOriginDataToReduceData[i]];
		}

		//delete //deletion is made by unique_ptr
		//for (int i = 0; i < nbSampleOfDataReduce; i++) {
		//	delete [] tabPartitionReduce[i];
		//}
		//delete [] tabPartitionReduce;
		//delete[] tabLabelReduce;
	}

	// compute _label
	recopyTabToVector(tabLabel.get(), _label, _nbSample);
	//delete [] tabLabel; //deletion done by unique_ptr
}

//------------
// Constructor
//------------
Label::Label(const Label & iLabel) {
	_nbSample = iLabel.getNbSample();
	_label = iLabel.getLabel();
}

//-----------
// Destructor
//-----------
Label::~Label() {
}

//--------------------
/// Comparison operator
//--------------------
bool Label::operator ==(const Label & label) const {
	if (_nbSample != label.getNbSample()) return false;
	for (int i = 0; i < _nbSample; i++) {
		if (_label[i] != label.getLabel()[i]) return false;
	}
	return true;
}

//----------
// editProba
//----------
void Label::edit(std::ostream & stream) const {
	stream.setf(ios::fixed, ios::floatfield);
	for (int i = 0; i < _nbSample; i++) {
		stream << _label[i] << endl;
	}
}

//---------
// getProba
//---------
int * Label::getTabLabel() const {
	int * res;
	recopyVectorToTab(_label, res);
	return res;
}

//---------
// get Error Rate
//---------
const float Label::getErrorRate(std::vector<int> const & label) const {
	if (_nbSample != (int) label.size()) {
		THROW(InputException, badNumberOfValuesInLabelInput);
	}

	float missClass = 0.0;
	for (int i = 0; i < _nbSample; i++) {
		if (_label[i] != label[i]) ++missClass;
	}
	return missClass / _nbSample;
}

//---------
// get getClassificationTab
//---------
int** Label::getClassificationTab(std::vector<int> const & label, int nbCluster) const {
	if (_nbSample != (int) label.size()) {
		THROW(InputException, badNumberOfValuesInLabelInput);
	}

	// memory allocation
	int** classTab = new int*[nbCluster];
	for (unsigned int i = 0; i < nbCluster; i++) {
		classTab[i] = new int[nbCluster];
	}
	// initialization
	for (unsigned int i = 0; i < nbCluster; i++)
		for (unsigned int j = 0; j < nbCluster; j++)
			classTab[i][j] = 0;

	// loop over labels
	for (int i = 0; i < _nbSample; i++) {
	  //cout<<_label[i]<<endl;
      if(label[i] > 0){
		++classTab[_label[i] - 1][label[i] - 1];
      }
	}

	return classTab;
}

// -----------
//input stream
// read labels between 1 and nbCluster
// -----------
void Label::input(std::ifstream & flux, int nbCluster) {
	int i = 0;
	int read;

	while (i < _nbSample && !flux.eof()) {
		flux >> read;
		if (read >= 1 && read <= nbCluster) {
			_label[i] = read;
		}
		else {
			THROW(InputException, badValueInLabelInput);
		}
		i++;
	}

	if (!flux.eof() && i != _nbSample) {
		THROW(InputException, badNumberOfValuesInLabelInput);
	}
}

// -----------
//input stream
// read labels between 1 and nbCluster
// -----------
void Label::input(const LabelDescription & labelDescription) {
	int i = 0;
	int readLabel;

	std::string labelFilename = labelDescription.getFileName();

	_nbSample = labelDescription.getNbSample();

	std::ifstream fi((labelFilename).c_str(), ios::in);
	if (!fi.is_open()) {
		THROW(InputException, wrongDataFileName);
	}

	while (i < _nbSample && !fi.eof()) {
		for (int j = 0; j < labelDescription.getNbColumn(); ++j) {
			if (fi.eof()) {
				THROW(InputException, endDataFileReach);
			}
			if (typeid (*(labelDescription.getColumnDescription(j)))
					== typeid (IndividualColumnDescription))
			{
				std::string stringTmp;
				fi >> stringTmp;
				//cout<<stringTmp<<endl;
			}
			else {
				fi >> readLabel;
				_label.push_back(readLabel);
			}
		}
		++i;
	}
	if (!fi.eof() && i != _nbSample) {
		THROW(InputException, badNumberOfValuesInLabelInput);
	}
}

}
