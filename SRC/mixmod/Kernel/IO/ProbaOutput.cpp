/***************************************************************************
                             SRC/mixmod/Kernel/IO/ProbaOutput.cpp  description
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

#include "mixmod/Kernel/IO/ProbaOutput.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/BinaryModel.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include <vector>

namespace XEM {

//------------
// Constructor
//------------
ProbaOutput::ProbaOutput() {
	_CVLabelAvailable = false;
	_tabLabel = NULL;
	_tabCVLabel = NULL;
	_tabPartition = NULL;
	_tabPostProba = NULL;
}

//------------
// Constructor
//------------
ProbaOutput::ProbaOutput(Model * model) {

	_CVLabelAvailable = false;
	_tabCVLabel = NULL;
	_nbCluster = model->getNbCluster();

	if (model == NULL) {
		THROW(OtherException, internalMixmodError);
	}

	bool binary = isBinary(model->getModelType()->_nameModel);
	if (!binary || (binary && !DATA_REDUCE)) {
		// gaussian case 
		_nbSample = model->getNbSample();
		_tabLabel = new int[_nbSample];
		_tabPartition = new int *[_nbSample];
		for (int i = 0; i < _nbSample; i++) {
			_tabPartition[i] = new int[_nbCluster];
		}
		model->getLabelAndPartitionByMAPOrKnownPartition(_tabLabel, _tabPartition);
		// copy
		_tabPostProba = copyTab(model->getPostProba(), _nbSample, _nbCluster);
	}
	else {
		const std::vector<int> & correspondenceOriginDataToReduceData = 
				dynamic_cast<BinaryModel*> (model)->getCorrespondenceOriginDataToReduceData();
		//binary case
		_nbSample = correspondenceOriginDataToReduceData.size();
		_tabLabel = new int[_nbSample];
		_tabPartition = new int *[_nbSample];
		for (int i = 0; i < _nbSample; i++) {
			_tabPartition[i] = new int[_nbCluster];
		}
		_tabPostProba = new float *[_nbSample];
		for (int i = 0; i < _nbSample; i++) {
			_tabPostProba[i] = new float[_nbCluster];
		}

		//label et partition on reduceData
		int nbSampleOfDataReduce = model->getNbSample();
		//int ** tabPartitionReduce = new int*[nbSampleOfDataReduce];
        std::unique_ptr<int*[], TabDeleter<int>>  tabPartitionReduce(new int*[nbSampleOfDataReduce],
                                                                             TabDeleter<int>(nbSampleOfDataReduce));
		for (int i = 0; i < nbSampleOfDataReduce; i++) {
			tabPartitionReduce[i] = new int[_nbCluster];
		}
		//int * tabLabelReduce = new int[nbSampleOfDataReduce];
        std::unique_ptr<int[]> tabLabelReduce(new int[nbSampleOfDataReduce]);
		model->getLabelAndPartitionByMAPOrKnownPartition(tabLabelReduce.get(), tabPartitionReduce.get());

		//float ** tabPostProbaReduce = NULL;
		// copy
		//tabPostProbaReduce = copyTab(model->getPostProba(), nbSampleOfDataReduce, _nbCluster);
        std::unique_ptr<float*[], TabDeleter<float>>  tabPostProbaReduce(copyTab(model->getPostProba(), nbSampleOfDataReduce, _nbCluster),
                                                                           TabDeleter<float>(nbSampleOfDataReduce));
		// convert labelReduce, partitionReduce, postProbaReduce to label, partition, postProba
		for (int i = 0; i < _nbSample; i++) {
			_tabLabel[i] = tabLabelReduce[correspondenceOriginDataToReduceData[i]];
			for (int k = 0; k < _nbCluster; k++) {
				int index = correspondenceOriginDataToReduceData[i];
				_tabPostProba[i][k] = tabPostProbaReduce[index][k];
				_tabPartition[i][k] = tabPartitionReduce[index][k];
			}
		}

		//delete
		//for (int i = 0; i < nbSampleOfDataReduce; i++) {
		//	delete [] tabPartitionReduce[i];
		//}
		//delete [] tabPartitionReduce;

		//for (int i = 0; i < nbSampleOfDataReduce; i++) {
		//	delete [] tabPostProbaReduce[i];
		//}
		//delete [] tabPostProbaReduce;

		//delete[] tabLabelReduce;
	}
}

ProbaOutput::ProbaOutput(ProbaOutput * iProbaOutput) {
	_nbSample = iProbaOutput->getNbSample();
	_nbCluster = iProbaOutput->getNbCluster();

	_tabPostProba = copyTab(iProbaOutput->getTabPostProba(), _nbSample, _nbCluster);
	_tabLabel = copyTab(iProbaOutput->getTabLabel(), _nbSample);

	_tabPartition = NULL;
	_tabCVLabel = NULL;
	_CVLabelAvailable = false;
}

//-----------
// Destructor
//-----------
ProbaOutput::~ProbaOutput() {
	if (_tabLabel) {
		delete[] _tabLabel;
		_tabLabel = NULL;
	}

	if (_tabCVLabel) {
		delete[] _tabCVLabel;
		_tabCVLabel = NULL;
	}

	int i;
	if (_tabPartition) {
		for (i = 0; i < _nbSample; i++) {
			delete[] _tabPartition[i];
			_tabPartition[i] = NULL;
		}
		delete [] _tabPartition;
		_tabPartition = NULL;
	}

	if (_tabPostProba) {
		for (i = 0; i < _nbSample; i++) {
			delete[] _tabPostProba[i];
			_tabPostProba[i] = NULL;
		}
		delete [] _tabPostProba;
		_tabPostProba = NULL;
	}
}

//-----------
// setCVLabel
//-----------
void ProbaOutput::setCVLabel(int * CVLabel) {
	_CVLabelAvailable = true;
	_tabCVLabel = new int[_nbSample];
	recopyTab(CVLabel, _tabCVLabel, _nbSample);
}

//--------------
// editPartition
//--------------
void ProbaOutput::editPartition(std::ofstream & oFile) {
	int k, i;

	int ** p_tabPartition;
	int * p_tabPartition_i;

	p_tabPartition = _tabPartition;
	for (i = 0; i < _nbSample; i++) {
		p_tabPartition_i = *p_tabPartition;
		for (k = 0; k < _nbCluster; k++) {
			oFile << p_tabPartition_i[k] << "\t";
		}
		oFile << endl;

		p_tabPartition++;
	}
}

//----------
// editLabel
//----------
void ProbaOutput::editLabel(std::ofstream & oFile) {
	int i;
	for (i = 0; i < _nbSample; i++)
		oFile << _tabLabel[i] << endl;

}

//----------
// editLabel (verbose: to standard output)
//----------
void ProbaOutput::editLabel() {
	int i;
	for (i = 0; i < _nbSample; i++)
		cout << _tabLabel[i] << endl;
}

//--------------
// editPostProba
//--------------
void ProbaOutput::editPostProba(std::ofstream & oFile) {
	oFile.setf(ios::fixed, ios::floatfield);
	editTab(_tabPostProba, _nbSample, _nbCluster, oFile, "\t", "");
}

//------------
// editCVLabel
//------------
void ProbaOutput::editCVLabel(std::ofstream & oFile) {
	if (_CVLabelAvailable) {
		int i;
		for (i = 0; i < _nbSample; i++)
			oFile << _tabCVLabel[i] << endl;
	}
}

ProbaOutput * ProbaOutput::clone() {
	ProbaOutput * newProbaOutput = new ProbaOutput(this);
	return newProbaOutput; // ceci n'etait pas fait...
}

}
