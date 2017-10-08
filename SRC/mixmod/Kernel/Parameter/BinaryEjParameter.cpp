/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/BinaryEjParameter.cpp  description
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
#include "mixmod/Kernel/Parameter/BinaryEjParameter.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/BinarySample.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Utilities/Random.h"

namespace XEM {

//--------------------
// Default constructor
//--------------------
BinaryEjParameter::BinaryEjParameter() {
	THROW(OtherException, wrongConstructorType);
}

//-------------------------------
// Constructor called by XEMModel
//-------------------------------
BinaryEjParameter::BinaryEjParameter(
		Model * iModel, ModelType * iModelType, int * tabNbModality) 
: BinaryParameter(iModel, iModelType, tabNbModality) 
{
	_scatter = new float[_pbDimension];
	for (int j = 0; j < _pbDimension; j++) {
		_scatter[j] = 0;
	}
}

//-----------------
// copy Constructor
//-----------------
BinaryEjParameter::BinaryEjParameter(const BinaryEjParameter * iParameter) 
: BinaryParameter(iParameter) 
{
	_scatter = new float[_pbDimension];
	float * iScatter = iParameter->getScatter();
	for (int j = 0; j < _pbDimension; j++) {
		_scatter[j] = iScatter[j];
	}
}

//------------------------
// reset to default values
//------------------------
void BinaryEjParameter::reset() {
	for (int j = 0; j < _pbDimension; j++) {
		_scatter[j] = 0;
	}
	BinaryParameter::reset();
}

//---------
// clone 
//---------
Parameter * BinaryEjParameter::clone() const {
	BinaryEjParameter * newParam = new BinaryEjParameter(this);
	return (newParam);
}

//-----------
// Destructor
//-----------
BinaryEjParameter::~BinaryEjParameter() {
	if (_scatter) {
		delete [] _scatter;
	}
}

//---------------------
/// Comparison operator
//---------------------
bool BinaryEjParameter::operator ==(const BinaryEjParameter & param) const {
	if (!BinaryParameter::operator==(param)) return false;
	for (int j = 0; j < _pbDimension; j++) {
		if (_scatter[j] != param.getScatter()[j]) return false;
	}
	return true;
}

//-----------
// getFreeParameter
//-----------
int BinaryEjParameter::getFreeParameter() const {
	int nbFreeParameter = _pbDimension;
	if (_freeProportion) {
		nbFreeParameter += _nbCluster - 1;
	}
	return nbFreeParameter;
}

//-------
// getPdf
//-------
float BinaryEjParameter::getPdf(int iSample, int kCluster) const {
	int j;
	float bernPdf = 1.0;
	BinaryData * data = _model->getBinaryData();
	BinarySample * curSample = (data->_matrix[iSample])->getBinarySample();

	for (j = 0; j < _pbDimension; j++) {
		// iSample have major modality ?//
		if (curSample->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf *= 1.0 - _scatter[j];
		}
		else {
			bernPdf *= _scatter[j] / (_tabNbModality[j] - 1.0);
		}
	}
	return bernPdf;
}

//----------
// getLogPdf
//----------
float BinaryEjParameter::getLogPdf(int iSample, int kCluster) const {
	int j;
	float bernPdf = 0.0;
	BinaryData * data = _model->getBinaryData();
	BinarySample * curSample = (data->_matrix[iSample])->getBinarySample();

	for (j = 0; j < _pbDimension; j++) {
		// iSample have major modality ?//
		if (curSample->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf += log(1.0 - _scatter[j]);
		}
		else {
			bernPdf += log(_scatter[j] / (_tabNbModality[j] - 1.0));
		}
	}
	return bernPdf;
}

//-------
// getPdf
//-------
/* Compute normal probability density function
	   for x vector and kCluster th cluster
 */
float BinaryEjParameter::getPdf(Sample * x, int kCluster) const {
	int j;
	float bernPdf = 1.0;
	BinarySample * binaryX = x->getBinarySample();

	for (j = 0; j < _pbDimension; j++) {
		// iSample have major modality ? //
		if (binaryX->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf *= 1.0 - _scatter[j];
		}
		else {
			bernPdf *= _scatter[j] / (_tabNbModality[j] - 1.0);
		}
	}
	return bernPdf;
}

//--------------------
// getlogLikelihoodOne (one cluster)
//--------------------
float BinaryEjParameter::getLogLikelihoodOne() const {
	int i;
	int j;
	float logLikelihoodOne = 0.0, pdf; //, * Scatter;
	std::unique_ptr<float[]> Scatter(new float[_pbDimension]);
	//int * Center = new int[_pbDimension];
    std::unique_ptr<int[]> Center(new int[_pbDimension]);
	//float * tabNbSampleInMajorModality = new float[_pbDimension];
	std::unique_ptr<float[]> tabNbSampleInMajorModality(new float[_pbDimension]);
	int nbSample = _model->getNbSample();
	BinaryData * data = _model->getBinaryData();

	// Compute Center fo One cluster //
	getTabCenterIfOneCluster(Center.get(), tabNbSampleInMajorModality.get());

	// Compute Scatter for One cluster //
	for (j = 0; j < _pbDimension; j++) {
		//Scatter[j] = 1- (tabNbSampleInMajorModality[j] / data->_weightTotal);
		Scatter[j] = 1 - ((tabNbSampleInMajorModality[j] 
				+ (1. / _tabNbModality[j])) / (data->_weightTotal + 1));
	}

	// Compute the log-likelihood for one cluster (k=1) //
	for (i = 0; i < nbSample; i++) {
      pdf = computePdfOneCluster(data->_matrix[i], Center.get(), Scatter.get(), _tabNbModality);
		logLikelihoodOne += log(pdf) * data->_weight[i];
	}

	//delete[] Center;
	//delete[] Scatter;
	//delete[] tabNbSampleInMajorModality;

	return logLikelihoodOne;
}

//----------------
// Compute scatter 
//----------------
void BinaryEjParameter::computeScatter() {
	int j, k;
	int i;
	float ej; // nb d'individus prenant la modalite majoritaire sur la variable j
	float ** tabCik = _model->getTabCik();

	BinaryData * data = _model->getBinaryData();
	Sample ** dataMatrix = data->getDataMatrix();
	BinarySample * curSample;
	float totalWeight = data->_weightTotal;
	int nbSample = _model->getNbSample();

	for (j = 0; j < _pbDimension; j++) {
		ej = 0.0;
		for (k = 0; k < _nbCluster; k++) {
			for (i = 0; i < nbSample; i++) {
				curSample = dataMatrix[i]->getBinarySample();
				if (curSample->getDataValue(j) == _tabCenter[k][j]) {
					ej += (tabCik[i][k] * data->_weight[i]);
				}
			} // end for i
		} // end for k
		//_scatter[j] = 1 - (ej / totalWeight);
		_scatter[j] = 1 - ((ej + ((_nbCluster * 1.) 
				/ _tabNbModality[j])) / (totalWeight + _nbCluster));
	} // end for j

	/*
  Version issue directement des formules mais moins rapide :
	for (j=0; j<_pbDimension; j++){
	ej = 0.0;
	for (k=0; k<_nbCluster; k++){
	for (i=0; i<nbSample; i++){
	curSample = (XEMBinarySample*)dataMatrix[i];
	if ( curSample->getDataValue(j) == _tabCenter[k][j] ){
	if (tabZikKnown[i])
	ej += (tabZik[i][k] * data->_weight[i]);
	else
	ej += (tabTik[i][k] * data->_weight[i]);
  }
  } // end for i
	ej += 1./_tabNbModality[j];
  } // end for k
	  //_scatter[j] = 1 - (ej / totalWeight);
	_scatter[j] = 1 - (ej  / (totalWeight + _nbCluster));
  } // end for j
	 */
}

//--------------------------
// Compute random scatter(s)
//--------------------------
void BinaryEjParameter::computeRandomScatter() {
	// tirage d'une valeur comprise entre 0 et 1/_tabNbModality[j] 
	for (int j = 0; j < _pbDimension; j++) {
		_scatter[j] = rnd() / _tabNbModality[j];
	}
}

//---------------
//recopy scatter from iParam 
//---------------
// Note : iParam must be a XEMBinaryEjParameter*
void BinaryEjParameter::recopyScatter(Parameter * iParam) {
	if (typeid (*iParam) != typeid (*this)) {
		THROW(OtherException, badBinaryParameterClass);
	}
	float * iScatter = ((BinaryEjParameter*) iParam)->getScatter();
	for (int j = 0; j < _pbDimension; j++) {
		_scatter[j] = iScatter[j];
	}
}

//---------------
//create scatter from Scatter Ekjh 
//---------------
//on fait une moyenne
void BinaryEjParameter::createScatter(float *** scatter) {
	int k, j, h;
	for (j = 0; j < _pbDimension; j++) {
		_scatter[j] = 0.0;
		for (k = 0; k < _nbCluster; k++) {
			h = _tabCenter[k][j];
			_scatter[j] += scatter[k][j][h - 1];
		}
		_scatter[j] /= _nbCluster;
	}
}

//------------
// editScatter (for debug)
//------------
void BinaryEjParameter::editScatter(int k) {
	int j, h;
	for (j = 0; j < _pbDimension; j++) {
		for (h = 1; h <= _tabNbModality[j]; h++) {
			if (h == _tabCenter[k][j]) {
				cout << "\t" << _scatter[j];
			}
			else {
				cout << "\t" << _scatter[j] / (_tabNbModality[j] - 1);
			}
		}
		cout << endl;
	}
}

// editScatter 
//------------
void BinaryEjParameter::editScatter(std::ofstream & oFile, int k, bool text) {
  int j, h;
  if (text) {
    oFile << "\t\t\tScattering : \n";
  }
  for (j = 0; j < _pbDimension; j++) {
    if (text) {
      oFile << "\t\t\t\t\t";
      ;
    }
    for (h = 1; h <= _tabNbModality[j]; h++) {
      if (h == _tabCenter[k][j])
				putFloatInStream(oFile, _scatter[j], "  ");
      else
				putFloatInStream(oFile, _scatter[j] / (_tabNbModality[j] - 1), "  ");
    }
    oFile << endl;
  }
}

// Read Scatter in input file
//---------------------------
void BinaryEjParameter::inputScatter(std::ifstream & fi, int k) {
	THROW(OtherException, internalMixmodError);
}

// Read Scatter in input containers
//---------------------------
void BinaryEjParameter::inputScatter(float *** scatters) {
	THROW(OtherException, internalMixmodError);
}

float *** BinaryEjParameter::scatterToArray() const {
	int k, j, h;
	float *** tabScatter = new float**[_nbCluster];
	for (k = 0; k < _nbCluster; k++) {
		tabScatter[k] = new float*[_pbDimension];
		for (j = 0; j < _pbDimension; j++) {
			tabScatter[k][j] = new float[_tabNbModality[j]];
			for (h = 1; h <= _tabNbModality[j]; h++) {
				if (h == _tabCenter[k][j]) {
					tabScatter[k][j][h - 1] = _scatter[j];
				}
				else {
					tabScatter[k][j][h - 1] = _scatter[j] / (_tabNbModality[j] - 1);
				}
			}
		}
	}
	return tabScatter;
}

}
