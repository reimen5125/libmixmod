/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/BinaryEParameter.cpp  description
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
#include "mixmod/Kernel/Parameter/BinaryEParameter.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/BinarySample.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Utilities/Random.h"

namespace XEM {

//--------------------
// Default constructor
//--------------------
BinaryEParameter::BinaryEParameter() {
	THROW(OtherException, wrongConstructorType);
}

//-------------------------------
// Constructor called by XEMModel
//-------------------------------
BinaryEParameter::BinaryEParameter(
		Model * iModel, ModelType * iModelType, int * tabNbModality) 
: BinaryParameter(iModel, iModelType, tabNbModality) 
{
	_scatter = 0;
}

//-----------------
// copy Constructor
//-----------------
BinaryEParameter::BinaryEParameter(const BinaryEParameter * iParameter) 
: BinaryParameter(iParameter) 
{
	_scatter = iParameter->getScatter();
}

//---------
// clone 
//---------
Parameter * BinaryEParameter::clone() const {
	BinaryEParameter * newParam = new BinaryEParameter(this);
	return (newParam);
}

//-----------
// Destructor
//-----------
BinaryEParameter::~BinaryEParameter() {
}

//---------------------
/// Comparison operator
//---------------------
bool BinaryEParameter::operator ==(const BinaryEParameter & param) const {
	if (!BinaryParameter::operator==(param)) return false;
	if (_scatter != param.getScatter()) return false;
	return true;
}

//------------------------
// reset to default values
//------------------------
void BinaryEParameter::reset() {
	_scatter = 0.0;
	BinaryParameter::reset();
}

//-----------
// getFreeParameter
//-----------
int BinaryEParameter::getFreeParameter() const {
	int nbFreeParameter = 1;
	if (_freeProportion) {
		nbFreeParameter += _nbCluster - 1;
	}
	return nbFreeParameter;
}

//-------
// getPdf
//-------
float BinaryEParameter::getPdf(int iSample, int kCluster) const {
	int j;
	float bernPdf = 1.0;
	BinaryData * data = _model->getBinaryData();
	BinarySample * curSample = (data->_matrix[iSample])->getBinarySample();


	for (j = 0; j < _pbDimension; j++) {
		// iSample have major modality ?//
		if (curSample->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf *= 1.0 - _scatter;
		}
		else {
			bernPdf *= _scatter / (_tabNbModality[j] - 1.0);
		}
	}
	return bernPdf;
}

//----------
// getLogPdf
//----------
float BinaryEParameter::getLogPdf(int iSample, int kCluster) const {
	int j;
	float bernPdf = 0.0;
	BinaryData * data = _model->getBinaryData();
	BinarySample * curSample = (data->_matrix[iSample])->getBinarySample();


	for (j = 0; j < _pbDimension; j++) {
		// iSample have major modality ?//
		if (curSample->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf += logf(1.0 - _scatter);
		}
		else {
			bernPdf += logf(_scatter / (_tabNbModality[j] - 1.0));
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
float BinaryEParameter::getPdf(Sample * x, int kCluster) const {
	int j;
	float bernPdf = 1.0;
	BinarySample * binaryX = x->getBinarySample();

	for (j = 0; j < _pbDimension; j++) {
		// iSample have major modality ? //
		if (binaryX->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf *= 1.0 - _scatter;
		}
		else {
			bernPdf *= _scatter / (_tabNbModality[j] - 1.0);
		}
	}
	return bernPdf;
}

//--------------------
// getlogLikelihoodOne (one cluster)
//--------------------
float BinaryEParameter::getLogLikelihoodOne() const {
	int i;
	int j;
	float logLikelihoodOne = 0.0, value, pdf, Scatter;
	//int * Center = new int[_pbDimension];
	//float* tabNbSampleInMajorModality = new float[_pbDimension];
	std::unique_ptr<int[]> Center(new int[_pbDimension]);
	std::unique_ptr<float[]> tabNbSampleInMajorModality(new float[_pbDimension]);
	int nbSample = _model->getNbSample();
	BinaryData * data = _model->getBinaryData();

	// Compute Center fo One cluster //
	getTabCenterIfOneCluster(Center.get(), tabNbSampleInMajorModality.get());

	// Compute Scatter for One cluster //
	value = 0.0;
	for (j = 0; j < _pbDimension; j++) {
		//value += tabNbSampleInMajorModality[j];
		value += tabNbSampleInMajorModality[j] + 1. / _tabNbModality[j];
	}
	Scatter = 1 - (value / ((data->_weightTotal + 1) * _pbDimension));

	// Compute the log-likelihood for one cluster (k=1) //
	//--------------------------------------------------//
	for (i = 0; i < nbSample; i++) {
      pdf = computePdfOneCluster(data->_matrix[i], Center.get(), Scatter, _tabNbModality);
		logLikelihoodOne += logf(pdf) * data->_weight[i];
	}

	//delete[] Center;
	//delete[] tabNbSampleInMajorModality;

	return logLikelihoodOne;
}

//----------------
// Compute scatter 
//----------------
void BinaryEParameter::computeScatter() {
	int j, k;
	int i;
	float ** tabCik = _model->getTabCik();

	BinaryData * data = _model->getBinaryData();
	Sample ** dataMatrix = data->getDataMatrix();
	BinarySample * curSample;
	float totalWeight = data->_weightTotal;
	int nbSample = _model->getNbSample();

	float e = 0.0; // nb d'individus prenant la modalite majoritaire 
	                // (pour toutes les classes, pour toutes les variables)
	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			for (i = 0; i < nbSample; i++) {
				curSample = dataMatrix[i]->getBinarySample();
				if (curSample->getDataValue(j) == _tabCenter[k][j]) {
					e += (tabCik[i][k] * data->_weight[i]);
				}
			} // end for i
			e += 1. / _tabNbModality[j];
		} // end for j
	} // end for k
	//_scatter = 1- (e / (totalWeight *_pbDimension));
	_scatter = 1 - (e / ((totalWeight + _nbCluster) * _pbDimension));
}

//----------------------------------------------------
// Compute scatter(s)  as if there was only one cluster
//---------------------------------------------------
/*void XEMBinaryEParameter::computeScatterIfOneCluster(float totalWeight, 
     float * tabNbSampleInMajorModality, float ** tabNbSamplePerModality){
  float value = 0.0;
  for (int j=0; j<_pbDimension; j++){
	cout<<"tabNbSampleInMajorModality[j] : "<<tabNbSampleInMajorModality[j]<<endl;
	value += tabNbSampleInMajorModality[j];
}
  _scatter = 1-  (value / (totalWeight * _pbDimension));
}
 */

//---------------------------
// Compute random scatter(s) 
//---------------------------
void BinaryEParameter::computeRandomScatter() {
	int minNbModality = _tabNbModality[0];
	for (int j = 1; j < _pbDimension; j++) {
		if (_tabNbModality[j] < minNbModality) {
			minNbModality = _tabNbModality[j];
		}
	}

	// tirage d'une valeur comprise entre 0 et 1/minNbModality
	_scatter = rnd() / minNbModality;
}

//---------------
//recopy scatter from iParam 
//---------------
// Note : iParam must be a BinaryEParameter*
void BinaryEParameter::recopyScatter(Parameter * iParam) {
	if (typeid (*iParam) != typeid (*this)) {
		THROW(OtherException, badBinaryParameterClass);
	}
	float iScatter = ((BinaryEParameter*) iParam)->getScatter();
	_scatter = iScatter;
}

//---------------
//create scatter from Scatter Ekjh 
//---------------
//on fait une moyenne
void BinaryEParameter::createScatter(float *** scatter) {
	_scatter = 0.0;
	int k, j, h;
	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			h = _tabCenter[k][j];
			_scatter += scatter[k][j][h - 1];
		}
	}
	_scatter /= (_nbCluster * _pbDimension);
}

//------------
// editScatter (for debug)
//------------
void BinaryEParameter::editScatter(int k) {
	int j, h;
	for (j = 0; j < _pbDimension; j++) {
		for (h = 1; h <= _tabNbModality[j]; h++) {
			if (h == _tabCenter[k][j]) {
				cout << "\t" << _scatter;
			}
			else {
				cout << "\t" << _scatter / (_tabNbModality[j] - 1);
			}
		}
		cout << endl;
	}
}

// editScatter 
//------------
void BinaryEParameter::editScatter(std::ofstream & oFile, int k, bool text) {
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
				putFloatInStream(oFile, _scatter, "  ");
      else
				putFloatInStream(oFile, _scatter / (_tabNbModality[j] - 1), "  ");
		}
		oFile << endl;
	}
}

// Read Scatter in input file
//---------------------------
void BinaryEParameter::inputScatter(std::ifstream & fi, int k) {
	THROW(OtherException, internalMixmodError);
}

// Read Scatter in input containers
//---------------------------
void BinaryEParameter::inputScatter(float *** scatters) {
	THROW(OtherException, internalMixmodError);
}

float *** BinaryEParameter::scatterToArray() const {
	float*** tabScatter = new float**[_nbCluster];
	int k, j, h;
	for (k = 0; k < _nbCluster; k++) {
		tabScatter[k] = new float*[_pbDimension];
		for (j = 0; j < _pbDimension; j++) {
			tabScatter[k][j] = new float[_tabNbModality[j] ];
			for (h = 1; h <= _tabNbModality[j]; h++) {
				if (h == _tabCenter[k][j]) {
					tabScatter[k][j][h - 1] = _scatter;
				}
				else {
					tabScatter[k][j][h - 1] = _scatter / (_tabNbModality[j] - 1);
				}
			}
		}
	}

	return tabScatter;
}

}
