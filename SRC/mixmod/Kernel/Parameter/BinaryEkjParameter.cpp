/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/BinaryEkjParameter.cpp  description
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
#include "mixmod/Kernel/Parameter/BinaryEkjParameter.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/BinarySample.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Utilities/Random.h"

namespace XEM {

//--------------------
// Default constructor
//--------------------
BinaryEkjParameter::BinaryEkjParameter() {
	THROW(OtherException, wrongConstructorType);
}

//-------------------------------
// Constructor called by XEMModel
//-------------------------------
BinaryEkjParameter::BinaryEkjParameter(
		Model * iModel, ModelType * iModelType, int * tabNbModality) 
: BinaryParameter(iModel, iModelType, tabNbModality) 
{
	_scatter = new float*[_nbCluster];
	for (int k = 0; k < _nbCluster; k++) {
		_scatter[k] = new float[_pbDimension];
		for (int j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = 0.0;
		}
	}
}

//-----------------
// copy Constructor
//-----------------
BinaryEkjParameter::BinaryEkjParameter(const BinaryEkjParameter * iParameter) 
: BinaryParameter(iParameter) 
{
	_scatter = new float*[_nbCluster];
	for (int k = 0; k < _nbCluster; k++) {
		_scatter[k] = new float[_pbDimension];
	}
	float ** iScatter = iParameter->getScatter();
	for (int k = 0; k < _nbCluster; k++) {
		for (int j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = iScatter[k][j];
		}
	}
}

//---------
// clone 
//---------
Parameter * BinaryEkjParameter::clone() const {
	BinaryEkjParameter * newParam = new BinaryEkjParameter(this);
	return (newParam);
}

//-----------
// Destructor
//-----------
BinaryEkjParameter::~BinaryEkjParameter() {
	if (_scatter) {
		for (int k = 0; k < _nbCluster; k++) {
			delete [] _scatter[k];
		}
	}
	delete [] _scatter;
	_scatter = NULL;
}

//---------------------
/// Comparison operator
//---------------------
bool BinaryEkjParameter::operator ==(const BinaryEkjParameter & param) const {
	if (!BinaryParameter::operator==(param)) return false;
	for (int k = 0; k < _nbCluster; k++) {
		for (int j = 0; j < _pbDimension; j++) {
			if (_scatter[k][j] != param.getScatter()[k][j]) return false;
		}
	}
	return true;
}

//------------------------
// reset to default values
//------------------------
void BinaryEkjParameter::reset() {
	int k, j;
	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = 0.0;
		}
	}
	BinaryParameter::reset();
}

//-----------
// getFreeParameter
//-----------
int BinaryEkjParameter::getFreeParameter() const {
	int nbFreeParameter = _pbDimension * _nbCluster;
	if (_freeProportion) {
		nbFreeParameter += _nbCluster - 1;
	}
	return nbFreeParameter;
}

//-------
// getPdf
//-------
float BinaryEkjParameter::getPdf(int iSample, int kCluster) const {
	//cout<<" XEMBinaryEkjParameter::getPdf"<<endl;
	int j;
	float bernPdf = 1.0;
	BinaryData * data = _model->getBinaryData();
	BinarySample * curSample = (data->_matrix[iSample])->getBinarySample();

	for (j = 0; j < _pbDimension; j++) {
		//cout<<"curSample :  "<<curSample->getDataValue(j)<<endl;
		//cout<<" _tabCenter[kCluster][j] :  "<< _tabCenter[kCluster][j]<<endl;
		// iSample have major modality ?//
		if (curSample->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf *= 1.0 - _scatter[kCluster][j];
		}
		else {
			bernPdf *= _scatter[kCluster][j] / (_tabNbModality[j] - 1.0);
		}
	}
	return bernPdf;
}

//----------
// getLogPdf
//----------
float BinaryEkjParameter::getLogPdf(int iSample, int kCluster) const {
	//cout<<" XEMBinaryEkjParameter::getPdf"<<endl;
	int j;
	float bernPdf = 0.0;
	BinaryData * data = _model->getBinaryData();
	BinarySample * curSample = (data->_matrix[iSample])->getBinarySample();

	for (j = 0; j < _pbDimension; j++) {
		//cout<<"curSample :  "<<curSample->getDataValue(j)<<endl;
		//cout<<" _tabCenter[kCluster][j] :  "<< _tabCenter[kCluster][j]<<endl;
		// iSample have major modality ?//
		if (curSample->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf += log(1.0 - _scatter[kCluster][j]);
		}
		else {
			bernPdf += log(_scatter[kCluster][j] / (_tabNbModality[j] - 1.0));
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
float BinaryEkjParameter::getPdf(Sample * x, int kCluster) const {
	int j;
	float bernPdf = 1.0;
	BinarySample * binaryX = x->getBinarySample();

	for (j = 0; j < _pbDimension; j++) {
		// iSample have major modality ? //
		if (binaryX->getDataValue(j) == _tabCenter[kCluster][j]) {
			bernPdf *= 1.0 - _scatter[kCluster][j];
		}
		else {
			bernPdf *= _scatter[kCluster][j] / (_tabNbModality[j] - 1.0);
		}
	}
	return bernPdf;
}

//--------------------
// getlogLikelihoodOne (one cluster)
//--------------------
float BinaryEkjParameter::getLogLikelihoodOne() const {
	int i;
	int j;
	float logLikelihoodOne = 0.0, pdf;//, * Scatter;
	//Scatter = new float[_pbDimension];
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
	for (j = 0; j < _pbDimension; j++)
		//Scatter[j] = 1- (tabNbSampleInMajorModality[j] / data->_weightTotal);
		Scatter[j] = 1 - ((tabNbSampleInMajorModality[j] 
				+ 1. / _tabNbModality[j]) / (data->_weightTotal + 1));

	// Compute the log-likelihood for one cluster (k=1) //
	//--------------------------------------------------//
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
void BinaryEkjParameter::computeScatter() {
	int j, k;
	int i;
	float ekj; // nb d'individus de la classe k prenant la modalite maj sur la variable j
	float * tabNk = _model->getTabNk();
	float ** tabCik = _model->getTabCik();

	BinaryData * data = _model->getBinaryData();
	Sample ** dataMatrix = data->getDataMatrix();
	BinarySample * curSample;
	int nbSample = _model->getNbSample();

	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			ekj = 0.0;
			for (i = 0; i < nbSample; i++) {
				curSample = dataMatrix[i]->getBinarySample();
				if (curSample->getDataValue(j) == _tabCenter[k][j]) {
					ekj += (tabCik[i][k] * data->_weight[i]);
				}
			}
			//_scatter[k][j] = 1 - (ekj / tabNk[k]);
			_scatter[k][j] = 1 - ((ekj + 1. / _tabNbModality[j]) / (tabNk[k] + 1));
		} // end for j
	} // end for k
}

//--------------------------------------------------
// Compute scatter(s)  as if there was only one cluster
//---------------------------------------------------
/*void XEMBinaryEkjParameter::computeScatterIfOneCluster(float totalWeight, 
      float * tabNbSampleInMajorModality, float ** tabNbSamplePerModality){
  for (int k=0; k<_nbCluster; k++){
	for (int j=0; j<_pbDimension; j++){
	  _scatter[k][j] = 1 - (tabNbSampleInMajorModality[j] / totalWeight);
	}
  }
}*/

//--------------------------
// Compute random scatter(s)
//--------------------------
void BinaryEkjParameter::computeRandomScatter() {
	for (int k = 0; k < _nbCluster; k++) {
		// tirage d'une valeur comprise entre 0 et 1./_tabNbModality[j] 
		for (int j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = rnd() / _tabNbModality[j];
		}
	}
}

//---------------
//recopy scatter from iParam 
//---------------
// Note : iParam must be a XEMBinaryEkjParameter*
void BinaryEkjParameter::recopyScatter(Parameter * iParam) {
	if (typeid (*iParam) != typeid (*this)) {
		THROW(OtherException, badBinaryParameterClass);
	}
	float ** iScatter = ((BinaryEkjParameter*) iParam)->getScatter();
	for (int k = 0; k < _nbCluster; k++) {
		for (int j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = iScatter[k][j];
		}
	}
}

//---------------
//create scatter from Scatter Ekjh 
//---------------
void BinaryEkjParameter::createScatter(float *** scatter) {
	int k, j, h;
	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			h = _tabCenter[k][j];
			_scatter[k][j] = scatter[k][j][h - 1];
		}
	}
}

//------------
// editScatter (for debug)
//------------
void BinaryEkjParameter::editScatter(int k) {
	int j, h;
	for (j = 0; j < _pbDimension; j++) {
		for (h = 1; h <= _tabNbModality[j]; h++) {
			if (h == _tabCenter[k][j]) {
				cout << "\t" << _scatter[k][j];
			}
			else {
				cout << "\t" << _scatter[k][j] / (_tabNbModality[j] - 1);
			}
		}
		cout << endl;
	}
}

// editScatter 
//------------
void BinaryEkjParameter::editScatter(std::ofstream & oFile, int k, bool text) {
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
				putFloatInStream(oFile, _scatter[k][j], "  ");
      else
				putFloatInStream(oFile, _scatter[k][j] / (_tabNbModality[j] - 1), "  ");
		}
		oFile << endl;
	}
}

// Read Scatter in input file
//---------------------------
void BinaryEkjParameter::inputScatter(std::ifstream & fi, int k) {
	THROW(OtherException, internalMixmodError);
}

// Read Scatter in input containers
//---------------------------
void BinaryEkjParameter::inputScatter(float *** scatters) {
	THROW(OtherException, internalMixmodError);
}

float *** BinaryEkjParameter::scatterToArray() const {
	int k, j, h;
	float *** tabScatter = new float**[_nbCluster];
	for (k = 0; k < _nbCluster; k++) {
		tabScatter[k] = new float*[_pbDimension];
		for (j = 0; j < _pbDimension; j++) {
			tabScatter[k][j] = new float[_tabNbModality[j]];
			for (h = 1; h <= _tabNbModality[j]; h++) {
				if (h == _tabCenter[k][j]) {
					tabScatter[k][j][h - 1] = _scatter[k][j];
				}
				else {
					tabScatter[k][j][h - 1] = _scatter[k][j] / (_tabNbModality[j] - 1);
				}
			}
		}
	}
	return tabScatter;
}

}
