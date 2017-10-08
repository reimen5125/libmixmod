/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/BinaryEkjhParameter.cpp  description
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
#include "mixmod/Kernel/Parameter/BinaryEkjhParameter.h"
#include "mixmod/Kernel/IO/BinaryData.h"
#include "mixmod/Kernel/IO/BinarySample.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Utilities/Random.h"

namespace XEM {

//--------------------
// Default constructor
//--------------------
BinaryEkjhParameter::BinaryEkjhParameter() {
	THROW(OtherException, wrongConstructorType);
}

//-------------------------------
// Constructor called by XEMModel
//-------------------------------
BinaryEkjhParameter::BinaryEkjhParameter(
		Model * iModel, ModelType * iModelType, int * tabNbModality) 
: BinaryParameter(iModel, iModelType, tabNbModality) 
{
	_scatter = new float**[_nbCluster];
	int k, j, h;
	for (k = 0; k < _nbCluster; k++) {
		_scatter[k] = new float*[_pbDimension];
		for (j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = new float[_tabNbModality[j]];
			for (h = 0; h < _tabNbModality[j]; h++) {
				_scatter[k][j][h] = 0.0;
			}
		}
	}
}

//-----------------
// copy Constructor
//-----------------
BinaryEkjhParameter::BinaryEkjhParameter(const BinaryEkjhParameter * iParameter) 
: BinaryParameter(iParameter) 
{
	int k, j, h;
	_scatter = new float**[_nbCluster];

	for (k = 0; k < _nbCluster; k++) {
		_scatter[k] = new float*[_pbDimension];
		for (j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = new float[_tabNbModality[j]];
		}
	}
	float *** iScatter = iParameter->getScatter();
	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			for (h = 0; h < _tabNbModality[j]; h++) {
				_scatter[k][j][h] = iScatter[k][j][h];
			}
		}
	}
}

//------------
// Constructor
// called by XEMStrategyType if USER partition
//------------
BinaryEkjhParameter::BinaryEkjhParameter(int iNbCluster, int iPbDimension, 
		ModelType * iModelType, int * tabNbModality, std::string & iFileName) 
: BinaryParameter(iNbCluster, iPbDimension, iModelType, tabNbModality) 
{
	_scatter = new float**[_nbCluster];
	for (int k = 0; k < _nbCluster; k++) {
		_scatter[k] = new float*[_pbDimension];
		for (int j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = new float[_tabNbModality[j]];
		}
	}

	if (iFileName.compare("") != 0) {
		std::ifstream paramFile(iFileName.c_str(), ios::in);
		if (!paramFile.is_open()) {
			THROW(InputException, wrongParamFileName);
		}
		input(paramFile);
		paramFile.close();
	}
}

//------------
// Constructor
//------------
BinaryEkjhParameter::BinaryEkjhParameter(
		int iNbCluster, 
		int iPbDimension, 
		ModelType * iModelType, 
		int * tabNbModality, 
		float * proportions, 
		float ** centers, 
		float *** scatters)
: BinaryParameter(iNbCluster, iPbDimension, iModelType, tabNbModality) 
{
	_scatter = new float**[_nbCluster];
	for (int k = 0; k < _nbCluster; k++) {
		_scatter[k] = new float*[_pbDimension];
		for (int j = 0; j < _pbDimension; j++) {
			_scatter[k][j] = new float[_tabNbModality[j]];
		}
	}
	input(proportions, centers, scatters);
}

//---------
// clone 
//---------
Parameter * BinaryEkjhParameter::clone() const {
	BinaryEkjhParameter * newParam = new BinaryEkjhParameter(this);
	return (newParam);
}

//-----------
// Destructor
//-----------
BinaryEkjhParameter::~BinaryEkjhParameter() {
	if (_scatter) {
		for (int k = 0; k < _nbCluster; k++) {
			for (int j = 0; j < _pbDimension; j++) {
				delete [] _scatter[k][j];
			}
			delete [] _scatter[k];
		}
		delete [] _scatter;
	}
	_scatter = NULL;
}

//---------------------
/// Comparison operator
//---------------------
bool BinaryEkjhParameter::operator ==(const BinaryEkjhParameter & param) const {
	if (!BinaryParameter::operator==(param)) return false;
	for (int k = 0; k < _nbCluster; k++) {
		for (int j = 0; j < _pbDimension; j++) {
			for (int h = 0; h < _tabNbModality[j]; h++) {
				if (_scatter[k][j][h] != param.getScatter()[k][j][h]) return false;
			}
		}
	}
	return true;
}

//------------------------
// reset to default values
//------------------------
void BinaryEkjhParameter::reset() {
	int k, j, h;
	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			for (h = 0; h < _tabNbModality[j]; h++) {
				_scatter[k][j][h] = 0.0;
			}
		}
	}
	BinaryParameter::reset();
}

//-----------
// getFreeParameter
//-----------
int BinaryEkjhParameter::getFreeParameter() const {
	int nbFreeParameter = 0;
	for (int j = 0; j < _pbDimension; j++) {
		nbFreeParameter += _tabNbModality[j] - 1;
	}
	nbFreeParameter *= _nbCluster;
	if (_freeProportion) {
		nbFreeParameter += _nbCluster - 1;
	}
	return nbFreeParameter;
}

//-------
// getPdf
//-------
float BinaryEkjhParameter::getPdf(int iSample, int kCluster) const {
	int j;
	float bernPdf = 1.0;
	BinaryData * data = _model->getBinaryData();
	BinarySample * curSample = (data->_matrix[iSample])->getBinarySample();
	int value;

	for (j = 0; j < _pbDimension; j++) {
		// iSample have major modality ?//
		value = curSample->getDataValue(j);
		if (value == _tabCenter[kCluster][j]) {
			bernPdf *= 1.0 - _scatter[kCluster][j][value - 1];
		}
		else {
			bernPdf *= _scatter[kCluster][j][value - 1];
		}
	}
	return bernPdf;
}

//----------
// getLogPdf
//----------
float BinaryEkjhParameter::getLogPdf(int iSample, int kCluster) const {
	int j;
	float bernPdf = 0.0;
	BinaryData * data = _model->getBinaryData();
	BinarySample * curSample = (data->_matrix[iSample])->getBinarySample();
	int value;

	for (j = 0; j < _pbDimension; j++) {
		// iSample have major modality ?//
		value = curSample->getDataValue(j);
		if (value == _tabCenter[kCluster][j]) {
			bernPdf += log(1.0 - _scatter[kCluster][j][value - 1]);
		}
		else {
			bernPdf += log(_scatter[kCluster][j][value - 1]);
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
float BinaryEkjhParameter::getPdf(Sample * x, int kCluster) const {
	int j;
	float bernPdf = 1.0;
	BinarySample * binaryX = x->getBinarySample();
	int value;

	for (j = 0; j < _pbDimension; j++) {
		value = binaryX->getDataValue(j);
		// iSample have major modality ? //
		if (value == _tabCenter[kCluster][j]) {
			bernPdf *= 1.0 - _scatter[kCluster][j][value - 1];
		}
		else {
			bernPdf *= _scatter[kCluster][j][value - 1];
		}
	}
	return bernPdf;
}

//--------------------
// getlogLikelihoodOne (one cluster)
//--------------------
float BinaryEkjhParameter::getLogLikelihoodOne() const {
	int i;
	int j, h;
	float logLikelihoodOne = 0.0, pdf; //, ** Scatter;
	//Scatter = new float*[_pbDimension];
    std::unique_ptr<float*[], TabDeleter<float>>  Scatter(new float*[_pbDimension], TabDeleter<float>(_pbDimension));
	for (j = 0; j < _pbDimension; j++) {
		Scatter[j] = new float[_tabNbModality[j]];
	}

	//int * Center = new int[_pbDimension];
    std::unique_ptr<int[]> Center(new int[_pbDimension]);
	//float * tabNbSampleInMajorModality = new float [_pbDimension];
    std::unique_ptr<float[]> tabNbSampleInMajorModality(new float[_pbDimension]);
	//float ** tabNbSamplePerModality = new float * [_pbDimension];
    std::unique_ptr<float*[], TabDeleter<float>>  tabNbSamplePerModality(new float * [_pbDimension], TabDeleter<float>(_pbDimension));
	for (j = 0; j < _pbDimension; j++) {
		tabNbSamplePerModality[j] = new float[_tabNbModality[j]];
	}
	int nbSample = _model->getNbSample();
	BinaryData * data = _model->getBinaryData();

	// Compute Center fo One cluster //
	getTabCenterIfOneCluster(Center.get(), tabNbSampleInMajorModality.get(), tabNbSamplePerModality.get());


	// Compute Scatter for One cluster //
	for (j = 0; j < _pbDimension; j++) {
		for (h = 0; h < _tabNbModality[j]; h++) {
			if (h + 1 == Center[j]) {
				//Scatter[j][h] = 1 - (tabNbSampleInMajorModality[j] / data->_weightTotal);
				Scatter[j][h] = 1 - (tabNbSampleInMajorModality[j] 
						+ (1. / _tabNbModality[j])) / (data->_weightTotal + 1);
			}
			else {
				//Scatter[j][h] = tabNbSamplePerModality[j][h] / data->_weightTotal;
				Scatter[j][h] = (tabNbSamplePerModality[j][h] 
						+ (1. / _tabNbModality[j])) / (data->_weightTotal + 1);
			}
		}
	}

	// Compute the log-likelihood for one cluster (k=1) //
	//--------------------------------------------------//
	for (i = 0; i < nbSample; i++) {
      pdf = computePdfOneCluster(data->_matrix[i], Center.get(), Scatter.get(), _tabNbModality);
		logLikelihoodOne += log(pdf) * data->_weight[i];
	}

	//delete[] Center;
	//for (j = 0; j < _pbDimension; j++) {
	//	delete[] Scatter[j];
	//}
	//delete [] Scatter;
	//for (j = 0; j < _pbDimension; j++) {
	//	delete[] tabNbSamplePerModality[j];
	//}
	//delete[] tabNbSamplePerModality;
	//delete[] tabNbSampleInMajorModality;

	return logLikelihoodOne;
}

//----------------
// Compute scatter 
//----------------
void BinaryEkjhParameter::computeScatter() {
	int j, k, h;
	int i;
	float valueScatter;
	float * tabNk = _model->getTabNk();
	float ** tabCik = _model->getTabCik();

	BinaryData * data = _model->getBinaryData();
	Sample ** dataMatrix = data->getDataMatrix();
	BinarySample * curSample;
	int nbSample = _model->getNbSample();
	int value;

	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			for (h = 0; h < _tabNbModality[j]; h++) {
				valueScatter = 0.0;
				for (i = 0; i < nbSample; i++) {
					curSample = dataMatrix[i]->getBinarySample();
					value = curSample->getDataValue(j);
					if (value == h + 1) {
						valueScatter += (tabCik[i][k] * data->_weight[i]);
					}
				}
				if (h + 1 == _tabCenter[k][j]) {
					//_scatter[k][j][h] = (tabNk[k] - valueScatter) / tabNk[k];
					_scatter[k][j][h] = 1 
							- ((valueScatter + 1. / _tabNbModality[j]) / (tabNk[k] + 1));
				}
				else {
					//_scatter[k][j][h] = valueScatter / tabNk[k];
					_scatter[k][j][h] = (valueScatter + 1. / _tabNbModality[j]) / (tabNk[k] + 1);
				}
			}// end for h
		} // end for j
	} // end for k
}

//--------------------------
// Compute random scatter(s)
//--------------------------
void BinaryEkjhParameter::computeRandomScatter() {
	for (int k = 0; k < _nbCluster; k++) {
		for (int j = 0; j < _pbDimension; j++) {
			float scatterKJOnMajorModality = rnd() / _tabNbModality[j];
			for (int h = 0; h < _tabNbModality[j]; h++) {
				if (h + 1 == _tabCenter[k][j]) {// on est sur le centre       
					_scatter[k][j][h] = scatterKJOnMajorModality;
				}
				else {
					_scatter[k][j][h] = scatterKJOnMajorModality / (_tabNbModality[j] - 1);
				}
			}
		}
	}
}

//---------------
//recopy scatter from iParam 
//---------------
// Note : iParam must be a XEMBinaryEkjhParameter*
void BinaryEkjhParameter::recopyScatter(Parameter * iParam) {
	if (typeid (*iParam) != typeid (*this)) {
		THROW(OtherException, badBinaryParameterClass);
	}
	float *** iScatter = ((BinaryEkjhParameter*) iParam)->getScatter();
	int k, j, h;
	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			for (h = 0; h < _tabNbModality[j]; h++) {
				_scatter[k][j][h] = iScatter[k][j][h];
			}
		}
	}
}

//---------------
//create scatter from Scatter Ekjh 
//---------------
void BinaryEkjhParameter::createScatter(float *** scatter) {
	int k, j, h;
	for (k = 0; k < _nbCluster; k++) {
		for (j = 0; j < _pbDimension; j++) {
			for (h = 0; h < _tabNbModality[j]; h++) {
				_scatter[k][j][h] = scatter[k][j][h];
			}
		}
	}
}

//------------
// editScatter (for debug)
//------------
void BinaryEkjhParameter::editScatter(int k) {
	int j, h;
	for (j = 0; j < _pbDimension; j++) {
		for (h = 0; h < _tabNbModality[j]; h++) {
			cout << "\t" << _scatter[k][j][h];
		}
		cout << endl;
	}
}

// editScatter 
//------------
void BinaryEkjhParameter::editScatter(std::ofstream & oFile, int k, bool text) {
	int j, h;
	if (text) {
		oFile << "\t\t\tScattering : \n";
	}
	for (j = 0; j < _pbDimension; j++) {
		if (text) {
			oFile << "\t\t\t\t\t";
			;
		}
		for (h = 0; h < _tabNbModality[j]; h++)
			putFloatInStream(oFile, _scatter[k][j][h], "  ");
		oFile << endl;
	}
}

// Read Scatter in input file
//---------------------------
void BinaryEkjhParameter::inputScatter(std::ifstream & fi, int k) {
  //Mise a jour de tabCenter, tabScatter et tabProportion
  int j, h;
  //for (k = 0; k < _nbCluster; k++) {
  //  // Proportions //
  //  fi >> _tabProportion[k];

  //  // Center  //
  //  for (j = 0; j < _pbDimension; j++) {
  //    fi >> _tabCenter[k][j];
  //  }

  // Scatter  //
  for (j = 0; j < _pbDimension; j++) {
    for (h = 0; h < _tabNbModality[j]; h++)
			_scatter[k][j][h] = getFloatFromStream(fi);
  }
  //} // end for k
}

// Read Scatter in input containers
//---------------------------
void BinaryEkjhParameter::inputScatter(float *** scatters) {

	// Scatter
	for (int k = 0; k < _nbCluster; k++) {
		for (int j = 0; j < _pbDimension; j++) {
			for (int h = 0; h < _tabNbModality[j]; h++) {
				_scatter[k][j][h] = scatters[k][j][h];
			}
		}
	} // end for k
}

float*** BinaryEkjhParameter::scatterToArray() const {
	float *** tabScatter = new float**[_nbCluster];
	int k, j;
	for (k = 0; k < _nbCluster; ++k) {
		tabScatter[k] = new float*[_pbDimension];
		for (j = 0; j < _pbDimension; ++j) {
			tabScatter[k][j] = new float [_tabNbModality[j]];
			recopyTab<float>(_scatter[k][j], tabScatter[k][j], _tabNbModality[j]);
		}
	}

	return tabScatter;
}

}
