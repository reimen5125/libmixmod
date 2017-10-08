/***************************************************************************
                             SRC/mixmod/Kernel/Parameter/GaussianEDDAParameter.cpp  description
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
#include "mixmod/Kernel/Parameter/GaussianEDDAParameter.h"
#include "mixmod/Kernel/Parameter/GaussianParameter.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Kernel/IO/GaussianSample.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/Parameter/GaussianGeneralParameter.h"
#include "mixmod/Utilities/Util.h"
#include "mixmod/Kernel/IO/GaussianData.h"
#include "mixmod/Utilities/Random.h"
#include "mixmod/Matrix/Matrix.h"
#include "mixmod/Matrix/DiagMatrix.h"
#include "mixmod/Matrix/SphericalMatrix.h"
#include "mixmod/Matrix/SymmetricMatrix.h"
#include "mixmod/Matrix/GeneralMatrix.h"

namespace XEM {

//-------------------- Constructors--------------

GaussianEDDAParameter::GaussianEDDAParameter() : GaussianParameter() {
	THROW(OtherException, wrongConstructorType);
}

//--------------
// constructor 
//-------------
GaussianEDDAParameter::GaussianEDDAParameter(Model * iModel, ModelType * iModelType) 
: GaussianParameter(iModel, iModelType) 
{
	_tabInvSqrtDetSigma = new float [_nbCluster];
	for (int k = 0; k < _nbCluster; k++) {
		_tabInvSqrtDetSigma[k] = 1.0;
	}
	_tabInvSigma = new Matrix* [_nbCluster];
	_tabSigma = new Matrix* [_nbCluster];
}

//------------------------------------------
// constructor called if USER_initialisation
//-------------------------------------------
GaussianEDDAParameter::GaussianEDDAParameter(
		int iNbCluster, int iPbDimension, ModelType * iModelType) 
: GaussianParameter(iNbCluster, iPbDimension, iModelType) 
{
	_tabInvSqrtDetSigma = new float [_nbCluster];
	for (int k = 0; k < _nbCluster; k++) {
		_tabInvSqrtDetSigma[k] = 0.0;
	}
	_tabInvSigma = new Matrix* [_nbCluster];
	_tabSigma = new Matrix* [_nbCluster];
}

//-----------------
// copy constructor
//-----------------
GaussianEDDAParameter::GaussianEDDAParameter(const GaussianEDDAParameter * iParameter) 
: GaussianParameter(iParameter) 
{
	int k;
	_tabInvSqrtDetSigma = new float[_nbCluster];
	float * iTabInvSqrtDetSigma = iParameter->getTabInvSqrtDetSigma();
	for (k = 0; k < _nbCluster; k++) {
		_tabInvSqrtDetSigma[k] = iTabInvSqrtDetSigma[k];
	}
	_tabInvSigma = new Matrix* [_nbCluster];
	_tabSigma = new Matrix* [_nbCluster];
}

/**************/
/* Destructor */
/**************/
GaussianEDDAParameter::~GaussianEDDAParameter() {
	if (_tabInvSqrtDetSigma) {
		delete[] _tabInvSqrtDetSigma;
		_tabInvSqrtDetSigma = NULL;
	}

	if (_tabInvSigma) {
		delete[] _tabInvSigma;
		_tabInvSigma = NULL;
	}

	if (_tabSigma) {
		delete[] _tabSigma;
		_tabSigma = NULL;
	}
}

//---------------------
/// Comparison operator
//---------------------
bool GaussianEDDAParameter::operator ==(const GaussianEDDAParameter & param) const {
	if (!GaussianParameter::operator==(param)) return false;
	for (int k = 0; k < _nbCluster; k++) {
		if (_tabSigma[k] != param.getTabSigma()[k]) return false;
	}
	return true;
}

//------------------------
// reset to default values
//------------------------
void GaussianEDDAParameter::reset() {
	for (int k = 0; k < _nbCluster; k++) {
		_tabInvSqrtDetSigma[k] = 0.0;
		*(_tabSigma[k]) = 1.0;
		*(_tabInvSigma[k]) = 1.0;
	}
	GaussianParameter::reset();
}

void GaussianEDDAParameter::getAllPdf(float** tabFik, float* tabProportion)const {

	GaussianData * data = _model->getGaussianData();
	int nbSample = _model->getNbSample();
	float * muk; // to store mu_k
	float InvPiInvDetProp; // 1/(pi)^(d/2) * 1/det(\Sigma_k) ^ (.5) * p_k
	Matrix * SigmakMoins1;
	float ** matrix = data->getYStore(); // to store x_i i=1,...,n

	float ** p_matrix; // pointeur : parcours de y
	float ** p_tabFik_i; // pointeur : parcours de tabFik

	float norme;

	// parcours
	int i; // parcours des individus
	int k; // cluster index

	float * xiMoinsMuk = data->getTmpTabOfSizePbDimension();
	float * xi;
	int p;

	for (k = 0; k < _nbCluster; k++) {
		InvPiInvDetProp = data->getInv2PiPow() * _tabInvSqrtDetSigma[k] * tabProportion[k];
		muk = _tabMean[k];
		SigmakMoins1 = _tabInvSigma[k];

		//computing
		p_matrix = matrix; // pointe sur la premiere ligne des donnees
		p_tabFik_i = tabFik;
		for (i = 0; i < nbSample; i++) {
			xi = *p_matrix;
			for (p = 0; p < _pbDimension; p++) {
				xiMoinsMuk[p] = xi[p] - muk[p];
			}
			// calcul de la norme (x-muk)'Sigma_k^{-1} (x-muk)
			norme = SigmakMoins1->norme(xiMoinsMuk); // virtual !!
			(*p_tabFik_i)[k] = InvPiInvDetProp * exp(-0.5 * norme);

			p_matrix++;
			p_tabFik_i++;
		} // end for i
	} // end for k
}

//-------
// getPdf
//-------
float GaussianEDDAParameter::getPdf(Sample * x, int kCluster)const {

	GaussianData * data = _model->getGaussianData();

	float normPdf = 0.0; // Normal pdf
	float invPi = data->getInv2PiPow(); // Inverse of (2*Pi)^(d/2)

	// Compute (x-m)*inv(S)*(x-m)'
	float * ligne = ((GaussianSample*) x)->getTabValue();

	Matrix * Sigma_kMoins1 = _tabInvSigma[kCluster];

	float norme;
	float * muk = _tabMean[kCluster];
	float * xiMoinsMuk = data->getTmpTabOfSizePbDimension();
	int p;

	for (p = 0; p < _pbDimension; p++) {
		xiMoinsMuk[p] = ligne[p] - muk[p];
	}

	// calcul de la norme (x-muk)'Sigma_k^{-1} (x-muk)
	norme = Sigma_kMoins1->norme(xiMoinsMuk); // virtual

	/* Compute normal density */
	normPdf = invPi * _tabInvSqrtDetSigma[kCluster] * exp(-0.5 * norme);

	return normPdf;
}

float GaussianEDDAParameter::getPdf(int iSample, int kCluster)const {

	GaussianData * data = _model->getGaussianData();
	float * x = (data->getYStore())[iSample];

	Matrix * Sigma_kCluster_moins1 = _tabInvSigma[kCluster];

	float norme;

	float * xiMoinsMuk = data->getTmpTabOfSizePbDimension();
	float * muk = _tabMean[kCluster];

	int p;

	for (p = 0; p < _pbDimension; p++) {
		xiMoinsMuk[p] = x[p] - muk[p];
	}

	// calcul de la norme (x-muk)'Sigma_k^{-1} (x-muk)
	norme = Sigma_kCluster_moins1->norme(xiMoinsMuk); // virtual

	float normPdf;

	// Compute normal density
	normPdf = data->getInv2PiPow() * _tabInvSqrtDetSigma[kCluster] * exp(-0.5 * norme);

	return normPdf;
}

void GaussianEDDAParameter::updateTabInvSigmaAndDet() {
	int k;
	float detSigma;
	for (k = 0; k < _nbCluster; k++) {
		NumericException error = NumericException(minDeterminantSigmaValueError);
		detSigma = _tabSigma[k]->determinant(error);
		_tabSigma[k]->inverse(_tabInvSigma[k]);
		_tabInvSqrtDetSigma[k] = 1.0 / sqrt(detSigma);
	}
}

//---------
// MAP step
//---------
// update _tabWk and _W. 
// No need for this algo but if CV is the criterion, these values
// must be updated
/*void XEMGaussianEDDAParameter::MAPStep(){
  computeTabWkW(); 
}*/


//--------------------
// init USER_PARTITION
//--------------------

void GaussianEDDAParameter::initForInitUSER_PARTITION(int & nbInitializedCluster, 
		bool * tabNotInitializedCluster, Partition * initPartition) 
{
	// init : tabSigma, _tabInvSigma, _tabInvSqrtDetSigma, 
	DiagMatrix * matrixDataVar = new DiagMatrix(_pbDimension, 0.0);
	computeGlobalDiagDataVariance(matrixDataVar);
	// Vaiance matrix initialization to variance matrix of data (Symmetric, Diag or Spherical)
	for (int k = 0; k < _nbCluster; k++) {
		(*_tabSigma[k]) = matrixDataVar;
	}
	updateTabInvSigmaAndDet();
	delete matrixDataVar;

	// _tabMean
	computeTabMeanInitUSER_PARTITION(
			nbInitializedCluster, tabNotInitializedCluster, initPartition);
}

void GaussianEDDAParameter::initUSER(Parameter * iParam) {
	/*
	  updated : _tabMean, _tabWk, _tabSigma, _tabProportion, _tabShape, 
	            _tabOrientation, _tabInvSigma, _tabInvSqrtDetSigma
	 */
	// recuperation
	// we got an XEMGaussianGeneralParameter 
	// because of the implementation in class XEMStrategyType
	// in gaussian cases the init parameters are allways General
	GaussianGeneralParameter * param = (GaussianGeneralParameter *) iParam->getGaussianParameter();
	float ** iTabMean = param->getTabMean();
	float * iTabProportion = param->getTabProportion();
	Matrix ** iTabWk = param->getTabWk();
	Matrix ** iTabSigma = param->getTabSigma();
	int k;
	for (k = 0; k < _nbCluster; k++) {
		recopyTab(iTabMean[k], _tabMean[k], _pbDimension);
		(* _tabWk[k]) = iTabWk[k];
		(* _tabSigma[k]) = iTabSigma[k];
		if (hasFreeProportion(_modelType->_nameModel)) {
			_tabProportion[k] = iTabProportion[k];
		}
		else {
			_tabProportion[k] = 1.0 / _nbCluster;
		}
	}
}

void GaussianEDDAParameter::input(std::ifstream & fi) {

	int j, k;
	float * muk;
	for (k = 0; k < _nbCluster; k++) {
		muk = _tabMean[k];
		// Proportions //
		_tabProportion[k] = getFloatFromStream(fi);

		// Center (mean) //
		for (j = 0; j < _pbDimension; j++)
			muk[j] = getFloatFromStream(fi);
		// Variance matrix //
		_tabSigma[k]->input(fi); // virtual method
	} // end for k
}

//For Heterogeneous model.
void GaussianEDDAParameter::input(std::ifstream & fi, int nbVariables_binary, std::vector< int > nbFactor) {

	int j, k;
	float * muk;
  float garbage;
  int sumNbFactor = 0;
  for (int l = 0; l < nbFactor.size(); l++) sumNbFactor += nbFactor[l]; 
  for (int t = 0; t < _nbCluster * (1 + nbVariables_binary + sumNbFactor); t++) { 
    garbage = getFloatFromStream(fi);
  }
	for (k = 0; k < _nbCluster; k++) {
		muk = _tabMean[k];
		// Proportions //
		_tabProportion[k] = getFloatFromStream(fi);

		// Center (mean) //
		for (j = 0; j < _pbDimension; j++)
			muk[j] = getFloatFromStream(fi);
		// Variance matrix //
		_tabSigma[k]->input(fi); // virtual method
	} // end for k
}

void GaussianEDDAParameter::input(float * proportions
		, float ** means
		, float *** variances
		) {
	for (int k = 0; k < _nbCluster; k++) {

		// Proportions //
		_tabProportion[k] = proportions[k];
		// Center (mean) //
		for (int j = 0; j < _pbDimension; j++) {
			_tabMean[k][j] = means[k][j];
		}
		// Variance matrix //
		_tabSigma[k]->input(variances[k]); // virtual method
	} // end for k
}

void GaussianEDDAParameter::updateForCV(Model * originalModel, CVBlock & CVBlock) {
	GaussianParameter::updateForCV(originalModel, CVBlock);
	computeTabSigma();
}

/*-------------------------------------------------------------------------------------------
	M step
	------
	
	already updated :
	- _model (_tabFik, _tabTik, _tabCik, _tabNk) if Estep is done before
	- _model->_tabCik, _model->_tabNk if USER_PARTITION is done before (Disciminant analysis)
	In all cases, only  _tabCik and _tabNk are needed 
	
	updated in this method :
	- in XEMGaussianParameter::MStep :
			-	_tabProportion
			- _tabMean
			- _tabWk 
			- _W
	- here :
			- _tabSigma
			- _tabInvSigma
			- _tabInvSqrtDetSigma
--------------------------------------------------------------------------------------------*/
void GaussianEDDAParameter::MStep() {
	GaussianParameter::MStep();
	computeTabSigma();
}

//--------------------------------------------
// initialize attributes before an InitRandom  
//--------------------------------------------
void GaussianEDDAParameter::initForInitRANDOM() {
	DiagMatrix * matrixDataVar = new DiagMatrix(_pbDimension, 0.0);
	computeGlobalDiagDataVariance(matrixDataVar);

	// Vaiance matrix initialization to variance matrix of data (Symmetric, Diag or Spherical)
	// init : tabSigma, _tabInvSigma, _tabInvSqrtDetSigma, 
	for (int k = 0; k < _nbCluster; k++) {
		(*_tabSigma[k]) = matrixDataVar;
	}
	updateTabInvSigmaAndDet();

	delete matrixDataVar;
}

//-----------------------------------------
// computeTik when underflow
// -model->_tabSumF[i] pour ith sample = 0
// i : 0 ->_nbSample-1
//-----------------------------------------
void GaussianEDDAParameter::computeTikUnderflow(int i, float ** tabTik) {

	GaussianData * data = _model->getGaussianData();
	int * lnFk = new int [_nbCluster];
	float * lnFkPrim = new float[_nbCluster];
	float * fkPrim = new float[_nbCluster];
	int k, k0;
	float lnFkMax, fkTPrim;
	float detSigma;

	float * ligne = (data->getYStore())[i];
	float norme;
	float ** p_tabMean;
	float * tabTik_i = tabTik[i];

	float * xiMoinsMuk = data->getTmpTabOfSizePbDimension();
	float * muk;

	p_tabMean = _tabMean;
	int p;
	for (k = 0; k < _nbCluster; k++) {
		NumericException error = NumericException(minDeterminantSigmaValueError);
		detSigma = (_tabSigma[k])->determinant(error); // virtual method
		muk = *p_tabMean;
		for (p = 0; p < _pbDimension; p++) {
			xiMoinsMuk[p] = ligne[p] - muk[p];
		}
		norme = (_tabInvSigma[k])->norme(xiMoinsMuk); // virtual method

		// compute lnFik[k]

		lnFk[k] = log(_tabProportion[k]) - data->getHalfPbDimensionLog2Pi() 
				- 0.5 * log(detSigma) - 0.5 * norme;

		p_tabMean++;
	} // end for k

	lnFkMax = lnFk[0];
	for (k = 1; k < _nbCluster; k++) {
		if (lnFk[k] > lnFkMax) {
			lnFkMax = lnFk[k];
		}
	}

	fkTPrim = 0.0;
	for (k = 0; k < _nbCluster; k++) {
		lnFkPrim[k] = lnFk[k] - lnFkMax;
		fkPrim[k] = exp(lnFkPrim[k]);
		fkTPrim += fkPrim[k];
	}

	// compute tabTik
	if (fkTPrim == 0) {
		// reset tabTik
		initToZero(tabTik_i, _nbCluster);
		k0 = GaussianParameter::computeClassAssigment(i);
		tabTik_i[k0] = 1.0;
	}
	else {
		for (k = 0; k < _nbCluster; k++) {
			tabTik_i[k] = fkPrim[k] / fkTPrim;
		}
	}

	delete [] lnFk;
	delete [] lnFkPrim;
	delete [] fkPrim;
}

//debug
void GaussianEDDAParameter::edit() {
	int k;
	for (k = 0; k < _nbCluster; k++) {
		cout << "\tcomponent : " << k << endl;
		cout << "\t\tproportion : " << _tabProportion[k] << endl;
		editTab(_tabMean + k, 1, _pbDimension, cout, " ", "\t\tmean : ");
		cout << "\t\tsigma : " << endl;
		_tabSigma[k]->edit(cout, "\t\t\t");

		cout << "\t\tWk : " << endl;
		_tabWk[k]->edit(cout, "\t\t\t");

		cout << "\t\tinvSigma : " << endl;
		_tabInvSigma[k]->edit(cout, "\t\t\t");

		cout << "\t\ttabInvSqrtDetSigma : " << _tabInvSqrtDetSigma[k] << endl;
	}

	cout << "\tW : " << endl;
	_W->edit(cout, "\t\t");
}

//------
// Edit
//-------
void GaussianEDDAParameter::edit(std::ofstream & oFile, bool text) {
	int k;
	if (text) {
		for (k = 0; k < _nbCluster; k++) {
			oFile << "\t\t\tComponent " << k + 1 << endl;
			oFile << "\t\t\t---------" << endl;
			oFile << "\t\t\tMixing proportion : " << _tabProportion[k] << endl;
			editTab(_tabMean + k, 1, _pbDimension, oFile, " ", "\t\t\tMean : ");

			oFile << "\t\t\tCovariance matrix : " << endl;
			_tabSigma[k]->edit(oFile, "\t\t\t\t");

			oFile << endl;
		}
		oFile << endl;
	}
	else {
		for (k = 0; k < _nbCluster; k++) {
			putFloatInStream(oFile, _tabProportion[k]);

			editTab(_tabMean + k, 1, _pbDimension, oFile, " ", "");
			_tabSigma[k]->edit(oFile, "");
			oFile << endl;

		}
		oFile << endl;

	}
}

void GaussianEDDAParameter::recopy(Parameter * otherParameter) {

	GaussianEDDAParameter * iParameter = 
			(GaussianEDDAParameter *) otherParameter->getGaussianParameter();
	recopyTab(iParameter->getTabMean(), _tabMean, _nbCluster, _pbDimension);
	//  _W->recopy(iParameter->getW());
	(*_W) = iParameter->getW();

	int k;
	Matrix ** iTabSigma = iParameter->getTabSigma();
	Matrix ** iTabInvSigma = iParameter->getTabInvSigma();
	Matrix ** iTabWk = iParameter->getTabWk();
	for (k = 0; k < _nbCluster; k++) {
		// _tabSigma[k]->recopy(iTabSigma[k]);
		(* _tabSigma[k]) = iTabSigma[k];
		// _tabInvSigma[k]->recopy(iTabInvSigma[k]);
		(*_tabInvSigma[k]) = iTabInvSigma[k];
		// _tabWk[k]->recopy(iTabWk[k]);  
		(* _tabWk[k]) = iTabWk[k];
	}
	recopyTab(iParameter->getTabInvSqrtDetSigma(), _tabInvSqrtDetSigma, _nbCluster);
}

}
