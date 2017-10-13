/***************************************************************************
                             SRC/mixmod/Matrix/SphericalMatrix.cpp  description
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

#include "mixmod/Matrix/SphericalMatrix.h"
#include "mixmod/Matrix/GeneralMatrix.h"
#include "mixmod/Matrix/DiagMatrix.h"

namespace XEM {

//------------
// Constructor
//------------
SphericalMatrix::SphericalMatrix() {
	THROW(OtherException, wrongConstructorType);
}

SphericalMatrix::SphericalMatrix(int pbDimension, float initValue) : Matrix(pbDimension) {
	_store = initValue;
}

SphericalMatrix::SphericalMatrix(SphericalMatrix * A) : Matrix(A) {
	_store = A->getStore();
}

//----------
//Destructor
//----------
SphericalMatrix::~SphericalMatrix() {
}

float SphericalMatrix::determinant(Exception& errorType) {
	float det;
#ifdef __APPLE__
	det = powf(_store, (int) _s_pbDimension);
#else
	det = powf(_store, (float) _s_pbDimension);
#endif
	if (det < minDeterminantValue) {
		throw NumericException(dynamic_cast<NumericException&> (errorType));
	}
	return det;
}

float* SphericalMatrix::getDiagonalStore() {
	THROW(OtherException, wrongMatrixType);
}

float* SphericalMatrix::getSymmetricStore() {
	THROW(OtherException, wrongMatrixType);
}

float* SphericalMatrix::getGeneralStore() {
	THROW(OtherException, wrongMatrixType);
}

float SphericalMatrix::getSphericalStore() {
	return (_store);
}

void SphericalMatrix::computeSVD(DiagMatrix* & S, GeneralMatrix* & O) {
	THROW(OtherException, nonImplementedMethod);
}

void SphericalMatrix::equalToMatrixMultiplyByFloat(Matrix* D, float d) {
	//THROW(XEMInputException,nonImplementedMethod);
	float store_D = D->putSphericalValueInStore(_store);
	_store = store_D*d;
}

void SphericalMatrix::inverse(Matrix * & Inv) {

	/*1.
	  if (Inv == NULL){
		Inv = new XEMSphericalMatrix(_s_pbDimension);
	  }
	  else{
		if (typeof(Inv) != XEMSphericalMatrix){
		  delete Inv;
		  Inv = new XEMSphericalMatrix(_s_pbDimension);
		}
	  }
	 */

	if (Inv == NULL) {
		Inv = new SphericalMatrix(_s_pbDimension);
	}

	// Inv = new XEMSphericalMatrix(_s_pbDimension);
	float store_Inv = 1.0 / _store; // = A->getSphericalStore();
	Inv->setSphericalStore(store_Inv); // virtual
}

void SphericalMatrix::compute_product_Lk_Wk(Matrix* Wk, float L) {
	THROW(OtherException, nonImplementedMethod);
}

float SphericalMatrix::computeTrace() {
	float trace = _s_pbDimension * _store;
	return trace;
}

void SphericalMatrix::addToValue(float a) {
	_store += a;
}

float SphericalMatrix::norme(float * xMoinsMean) {
	int p;
	float termesDiag = 0.0;
	float xMoinsMean_p;

	for (p = 0; p < _s_pbDimension; p++) {
		xMoinsMean_p = xMoinsMean[p];
		termesDiag += xMoinsMean_p * xMoinsMean_p;
	}
	termesDiag *= _store;
	return termesDiag;
}

// (this) will be A / d
void SphericalMatrix::equalToMatrixDividedByFloat(Matrix * A, float d) {
	_store = (A->getSphericalStore()) / d;
}

void SphericalMatrix::compute_as__multi_O_S_O(float multi, GeneralMatrix* & O, DiagMatrix *& S) {
	THROW(OtherException, nonImplementedMethod);
}

float SphericalMatrix::trace_this_O_Sm1_O(GeneralMatrix* & O, DiagMatrix* & S) {
	THROW(OtherException, nonImplementedMethod);
}

float SphericalMatrix::compute_trace_W_C(Matrix * C) {
	THROW(OtherException, nonImplementedMethod);
}

void SphericalMatrix::computeShape_as__diag_Ot_this_O(
		DiagMatrix* & Shape, GeneralMatrix* & Ori, float diviseur) 
{
	THROW(OtherException, nonImplementedMethod);
}

// add :  cik * xMoinsMean * xMoinsMean'  to this
void SphericalMatrix::add(float * xMoinsMean, float cik) {

	int p;
	float xMoinsMean_p, tmp;

	tmp = 0.0;
	for (p = 0; p < _s_pbDimension; p++) {
		xMoinsMean_p = xMoinsMean[p];
		tmp += xMoinsMean_p * xMoinsMean_p;
	}//end for p

	_store += tmp / _s_pbDimension * cik;
}

// add : diag( cik * xMoinsMean * xMoinsMean' )  to this
void SphericalMatrix::addDiag(float * xMoinsMean, float cik) {

	int p;
	float xMoinsMean_p, tmp;

	tmp = 0.0;
	for (p = 0; p < _s_pbDimension; p++) {
		xMoinsMean_p = xMoinsMean[p];
		tmp += xMoinsMean_p * xMoinsMean_p;
	}//end for p

	_store += tmp / _s_pbDimension * cik;
}

// set the value of (d x Identity) to this  
void SphericalMatrix::operator=(const float& d) {
	_store = d;
}

// divide each element by d
void SphericalMatrix::operator/=(const float& d) {
	_store /= d;
}

// multiply each element by d
void SphericalMatrix::operator*=(const float& d) {
	_store *= d;
}

//add M to this
void SphericalMatrix::operator+=(Matrix* M) {
	M->addSphericalValueInStore(_store);
}

void SphericalMatrix::operator=(Matrix* M) {
	M->putSphericalValueInStore(_store);
}

float SphericalMatrix::putSphericalValueInStore(float & store) {
	store = _store;
	return (store);
}

float SphericalMatrix::addSphericalValueInStore(float & store) {
	store += _store;
	return (store);
}

float* SphericalMatrix::putDiagonalValueInStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

float* SphericalMatrix::addDiagonalValueInStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

float* SphericalMatrix::putSymmetricValueInStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

float* SphericalMatrix::addSymmetricValueInStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

float* SphericalMatrix::putGeneralValueInStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

float* SphericalMatrix::addGeneralValueInStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

void SphericalMatrix::input(std::ifstream & fi) {
	int p, q;
	float garbage;

	for (p = 0; p < _s_pbDimension; p++) {
		// useless because all are 0
		for (q = 0; q < _s_pbDimension; q++) {
			if (p == 0 && q == 0)
				_store = getFloatFromStream(fi);
      else {
				// we don't use the value (garbage)
				getFloatFromStream(fi);
			}
		}
	}
}

void SphericalMatrix::input(float ** variances) {
	int p, q;

	for (p = 0; p < _s_pbDimension; p++) {
		// useless because all are 0
		for (q = 0; q < _s_pbDimension; q++) {
			if (p == 0 && q == 0) {
				_store = variances[p][q];
			}
		}
	}
}

float SphericalMatrix::detDiag(Exception& errorType) {
	return determinant(errorType);
}

float** SphericalMatrix::storeToArray() const {

	int i, j;
	float** newStore = new float*[_s_pbDimension];
	for (i = 0; i < _s_pbDimension; ++i) {
		newStore[i] = new float[_s_pbDimension];

		for (j = 0; j < _s_pbDimension; ++j) {
			if (i == j) {
				newStore[i][j] = _store;
			}
			else {
				newStore[i][j] = 0;
			}
		}
	}

	return newStore;
}

}
