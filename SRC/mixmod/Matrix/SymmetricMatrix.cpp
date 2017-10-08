/***************************************************************************
                             SRC/mixmod/Matrix/SymmetricMatrix.cpp  description
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

#include "mixmod/Matrix/SymmetricMatrix.h"
#include "mixmod/Matrix/DiagMatrix.h"
#include "mixmod/Matrix/GeneralMatrix.h"

namespace XEM {

//------------
// Constructor
//------------
SymmetricMatrix::SymmetricMatrix() {
	_value = NULL;
	_store = NULL;
	THROW(OtherException, wrongConstructorType);
}

SymmetricMatrix::SymmetricMatrix(int pbDimension, float d) : Matrix(pbDimension) {
	_value = new MATH::SymmetricMatrix(_s_pbDimension);
	_store = _value->Store();

	_s_storeDim = _s_pbDimension * (_s_pbDimension + 1) / 2;
	(*this) = d;
}

// copy constructor
SymmetricMatrix::SymmetricMatrix(SymmetricMatrix * A) : Matrix(A) {
	_value = new MATH::SymmetricMatrix(_s_pbDimension);
	_store = _value->Store();

	_s_storeDim = _s_pbDimension * (_s_pbDimension + 1) / 2;

	recopyTab(A->getStore(), _store, _s_storeDim);
}

//----------
//Destructor
//----------
SymmetricMatrix::~SymmetricMatrix() {
	if (_value) {
		delete _value;
	}
	_value = NULL;
}

float SymmetricMatrix::determinant(Exception& errorType) {
	float det = 0;
	try {
		det = _value->LogDeterminant();
	}
	catch (...) {
		throw errorType;
	}
	if (det < minDeterminantValue) {
		throw NumericException(dynamic_cast<NumericException&> (errorType));
	}
	return det;
}

void SymmetricMatrix::equalToMatrixMultiplyByFloat(Matrix* D, float d) {
	THROW(OtherException, nonImplementedMethod);
}

float* SymmetricMatrix::getDiagonalStore() {
	THROW(OtherException, wrongMatrixType);
}

float* SymmetricMatrix::getSymmetricStore() {
	return (_store);
}

float* SymmetricMatrix::getGeneralStore() {
	return (_store);
}

float SymmetricMatrix::getSphericalStore() {
	THROW(OtherException, wrongMatrixType);
}

void SymmetricMatrix::compute_M_tM(float* V, int l) {
	int indice1 = l - 1, indice2;
	int indiceStoreGammak = _s_storeDim - 1;
	int dim = l / _s_pbDimension;
	while (indice1 > 0) {
		for (int j = 0; j < dim; j++) {
			_store[indiceStoreGammak] += V[(indice1 - j)] * V[(indice1 - j)];
		}
		indiceStoreGammak -= 1;
		indice2 = indice1 - dim;
		while (indice2 > 0) {
			for (int j = 0; j < dim; j++) {
				_store[indiceStoreGammak] += V[(indice1 - j)] * V[(indice2 - j)];
			}
			indice2 -= dim;
			indiceStoreGammak -= 1;
		}
		indice1 -= dim;
	}
}

void SymmetricMatrix::compute_product_Lk_Wk(Matrix* Wk, float L) {
	float * Wk_store;
	Wk_store = Wk->getSymmetricStore();
	for (int p = 0; p < _s_storeDim; p++) {
		_store[p] += Wk_store[p] / L;
	}
}

float SymmetricMatrix::compute_trace_W_C(Matrix * C) {
	float tabLambdak_k = 0.0;
	float termesHorsDiag;
	int p, q, r;
	float * C_store = C->getSymmetricStore();
	termesHorsDiag = 0.0;
	for (p = 0, r = 0; p < _s_pbDimension; p++, r++) {
		for (q = 0; q < p; q++, r++) {
			termesHorsDiag += _store[r] * C_store[r];
		}
		tabLambdak_k += _store[r] * C_store[r];
	}
	tabLambdak_k += 2.0 * termesHorsDiag;
	return tabLambdak_k;
}

void SymmetricMatrix::inverse(Matrix * & Inv) {
	//cout<<"Inv Symm :  "<<Inv<<endl;
	if (Inv == NULL) {
		Inv = new SymmetricMatrix(_s_pbDimension);
	}

	MATH::SymmetricMatrix* value_Inv = _value->Inverse();

	Inv->setSymmetricStore(value_Inv->Store());
	//cout<<"Inv Symm :  "<<Inv<<endl;
	delete value_Inv;
}

float SymmetricMatrix::norme(float * xMoinsMean) {
	int p, q, r;
	float termesHorsDiag = 0.0;
	float termesDiag = 0.0;
	float xMoinsMean_p;

	for (p = 0, r = 0; p < _s_pbDimension; p++, r++) {
		xMoinsMean_p = xMoinsMean[p];
		for (q = 0; q < p; q++, r++) {
			termesHorsDiag += xMoinsMean_p * xMoinsMean[q] * _store[r];
		}
		termesDiag += xMoinsMean_p * xMoinsMean_p * _store[r];
	}
	termesDiag += 2.0 * termesHorsDiag;
	return termesDiag;
}

float SymmetricMatrix::putSphericalValueInStore(float & store) {
	int p, r;
	int increment = 2;
	store = 0.0;

	for (p = 0, r = 0; p < _s_pbDimension; p++) {
		store += _store[r];
		r += increment;
		increment++;
	}
	store /= _s_pbDimension;
	return (store);

}

float SymmetricMatrix::addSphericalValueInStore(float & store) {
	int p, r;
	int increment = 2;
	for (p = 0, r = 0; p < _s_pbDimension; p++) {
		store += _store[r];
		r += increment;
		increment++;
	}
	store /= _s_pbDimension;
	return (store);
}

float* SymmetricMatrix::putDiagonalValueInStore(float * store) {
	int p, r;
	int increment = 2;
	for (p = 0, r = 0; p < _s_pbDimension; p++) {
		store[p] = _store[r];
		r += increment;
		increment++;
	}
	return (store);
}

float* SymmetricMatrix::addDiagonalValueInStore(float * store) {
	int p, r;
	int increment = 2;
	for (p = 0, r = 0; p < _s_pbDimension; p++) {
		store[p] += _store[r];
		r += increment;
		increment++;
	}
	return (store);
}

float* SymmetricMatrix::putSymmetricValueInStore(float * store) {
	for (int p = 0; p < _s_storeDim; p++) {
		store[p] = _store[p];
	}
	return (store);
}

float* SymmetricMatrix::addSymmetricValueInStore(float * store) {
	for (int p = 0; p < _s_storeDim; p++) {
		store[p] += _store[p];
	}
	return (store);
}

float* SymmetricMatrix::putGeneralValueInStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

float* SymmetricMatrix::addGeneralValueInStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

// (this) will be A / d
void SymmetricMatrix::equalToMatrixDividedByFloat(Matrix * A, float d) {
	A->putSymmetricValueInStore(_store);

	int p;
	for (p = 0; p < _s_storeDim; p++) {
		_store[p] /= d;
	}
}

// add :  cik * xMoinsMean * xMoinsMean'  to this
void SymmetricMatrix::add(float * xMoinsMean, float cik) {

	int p, q, r;
	float xMoinsMean_p;

	for (p = 0, r = 0; p < _s_pbDimension; p++, r++) {
		xMoinsMean_p = xMoinsMean[p];
		for (q = 0; q < p; q++, r++) {
			_store[r] += cik * xMoinsMean_p * xMoinsMean[q];
		} // end for q
		_store[r] += cik * xMoinsMean_p * xMoinsMean_p;
	}//end for p
}

// add : diag( cik * xMoinsMean * xMoinsMean' )  to this
/*void SymmetricMatrix::addDiag(float * xMoinsMean, float cik){
  int p,q,r;
  float xMoinsMean_p;

  for(p=0,r=0 ; p<_s_pbDimension ; p++,r++){
	xMoinsMean_p = xMoinsMean[p];
	for(q=0 ; q<p ; q++,r++) ;
	_store[r]  +=  cik * xMoinsMean_p * xMoinsMean_p;
  }  
}*/

// set the value of (d x Identity) to this  
void SymmetricMatrix::operator=(const float& d) {

	int p, q, r;
	for (p = 0, r = 0; p < _s_pbDimension; p++, r++) {
		for (q = 0; q < p; q++, r++) {
			_store[r] = 0.0;
		}
		_store[r] = d;
	}
}

// divide each element by d
void SymmetricMatrix::operator/=(const float& d) {
	int p;
	for (p = 0; p < _s_storeDim; p++) {
		_store[p] /= d;
	}
}

// multiply each element by d
void SymmetricMatrix::operator*=(const float& d) {
	int p;
	for (p = 0; p < _s_storeDim; p++) {
		_store[p] *= d;
	}
}

void SymmetricMatrix::operator=(Matrix* M) {
	M->putSymmetricValueInStore(_store);
}

//add M to this
void SymmetricMatrix::operator+=(Matrix* M) {
	M->addSymmetricValueInStore(_store);
}

// compute Shape as diag(Ot . this . O ) / diviseur
void SymmetricMatrix::computeShape_as__diag_Ot_this_O(
		DiagMatrix* & Shape, GeneralMatrix* & Ori, float diviseur) 
{
	int i_index, j_index;

	int p, q, r, j;
	float * O_store = Ori->getStore();
	float * Shape_store = Shape->getStore();

	float termesDiag, termesHorsDiag;
	float tmp;

	for (j = 0; j < _s_pbDimension; j++) {
		// computation of the [j,j] term of the diagonal

		//-- reset
		termesDiag = 0.0;
		termesHorsDiag = 0.0;

		i_index = j;
		for (p = 0, r = 0; p < _s_pbDimension; p++, r++, i_index += _s_pbDimension) {
			j_index = j;
			for (q = 0; q < p; q++, r++, j_index += _s_pbDimension) {
				termesHorsDiag += O_store[i_index] * O_store[j_index] * _store[r];
			}
			tmp = O_store[i_index];
			termesDiag += tmp * tmp * _store[r];
		}

		termesDiag += 2.0 * termesHorsDiag;
		termesDiag /= diviseur;
		Shape_store[j] = termesDiag;
	}
}

// compute this as : multi * (O * S * O' )
void SymmetricMatrix::compute_as__multi_O_S_O(float multi, GeneralMatrix* & O, DiagMatrix* & S) {

	int i_index = 0;
	int j_index;
	int p, q, r, l;

	float * O_store = O->getStore();
	float * S_store = S->getStore();
	float tmp;

	for (p = 0, r = 0; p < _s_pbDimension; p++, i_index += _s_pbDimension) {
		j_index = 0;
		for (q = 0; q <= p; q++, r++, j_index += _s_pbDimension) {
			// compute this[i,j] = \multi * sum_{l} ( O[i,l] * 0[j,l] * S|l] ) 
			tmp = 0.0;
			for (l = 0; l < _s_pbDimension; l++) {
				tmp += O_store[i_index + l] * O_store[j_index + l] * S_store[l];
			}
			tmp *= multi;
			_store[r] = tmp;
		}
	}
}

// compute this as : (O * S * O' )
void SymmetricMatrix::compute_as_O_S_O(GeneralMatrix* & O, float* & S_store) {

	int i_index = 0;
	int j_index;
	int p, q, r, l;

	for (int i = 0; i < _s_storeDim; i++) {
		_store[i] = 0;
	}

	float * O_store = O->getStore();
	float tmp;
	for (p = 0, r = 0; p < _s_pbDimension; p++, i_index += _s_pbDimension) {
		j_index = 0;
		for (q = 0; q <= p; q++, r++, j_index += _s_pbDimension) {
			tmp = 0.0;
			for (l = 0; l < _s_pbDimension; l++) {
				tmp += O_store[i_index + l] * O_store[j_index + l] * S_store[l];
			}
			_store[r] = tmp;
		}
	}
}

//compute trace of a symmetric matrix
float SymmetricMatrix::computeTrace() {
	int i;
	int indice = 0;
	float trace = 0.0;


	i = 0;
	while (indice < _s_storeDim) {
		trace += _store[indice];
		i++;
		indice += i + 1;
	}
	return trace;
}

void SymmetricMatrix::computeSVD(DiagMatrix* & S, GeneralMatrix* & O) {
	int dim = O->getPbDimension();
	MATH::DiagonalMatrix * tabShape_k = new MATH::DiagonalMatrix(dim);
	MATH::Matrix * tabOrientation_k = new MATH::Matrix(dim, dim);
	_value->computeSVD(tabShape_k, tabOrientation_k);

	float * storeS = S->getStore();
	float * storeO = O->getStore();

	float * storeTabShape_k = tabShape_k->Store();
	float * storeTabOrientation_k = (*tabOrientation_k).Store();

	recopyTab(storeTabShape_k, storeS, dim);
	recopyTab(storeTabOrientation_k, storeO, dim * dim);

	delete tabShape_k;
	delete tabOrientation_k;
}

void SymmetricMatrix::compute_as_M_tM(GeneralMatrix* M, int d) {

	int indiceStoreM1 = 0, indiceStoreM2;
	int indice = 0;
	int k1 = 0, k2;
	int DimStoreM = _s_pbDimension*_s_pbDimension;
	float * storeM = M->getStore();

	for (int i = 0; i < _s_storeDim; i++) {
		_store[i] = 0;
	}

	while (indiceStoreM1 < DimStoreM) {
		k2 = k1;
		indiceStoreM2 = indiceStoreM1;
		while (indiceStoreM2 < DimStoreM) {

			for (int j = 0; j < d; j++) {
				// attention vecteur contenant la matrice triangulaire supérieure
				_store[indice] += storeM[(indiceStoreM1 + j)] * storeM[(indiceStoreM2 + j)];
			}
			indiceStoreM2 = (k2 + 1) * _s_pbDimension;
			k2 += 1;
			indice += 1;
		}
		indiceStoreM1 = (k1 + 1) * _s_pbDimension;
		k1 += 1;
	}
}

void SymmetricMatrix::compute_as_M_V(SymmetricMatrix* M, float * V) {

	for (int i = 0; i < _s_pbDimension; i++) {
		_store[i] = 0;
	}
	int indiceV = 0, k = 0, indice = 0;
	int indiceM = 0;
	float* storeM = M->getStore();

	while (indice < _s_pbDimension) {
		for (int i = 0; i < (_s_pbDimension - k); i++) {
			_store[indice] += V[indiceV + i] * storeM[indiceM + i];
		}

		for (int j = 1; j < (_s_pbDimension - k); j++) {
			_store[indice + j] += V[indiceV] * storeM[indiceM + j];
		}
		indiceM += (_s_pbDimension - k);
		k += 1;
		indiceV += 1;
		indice += 1;
	}
}

// compute M as : M = ( O * S^{-1} * O' ) * this
void SymmetricMatrix::compute_M_as__O_Sinverse_Ot_this(
		GeneralMatrix & M, GeneralMatrix* & O, DiagMatrix* & S) 
{
	float * M_store = M.getStore();
	float * O_store = O->getStore();
	float * S_store = S->getStore();

	int i, j, l, p, r;
	int O1_index = 0;
	int O2_index;
	int r_decalage;
	float tmp, omega;
	int fillindex = 0;

	for (i = 0; i < _s_pbDimension; i++) {
		for (j = 0; j < _s_pbDimension; j++, fillindex++) {
			// filling tabMtmpk_store[i,j]

			tmp = 0.0;
			r_decalage = _s_pbDimension - j + 1;

			r = j;
			p = 0;
			O2_index = 0;


			while (p < j) {
				omega = 0.0;
				for (l = 0; l < _s_pbDimension; l++) {

					omega += O_store[O1_index + l] * O_store[O2_index + l] / S_store[l];
				}
				tmp += omega * _store[r];
				r += r_decalage;
				r_decalage--;
				p++;
				O2_index += _s_pbDimension;
			}

			while (p < _s_pbDimension) {
				omega = 0.0;

				for (l = 0; l < _s_pbDimension; l++) {
					omega += O_store[O1_index + l] * O_store[O2_index + l] / S_store[l];
				}
				tmp += omega * _store[r];
				r++;
				O2_index += _s_pbDimension;
				p++;
			}


			M_store[fillindex] = tmp;
		}
		O1_index += _s_pbDimension;
	}
}

void SymmetricMatrix::input(std::ifstream & fi) {
	int i, j, r = 0;
	float garbage;

	for (i = 0; i < _s_pbDimension; i++) {
    for (j = 0; j < i + 1; j++) {
			_store[r] = getFloatFromStream(fi);
      r++;
    }
    for (j = i + 1; j < _s_pbDimension; j++) {
			// we don't need values (all are 0 ?)
			getFloatFromStream(fi);
		}
	}
}

void SymmetricMatrix::input(float ** variances) {
	int i, j, r = 0;

	for (i = 0; i < _s_pbDimension; i++) {
		for (j = 0; j < i + 1; j++) {
			_store[r] = variances[i][j];
			r++;
		}
		for (j = i + 1; j < _s_pbDimension; j++) {
		}
	}
}

// gives : det(diag(this))
float SymmetricMatrix::detDiag(Exception& errorType) {
	int p, q, r;
	float det = 1.0;

	for (p = 0, r = 0; p < _s_pbDimension; p++, r++) {
		for (q = 0; q < p; q++, r++);
		det *= _store[r];
	}
	if (det < minDeterminantValue)
		throw errorType;
	return det;
}

// trace( this * O * S^{-1} * O' )
float SymmetricMatrix::trace_this_O_Sm1_O(GeneralMatrix* & O, DiagMatrix* & S) {
	float * O_store = O->getStore();
	float * S_store = S->getStore();

	float trace = 0.0;
	float termesHorsDiag = 0.0;
	float tmp, tmp2;

	int i_index = 0;
	int j_index;
	int p, q, r, l;

	for (p = 0, r = 0; p < _s_pbDimension; p++, r++, i_index += _s_pbDimension) {
		j_index = 0;
		for (q = 0; q < p; q++, r++, j_index += _s_pbDimension) {
			tmp = 0.0;

			for (l = 0; l < _s_pbDimension; l++) {
				tmp += O_store[i_index + l] * O_store[j_index + l] / S_store[l];
			}
			tmp *= _store[r];
			termesHorsDiag += tmp;
		}
		tmp = 0.0;
		for (l = 0; l < _s_pbDimension; l++) {
			tmp2 = O_store[i_index + l];
			tmp += tmp2 * tmp2 / S_store[l];
		}
		tmp *= _store[r];
		trace += tmp;
	}
	trace += 2.0 * termesHorsDiag;

	return trace;
}

float** SymmetricMatrix::storeToArray() const {

	int i, j, k = (_s_storeDim - 1);
	float** newStore = new float*[_s_pbDimension];

	for (i = 0; i < _s_pbDimension; ++i) {
		newStore[i] = new float[_s_pbDimension];
	}
	for (i = (_s_pbDimension - 1); i>-1; --i) {
		newStore[i][i] = _store[k];
		k--;
		for (j = (i - 1); j>-1; --j) {
			newStore[i][j] = _store[k];
			newStore[j][i] = _store[k];
			k--;
		}

	}

	return newStore;
}

}
