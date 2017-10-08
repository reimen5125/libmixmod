/***************************************************************************
                             SRC/mixmod/Matrix/SphericalMatrix.h  description
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
#ifndef XEMSPHERICALMATRIX_H
#define XEMSPHERICALMATRIX_H

#include "mixmod/Matrix/Matrix.h"

namespace XEM {

// pre-declaration
class GeneralMatrix;
class DiagMatrix;

/**
  @brief class XEMSphericalMatrix
  @author F Langrognet & A Echenim
 */

class SphericalMatrix : public Matrix {

public:

	/// Default constructor
	SphericalMatrix();

	/// contructor
	/// default initialisation : Id
	SphericalMatrix(int pbDimension, float initValue = 1.0);

	SphericalMatrix(SphericalMatrix * A);

	/// Desctructor
	virtual ~SphericalMatrix();

	/// compute determinant of spherical matrix
	float determinant(Exception& errorType);

	/// return store of spherical matrix
	float getStore();


	/// inverse spherical matrix
	void inverse(Matrix * & A);

	void compute_product_Lk_Wk(Matrix* Wk, float L);

	/// add a to the value of this
	void addToValue(float a);

	/// compute (x - mean)' this (x - mean) 
	float norme(float * xMoinsMean);

	/// (this) will be A / d
	void equalToMatrixDividedByFloat(Matrix * A, float d);

	/// (this) will be A * d
	void equalToMatrixMultiplyByFloat(Matrix*D, float d);

	/// compute trace of spherical matrix
	float computeTrace();


	/// add :  cik * xMoinsMean * xMoinsMean'  to this
	void add(float * xMoinsMean, float cik);

	/// add : diag( cik * xMoinsMean * xMoinsMean' )  to this
	void addDiag(float * xMoinsMean, float cik);

	/// this  = d * Identity
	void operator=(const float& d);

	/// this = this / (d * Identity)
	void operator/=(const float& d);

	/// this = this *  (d * Identity)
	void operator*=(const float& d);
	/// this = this + matrix
	void operator+=(Matrix* M);
	/// this = matrix
	void operator=(Matrix* M);


	/// read spherical matrix store in file
	void input(std::ifstream & fi);
	virtual void input(float ** variances);

	/// return store of a spherical matrix
	float putSphericalValueInStore(float & store);
	/// add store of a spherical matrix
	float addSphericalValueInStore(float & store);

	float getSphericalStore();

	/// Return store of a diagonal matrix
	float* putDiagonalValueInStore(float * store);
	/// Add store of a diagonal matrix in a diagonal one
	float* addDiagonalValueInStore(float * store);

	float* getDiagonalStore();

	/// Return store of a diagonal matrix
	float* putSymmetricValueInStore(float * store);
	/// Add store of a diagonal matrix in a diagonal one
	float* addSymmetricValueInStore(float * store);

	float* getSymmetricStore();

	/// Return store of a diagonal matrix
	float* putGeneralValueInStore(float * store);
	/// Add store of a diagonal matrix in a diagonal one
	float* addGeneralValueInStore(float * store);

	float* getGeneralStore();


	/// compute general matrix SVD decomposition
	void computeSVD(DiagMatrix* & S, GeneralMatrix* & O);

	void compute_as__multi_O_S_O(float multi, GeneralMatrix* & O, DiagMatrix *& S);
	float trace_this_O_Sm1_O(GeneralMatrix* & O, DiagMatrix* & S);
	float compute_trace_W_C(Matrix * C);
	void computeShape_as__diag_Ot_this_O(DiagMatrix* & Shape, GeneralMatrix* & Ori, float diviseur = 1.0);
	/// gives : det(diag(this))
	float detDiag(Exception& errorType);

	void setSymmetricStore(float * store);
	void setGeneralStore(float * store);
	void setDiagonalStore(float * store);
	void setSphericalStore(float store);
	float** storeToArray() const;

protected:

	float _store;
};

inline float SphericalMatrix::getStore() {
	return _store;
}

inline void SphericalMatrix::setSymmetricStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

inline void SphericalMatrix::setSphericalStore(float store) {
	_store = store;
}

inline void SphericalMatrix::setGeneralStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

inline void SphericalMatrix::setDiagonalStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

}

#endif
