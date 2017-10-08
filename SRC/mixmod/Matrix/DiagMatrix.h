/***************************************************************************
                             SRC/mixmod/Matrix/DiagMatrix.h  description
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
#ifndef XEMDIAGMATRIX_H
#define XEMDIAGMATRIX_H

#include "mixmod/Matrix/Matrix.h"

namespace XEM {

// pre-declaration
class GeneralMatrix;

/**
  @brief class XEMDiagMatrix
  @author F Langrognet & A Echenim
 */

class DiagMatrix : public Matrix {

public:

	/// Default constructor
	DiagMatrix();

	/// contructor : d*Id
	/// default value  = Id
	DiagMatrix(int pbDimension, float d = 1.0);

	DiagMatrix(DiagMatrix * A);

	/// Desctructor
	virtual ~DiagMatrix();

	/// compute determinant of diagonal matrix
	float determinant(Exception& errorType);
	/// return store of diagonal matrix
	float * getStore();
	/// compute inverse of diagonal matrix
	void inverse(Matrix * & A);

	void compute_product_Lk_Wk(Matrix* Wk, float L);

	/// compute (x - mean)' this (x - mean) 
	float norme(float * xMoinsMean);

	/// (this) will be A / d
	void equalToMatrixDividedByFloat(Matrix * A, float d);

	/// this = matrix * d
	void equalToMatrixMultiplyByFloat(Matrix*D, float d);

	///compute singular vector decomposition
	void computeSVD(DiagMatrix* & S, GeneralMatrix* & O);

	/// compute trace of general matrix
	float computeTrace();

	/// add :  cik * xMoinsMean * xMoinsMean'  to this
	void add(float * xMoinsMean, float cik);

	// add : diag( cik * xMoinsMean * xMoinsMean' )  to this
	//void addDiag(float * xMoinsMean, float cik);

	/// set the value of (d x Identity) to this  
	void operator=(const float& d);
	/// this = this / (d * Identity)
	void operator/=(const float& d);
	/// this = this * (d * Identity)
	void operator*=(const float& d);
	/// this = this + matrix
	void operator+=(Matrix* M);
	/// this = matrix
	void operator=(Matrix* M);

	/// Return store of a spherical matrix in a diagonal one
	float putSphericalValueInStore(float & store);
	/// Add store of a spherical matrix in a diagonal one
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

	/// read general matrix in an input file
	void input(std::ifstream & fi);
	virtual void input(float ** variances);

	///set store
	void setSymmetricStore(float * store);
	void setGeneralStore(float * store);
	void setDiagonalStore(float * store);
	void setSphericalStore(float store);
	float** storeToArray() const;

	/// gives : det(diag(this))
	float detDiag(Exception& errorType);

	void compute_as__multi_O_S_O(float multi, GeneralMatrix* & O, DiagMatrix *& S);
	float trace_this_O_Sm1_O(GeneralMatrix* & O, DiagMatrix* & S);
	float compute_trace_W_C(Matrix * C);
	void computeShape_as__diag_Ot_this_O(DiagMatrix* & Shape, GeneralMatrix* & Ori, float diviseur = 1.0);

	///sort diagonal matrix in decreasing order
	void sortDiagMatrix();
	
protected:

	float * _store;

};

inline float * DiagMatrix::getStore() {
	return _store;
}

inline void DiagMatrix::setSymmetricStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

inline void DiagMatrix::setGeneralStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

inline void DiagMatrix::setDiagonalStore(float * store) {
	//_store = store;
	recopyTab(store, _store, _s_pbDimension);
}

inline void DiagMatrix::setSphericalStore(float store) {
	THROW(OtherException, wrongMatrixType);
}

}

#endif
