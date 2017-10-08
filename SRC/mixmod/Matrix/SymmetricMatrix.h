/***************************************************************************
                             SRC/mixmod/Matrix/SymmetricMatrix.h  description
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
#ifndef XEMSYMMETRICMATRIX_H
#define XEMSYMMETRICMATRIX_H

#include "mixmod/Matrix/Matrix.h"

namespace XEM {

/**
  @brief class GeneralMatrix
  @author F Langrognet & A Echenim
 */

class DiagMatrix;

class SymmetricMatrix : public Matrix {

public:

	/// Default constructor
	SymmetricMatrix();

	/// constructor : d*Id
	/// default value : Id
	SymmetricMatrix(int pbDimension, float d = 1.0);

	SymmetricMatrix(SymmetricMatrix * A);

	/// Destructor
	virtual ~SymmetricMatrix();

	/// Compute determinant of symmetric matrix
	float determinant(Exception& errorType);

	/// Return store of symmetric matrix
	float * getStore();

	/// Return newmat symmetric matrix
	MATH::SymmetricMatrix * getValue();

	/// Return dimension of store
	int getStoreDim();

	/// Inverse symmetric matrix
	void inverse(Matrix * & A);

	/// compute (x - mean)' this (x - mean) 
	float norme(float * xMoinsMean);

	/// compute : this =  A / d
	void equalToMatrixDividedByFloat(Matrix * A, float d);

	/// compute : this =  A * d
	void equalToMatrixMultiplyByFloat(Matrix*D, float d);

	/// add :  cik * xMoinsMean * xMoinsMean'  to this
	void add(float * xMoinsMean, float cik);

	// add : diag( cik * xMoinsMean * xMoinsMean' )  to this
	//void addDiag(float * xMoinsMean, float cik);

	/// Return store of a spherical matrix in a symmetric one
	float putSphericalValueInStore(float & store);
	/// Add store of a spherical matrix in a symmetric one
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


	/// this =  (d * Identity)
	void operator=(const float& d);
	/// this =  this /  (Identity * d)
	void operator/=(const float& d);
	/// this =  this * ( Identity * d)
	void operator*=(const float& d);
	/// this =  this + matrix
	void operator+=(Matrix* M);
	/// this = matrix
	void operator=(Matrix* M);

	
	/// read symmetric matrix store in file
	void input(std::ifstream & fi);
	virtual void input(float ** variances);

	/* ///compute SVD decomposition for a symmetric matrix
	 void computeSVD(XEMDiagMatrix* & S, XEMGeneralMatrix* & O);*/

	/// compute Shape as diag(Ot . this . O ) / diviseur
	void computeShape_as__diag_Ot_this_O(DiagMatrix* & Shape, GeneralMatrix* & Ori, float diviseur = 1.0);
	//  float trace_this_O_Sm1_O(XEMGeneralMatrix* & O, XEMDiagMatrix* & S);

	/// compute this as : multi * (O * S * O' )
	void compute_as__multi_O_S_O(float multi, GeneralMatrix* & O, DiagMatrix* & S);

	/// compute this as O*S*O'
	void compute_as_O_S_O(GeneralMatrix* & O, float* & S_store);
	/// compute trace of this
	float computeTrace();

	/// compute this as M * M'
	void compute_as_M_tM(GeneralMatrix* M, int d);

	/// compute this as matrix * vector
	void compute_as_M_V(SymmetricMatrix* M, float * V);

	/// compute this as float * matrix
	void compute_product_Lk_Wk(Matrix* Wk, float L);

	/// copute trace of W * C
	float compute_trace_W_C(Matrix * C);

	/// compute M as : M = ( O * S^{-1} * O' ) * this
	void compute_M_as__O_Sinverse_Ot_this(GeneralMatrix & M, GeneralMatrix* & O, DiagMatrix* & S);

	/// compute this as vector * vector'
	void compute_M_tM(float* V, int l);

	/// gives : det(diag(this))
	float detDiag(Exception& errorType);

	/// trace( this * O * S^{-1} * O' )
	float trace_this_O_Sm1_O(GeneralMatrix* & O, DiagMatrix* & S);

	///set store

	void setSymmetricStore(float * store);
	void setGeneralStore(float * store);
	void setDiagonalStore(float * store);
	void setSphericalStore(float store);
	float** storeToArray() const;

protected:

	/// Symmetric matrix as in mathematical library (if specified)
	MATH::SymmetricMatrix * _value;

	/// store of matrix
	float * _store;
	
	/// dimension of store
	int _s_storeDim;
};

// TODO static :
// int XEMGeneralMatrix::_s_storeDim = 0;

inline float * SymmetricMatrix::getStore() {
	return _store;
}

inline MATH::SymmetricMatrix * SymmetricMatrix::getValue() {
	return _value;
}

inline int SymmetricMatrix::getStoreDim() {
	return _s_storeDim;
}

inline void SymmetricMatrix::setSymmetricStore(float * store) {
	// _store = store;
	recopyTab(store, _store, _s_storeDim);
}

inline void SymmetricMatrix::setSphericalStore(float store) {
	THROW(OtherException, wrongMatrixType);
}

inline void SymmetricMatrix::setGeneralStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

inline void SymmetricMatrix::setDiagonalStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

/* TODO static
inline void XEMGeneralMatrix::initiate(){
  _s_storeDim = _s_pbDimension * _s_pbDimension / 2; 
} 
 */

}

#endif
