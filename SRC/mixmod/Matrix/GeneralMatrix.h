/***************************************************************************
                             SRC/mixmod/Matrix/GeneralMatrix.h  description
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
#ifndef XEMGENERALMATRIX_H
#define XEMGENERALMATRIX_H

#include "mixmod/Matrix/Matrix.h"

namespace XEM {

// pre-declaration
class DiagMatrix;

/**
  @brief class GeneralMatrix
  @author F Langrognet & A Echenim
 */

class GeneralMatrix : public Matrix {

public:

	/// Default constructor
	GeneralMatrix();

	/// constructor : d*Id
	///default value = Id
	GeneralMatrix(int pbDimension, float d = 1.0);

	GeneralMatrix(GeneralMatrix * A);

	/// Destructor
	virtual ~GeneralMatrix();

	/// compute determinant of general matrix
	float determinant(Exception& errorType);

	/// return store of general matrix
	float * getStore();

	/// return newmat general matrix
	MATH::Matrix * getValue();

	/// return dimension of store
	int getStoreDim();

	/// inverse general matrix
	void inverse(Matrix * & A);

	void compute_product_Lk_Wk(Matrix* Wk, float L);

	/// compute (x - mean)' this (x - mean) 
	float norme(float * xMoinsMean);

	/// this =  A / d
	void equalToMatrixDividedByFloat(Matrix * A, float d);
	/// this =   A * d
	void equalToMatrixMultiplyByFloat(Matrix*D, float d);


	/// add :  cik * xMoinsMean * xMoinsMean'  to this
	void add(float * xMoinsMean, float cik);

	// add : diag( cik * xMoinsMean * xMoinsMean' )  to this
	//void addDiag(float * xMoinsMean, float cik);

	/// return store of a spherical matrix in a general one
	float putSphericalValueInStore(float & store);
	/// add store of a spherical matrix in a general one
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

	/// this =  (d x Identity)
	void operator=(const float& d);

	/// this =  this / (d x Identity)
	void operator/=(const float& d);

	/// this =  this * (d x Identity)
	void operator*=(const float& d);

	/// this =  this + matrix
	void operator+=(Matrix* M);

	/// this =  matrix
	void operator=(Matrix* M);


	/// edit general matrix
	void edit(std::ostream& flux, std::string before, std::string sep, int dim);

	/// read general matrix from input file
	void input(std::ifstream & fi);
	/// read general matrix from input file
	void input(std::ifstream & fi, int dim);
	virtual void input(float ** variances);

	
	/// compute general matrix SVD decomposition
	void computeSVD(DiagMatrix* & S, GeneralMatrix* & O);

	/// compute Shape as diag(Ot . this . O ) / diviseur
	void computeShape_as__diag_Ot_this_O(DiagMatrix* & Shape, GeneralMatrix* & Ori, float diviseur = 1.0);

	/// compute this as : multi * (O * S * O' )
	void compute_as__multi_O_S_O(float multi, GeneralMatrix* & O, DiagMatrix *& S);

	/// compute this as O * S *O'
	void compute_as_O_S_O(GeneralMatrix* & O, float* & S_store);

	/// compute trace of this
	float computeTrace();

	/// compute this as matrix * matrix'
	void compute_as_M_tM(GeneralMatrix* M, int d);

	/// compute this as matrix * vector
	void compute_as_M_V(GeneralMatrix* M, float * V);
	/// compute this as vector multiplied by matrix
	void multiply(float * V, int nk, GeneralMatrix * Q);

	/// compute M as : M = ( O * S^{-1} * O' ) * this
	void compute_M_as__O_Sinverse_Ot_this(GeneralMatrix & M, GeneralMatrix* & O, DiagMatrix* & S);
	float compute_trace_W_C(Matrix * C);
	//  void computeShape_as__diag_Ot_this_O(XEMDiagMatrix* & Shape, XEMGeneralMatrix* & Ori, float diviseur = 1.0);
	/// gives : det(diag(this))
	float detDiag(Exception& errorType);

	/// trace( this * O * S^{-1} * O' )
	float trace_this_O_Sm1_O(GeneralMatrix* & O, DiagMatrix* & S);

	//void refreshStore();  

	void setSymmetricStore(float * store);
	void setGeneralStore(float * store);
	void setDiagonalStore(float * store);
	void setSphericalStore(float store);

	float** storeToArray() const;

protected:

	// General matrix as in mathematical library
	MATH::Matrix * _value;

	float * _store;

	int _s_storeDim;
};

// TODO static :
// int XEMGeneralMatrix::_s_storeDim = 0;

inline float * GeneralMatrix::getStore() {
	return _store;
}

inline MATH::Matrix * GeneralMatrix::getValue() {
	return _value;
}

inline int GeneralMatrix::getStoreDim() {
	return _s_storeDim;
}

inline void GeneralMatrix::setSymmetricStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

inline void GeneralMatrix::setSphericalStore(float store) {
	THROW(OtherException, wrongMatrixType);
}

inline void GeneralMatrix::setGeneralStore(float * store) {
	//_store = store;
	recopyTab(store, _store, _s_storeDim);
}

inline void GeneralMatrix::setDiagonalStore(float * store) {
	THROW(OtherException, wrongMatrixType);
}

/* TODO static
inline void GeneralMatrix::initiate(){
  _s_storeDim = _s_pbDimension * _s_pbDimension / 2; 
} 
 */

}

#endif
