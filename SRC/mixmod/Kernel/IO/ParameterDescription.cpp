/***************************************************************************
                             SRC/mixmod/Kernel/IO/ParameterDescription.cpp  description
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

#include "mixmod/Kernel/IO/ParameterDescription.h"
#include "mixmod/Kernel/Parameter/GaussianGeneralParameter.h"
#include "mixmod/Kernel/Parameter/BinaryEkjhParameter.h"
#include "mixmod/Kernel/Parameter/CompositeParameter.h"
#include "mixmod/Kernel/Model/Model.h"
#include "mixmod/Kernel/Model/ModelType.h"
#include "mixmod/Kernel/IO/ModelOutput.h"
#include "mixmod/Kernel/Parameter/Parameter.h"
#include "mixmod/Kernel/IO/Data.h"

namespace XEM {

//------------
// Constructor by default
//------------
ParameterDescription::ParameterDescription() {
	_parameter = NULL;
}

//-------------------------------------
// Constructor after an estimation->run
//--------------------------------------
ParameterDescription::ParameterDescription(Model* iEstimation) {

	if (iEstimation) {
		_infoName = "Parameter";
		//_nbSample = iEstimation->getNbSample();
		_nbCluster = iEstimation->getNbCluster();
		_nbVariable = iEstimation->getData()->_pbDimension;
		_format = FormatNumeric::defaultFormatNumericFile;
		_filename = "";
		_modelType = new ModelType(*iEstimation->getModelType());
		_parameter = iEstimation->getParameter()->clone();
		if (isBinary(_modelType->_nameModel)) {
			BinaryParameter * bParameter = 
					dynamic_cast<BinaryParameter*> (iEstimation->getParameter());
			recopyTabToVector(bParameter->getTabNbModality(), _nbFactor, _nbCluster);
		}
		saveNumericValues(_filename);
	}
	else {
		THROW(OtherException, nullPointerError);
	}
}

//-------------------------------------
// Constructor after an estimation->run
//--------------------------------------
ParameterDescription::ParameterDescription(ModelOutput* iEstimation) {

	if (iEstimation) {
		_infoName = "Parameter";
		//_nbSample = iEstimation->getNbSample();
		_nbCluster = iEstimation->getNbCluster();
		_nbVariable = iEstimation->getParameterDescription()->getNbVariable();
		_format = FormatNumeric::defaultFormatNumericFile;
		_filename = "";
		_modelType = new ModelType(*iEstimation->getParameterDescription()->getModelType());
		_parameter = iEstimation->getParameterDescription()->getParameter()->clone();
		if (isBinary(_modelType->_nameModel)) {
			BinaryParameter * bParameter = dynamic_cast<BinaryParameter*> 
					(iEstimation->getParameterDescription()->getParameter());
			recopyTabToVector(bParameter->getTabNbModality(), _nbFactor, _nbCluster);
		}
	}
	else {
		THROW(OtherException, nullPointerError);
	}
}

// ---------------------------
//constructor by initialization for Binary
// ---------------------------
ParameterDescription::ParameterDescription(
		int nbCluster, 
		int nbVariable, 
		std::vector< int > nbFactor, 
		FormatNumeric::FormatNumericFile format, 
		std::string filename, 
		std::string infoName, 
		ModelName& modelName) 
{
	_infoName = "Parameter";
	_nbVariable = nbVariable;
	_filename = filename;
	_nbCluster = nbCluster;
	_format = format;
	_nbFactor = nbFactor;
	_modelType = new ModelType(modelName);
	std::ifstream fi(filename.c_str(), ios::in);
	if (!fi.is_open()) {
		THROW(InputException, wrongLabelFileName);
	}
	int * tabNbFactor = new int[_nbVariable];
	recopyVectorToTab(nbFactor, tabNbFactor);
	// create _parameter : always a XEMBinaryEkjhParameter is created
	_parameter = new BinaryEkjhParameter(
			nbCluster, _nbVariable, _modelType, tabNbFactor, filename);
}

// -----------------------------------------
//constructor by initialization for Gaussian
// ----------------------------------------
ParameterDescription::ParameterDescription(
		int nbCluster, 
		int nbVariable, 
		FormatNumeric::FormatNumericFile format, 
		std::string filename, 
		std::string infoName, 
		ModelName& modelName) 
{
	_infoName = "Parameter";
	_nbVariable = nbVariable;
	_filename = filename;
	_nbCluster = nbCluster;
	_format = format;
	//_nbFactor  empty
	_modelType = new ModelType(modelName);
	std::ifstream fi(filename.c_str(), ios::in);
	if (!fi.is_open()) {
		THROW(InputException, wrongLabelFileName);
	}
	// create _parameter : always a XEMGaussianGeneralParameter is created
	_parameter = new GaussianGeneralParameter(nbCluster, _nbVariable, _modelType, filename);
}

// ---------------------------
//constructor by initialization for Composite
// ---------------------------
ParameterDescription::ParameterDescription(
		int nbCluster,
		int nbVariable_binary,
		int nbVariable_gaussian,
		std::vector< int > nbFactor,
		FormatNumeric::FormatNumericFile format,
		std::string filename,
		std::string infoName,
		ModelName& modelName)
{
	_infoName = "Parameter";
	_nbVariable = nbVariable_binary + nbVariable_gaussian;
	_filename = filename;
	_nbCluster = nbCluster;
	_format = format;
	_nbFactor = nbFactor;
	_modelType = new ModelType(modelName);
	std::ifstream fi(filename.c_str(), ios::in);
	if (!fi.is_open()) {
		THROW(InputException, wrongLabelFileName);
	}
	int * tabNbFactor = new int[nbVariable_binary];
	recopyVectorToTab(nbFactor, tabNbFactor);

	GaussianGeneralParameter* gaussian_parameter =
		new GaussianGeneralParameter(nbCluster, nbVariable_gaussian, _modelType, filename, nbVariable_binary, nbFactor);
	BinaryEkjhParameter* binary_parameter =
		new BinaryEkjhParameter(
			nbCluster, nbVariable_binary, _modelType, tabNbFactor, filename);

	_parameter = new CompositeParameter(gaussian_parameter, binary_parameter, _modelType);
}

//constructor using XEMParameter
ParameterDescription::ParameterDescription(Parameter * iparam) {
	_parameter = iparam->clone();
	_infoName = "Parameter";
	_nbCluster = iparam->getNbCluster();
	_nbVariable = iparam->getPbDimension();
	_format = FormatNumeric::defaultFormatNumericFile;
	_filename = "";
	_modelType = new ModelType(*iparam->getModelType());
}

//constructor for binary data
ParameterDescription::ParameterDescription(
	int nbCluster,
	int nbVariable,
	ModelName& modelName,
	float * proportions,
	float ** centers,
	float *** scatters,
	std::vector< int> nbFactor)
{
	_infoName = "Parameter";
	_nbVariable = nbVariable;
	_filename = "";
	_nbCluster = nbCluster;
	_format = FormatNumeric::defaultFormatNumericFile;
	int* tabNbFactor = new int[nbVariable];
	recopyVectorToTab(nbFactor, tabNbFactor);
	_modelType = new ModelType(modelName);

	// create _parameter : always a XEMBinaryEkjhParameter is created
	_parameter = new BinaryEkjhParameter(nbCluster, _nbVariable , _modelType,
		tabNbFactor, proportions, centers, scatters);

	//TODO ?
	//delete[] tabNbFactor;
}

//constructor for Gaussian data
ParameterDescription::ParameterDescription(
	int nbCluster,
	int nbVariable,
	ModelName& modelName,
	float * proportions,
	float ** means,
	float *** variances)
{
	_infoName = "Parameter";
	_nbVariable = nbVariable;
	_filename = "";
	_nbCluster = nbCluster;
	_format = FormatNumeric::defaultFormatNumericFile;
	//_nbFactor is empty
	_modelType = new ModelType(modelName);

	// create _parameter : always a XEMGaussianGeneralParameter is created
	_parameter = new GaussianGeneralParameter(nbCluster, _nbVariable,
		_modelType, proportions, means, variances);
}

//constructor for Heterogeneous data
ParameterDescription::ParameterDescription(
	int nbCluster,
	int nbBinaryVariable,
	int nbGaussianVariable,
	ModelName& modelName,
	float * proportions,
	float ** centers,
	float *** scatters,
	float ** means,
	float *** variances,
	std::vector< int> nbFactor)
{
	_infoName = "Parameter";
	_nbVariable = nbBinaryVariable + nbGaussianVariable;
	_filename = "";
	_nbCluster = nbCluster;
	_format = FormatNumeric::defaultFormatNumericFile;
	int * tabNbFactor = new int[nbBinaryVariable];
	recopyVectorToTab(nbFactor, tabNbFactor);
	ModelType* _binarymodelType = new ModelType(getBinaryModelNamefromHeterogeneous(modelName));
	ModelType* _gaussianmodelType = new ModelType(getBinaryModelNamefromHeterogeneous(modelName));
	_modelType = new ModelType(modelName);
	// create _parameter : always a XEMGaussianGeneralParameter is created
	Parameter* g_parameter = new GaussianGeneralParameter(
			nbCluster, nbGaussianVariable, _gaussianmodelType, proportions, means, variances);
	// create _parameter : always a XEMBinaryEkjhParameter is created
	Parameter* b_paramemter = new BinaryEkjhParameter(
			nbCluster, nbBinaryVariable, _binarymodelType, 
			tabNbFactor, proportions, centers, scatters);
	_parameter = new CompositeParameter(g_parameter, b_paramemter, _modelType);

	//release memory
	delete g_parameter;
	delete b_paramemter;
	delete _binarymodelType;
	delete _gaussianmodelType;
}

//------------
// Destructor 
//------------
ParameterDescription::~ParameterDescription() {
	if (_modelType) delete _modelType;
	if (_parameter) delete _parameter;
}

/// Comparison operator
bool ParameterDescription::operator==(ParameterDescription & paramDescription) const {
	if (_infoName != paramDescription.getInfoName()) return false;
	if (_nbVariable != paramDescription.getNbVariable()) return false;
	if (_filename != paramDescription.getFilename()) return false;
	if (_nbCluster != paramDescription.getNbCluster()) return false;
	if (_format != paramDescription.getFormat()) return false;
	if (!(_modelType == paramDescription.getModelType())) return false;
	for (unsigned int i = 0; i < _nbFactor.size(); ++i) {
		if (_nbFactor[i] != paramDescription.getTabNbFactor()[i]) return false;
	}
	if (!(_parameter == paramDescription.getParameter())) return false;
	return true;
}

//--------
// ostream
//--------
void ParameterDescription::saveNumericValues(std::string fileName) {
	//if (_filename==""){
	std::ofstream fo(fileName.c_str(), ios::out);
	_parameter->edit(fo);
	_filename = fileName;
	//}
	/* else : if _fileName!="", paprameterDescription has been created by a XML file.
	In this case, the numeric file already exists. 
	 */
}

}
