/***************************************************************************
                             EXAMPLES/clustering/clustering_4.cpp  description
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

#include <mixmod.h>
#include "mixmod/Utilities/ExampleDataUtil.h"

//--------------------------
//--------------------------
// Clustering example
// Quantitative Data
//--------------------------
//--------------------------
void clustering_4_example() {

    cout<<"-----------------------------------------------------------------------"<<endl;
    cout<<"Clustering example : "<<endl;
    cout<<"  - iris data (quantitative data)"<<endl;
    cout<<"  - strategy initialisation : SMALL_EM"<<endl;
    cout<<"     - nbTryInInit = 20"<<endl;
    cout<<"     - nbIterationInInit = 10"<<endl;
    cout<<"-----------------------------------------------------------------------"<<endl<<endl;

    // will be executed from BIN folder
    string filename("../EXAMPLES/data/iris.train");

    
    // 1. create XEM::ClusteringInput
    //--------------------------

    // nbCluster contains the numbers of clusters to be tested.
    // Here we assume that there are 2, 3 or 4 clusters.
    vector<int64_t> nbCluster;
    nbCluster.push_back (2);
    nbCluster.push_back (3);
    nbCluster.push_back (4);

    // XEM::ClusteringInput
    XEM::ClusteringInput* cInput = XEM::getClusteringInput (filename, nbCluster);

    // Initialisation example  :
    // - SMALL_EM (defaut initialisation)
    //   - nbTryInInit = 20 // number of SMALL_EM runs
    //   - nbIterationInInit = 10; // maximum number of iterations for every SMALL_EM
    cInput->getStrategy()->setNbTryInInit ( 20 );
    cInput->getStrategy()->setNbIterationInInit ( 10 );

    // Finalize
    cInput->finalize();

    // edit cInput
    cout<<"-----------------------------------------------------------------------"<<endl;
    cout<<" Input : "<<endl;
    cout<<"-----------------------------------------------------------------------"<<endl<<endl;
    cInput->edit ( cout );
    cout<<endl;


    // 2. XEM::ClusteringMain
    //---------------------

    // 3.1. Create XEM::ClusteringMain
    XEM::ClusteringMain cMain ( cInput );

    // 2.2. run XEM::ClusteringMain
    cMain.run();

    //2.3. Create a new XEM::ClusteringOutput object
    // all Outputs are in XEM::ClusteringOutput
    XEM::ClusteringOutput * cOutput = cMain.getOutput();


    // 3. post treatment example
    //------------------------------
    // sort (using BIC - the default criterion)
    cOutput->sort ( XEM::BIC );


    if ( cOutput->atLeastOneEstimationNoError() ) {
        // get the best XEM::ClusteringModelOutput
        XEM::ClusteringModelOutput * cMOutput = cOutput->getClusteringModelOutput().front();
        XEM::ParameterDescription * paramDescription = cMOutput->getParameterDescription();

        cout<<"-----------------------------------------------------------------------"<<endl;
        cout<<"Best model is "<<endl;
        cout<<" - nbCluster : "<<paramDescription->getNbCluster() <<endl;
        cout<<"-----------------------------------------------------------------------"<<endl<<endl;

        cout<<"-----------------------------------------------------------------------"<<endl;
        cout<<"Parameters display"<<endl;

        XEM::Parameter * param = paramDescription->getParameter();
        // print out parameters
        param->edit();
        // print out criterion values
        for ( int64_t iCriterion=0; iCriterion < cInput->getCriterionName().size(); iCriterion++ )
            cMOutput->getCriterionOutput ( cInput->getCriterionName ( iCriterion ) ).editTypeAndValue ( std::cout );
    }
    cout<<"-----------------------------------------------------------------------"<<endl;

    
    // release memory
    delete cInput;
}
