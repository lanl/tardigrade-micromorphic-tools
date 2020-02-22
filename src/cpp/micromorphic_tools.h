/*!
 * Micromorphic Tools.h
 * ====================
 *
 * A collection of tools and utilities which seek to make the 
 * implementation of micromorphic continuum mechanics based constitutive 
 * theories easier.
 *
 */

#ifndef MICROMORPHIC_TOOLS_H
#define MICROMORHPIC_TOOLS_H

#include<error_tools.h>
#define USE_EIGEN
#include<vector_tools.h>
#include<constitutive_tools.h>

namespace micromorphicTools{

    typedef double variableType;
    typedef std::vector< variableType > variableVector;
    typedef std::vector< variableVector > variableMatrix;

    typedef double parameterType;
    typedef std::vector< parameterType > parameterVector;
    typedef std::vector< parameterVector > parameterMatrix;

    typedef double constantType;
    typedef std::vector< constantType > constantVector;
    typedef std::vector< constantVector > constantMatrix;

    typedef errorTools::Node errorNode;
    typedef errorNode* errorOut;

    errorOut computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
		         variableVector &Psi );

    errorOut computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
		         variableVector &Psi, variableMatrix &dPsidF, variableMatrix &dPsidXi );

    errorOut computeMicroStrain( const variableVector &Psi, variableVector &microStrain );

    errorOut computeMicroStrain( const variableVector &Psi, variableVector &microStrain,
		                 variableMatrix &dMicroStraindPsi );

    errorOut pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
		                              const variableVector &deformationGradient,
					      variableType &detF, variableVector &microStress );

    errorOut pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
		                              const variableVector &deformationGradient,
					      variableVector &microStress );

    errorOut pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
		                              const variableVector &deformationGradient,
					      variableVector &microStress, 
                                              variableMatrix &dMicroStressdReferenceMicroStress,
                                              variableMatrix &dMicroStressdDeformationGradient );

    errorOut pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress );

    errorOut pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress,
                                           variableMatrix &dHigherOrderStressdReferenceHigherOrderStress,
                                           variableMatrix &dHigherOrderStressdDeformationGradient,
                                           variableMatrix &dHigherOrderStressdMicroDeformation );
    
    errorOut computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress );
    
    errorOut computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress,
                                                 variableMatrix &dDeviatoricHigherOrderStressdHigherOrderStress);

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &higherOrderStress,
                                                          const variableVector &rightGreenLagrangeDeformation,
                                                          variableVector &deviatoricHigherOrderStress );
    
    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &higherOrderStress,
                                                          const variableVector &rightGreenLagrangeDeformation,
                                                          variableVector &deviatoricHigherOrderStress,
                                                          variableMatrix &dDeviatoricHigherOrderStressdHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG );
}

#endif
