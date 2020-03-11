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

    errorOut computeGamma( const variableVector &deformationGradient, const variableVector &gradXi,
                           variableVector &Gamma );

    errorOut computeGamma( const variableVector &deformationGradient, const variableVector &gradXi,
                           variableVector &Gamma, variableMatrix &dGammadF, variableMatrix &dGammadGradXi );

    errorOut computeMicroStrain( const variableVector &Psi, variableVector &microStrain );

    errorOut computeMicroStrain( const variableVector &Psi, variableVector &microStrain,
                                 variableMatrix &dMicroStraindPsi );


    errorOut pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress );

    errorOut pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress, 
                                   variableMatrix &dcauchyStressdPK2Stress,
                                   variableMatrix &dcauchyStressdDeformationGradient );

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
                                           variableType &detF,
                                           variableVector &higherOrderStress );

    errorOut pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress,
                                           variableMatrix &dHigherOrderStressdReferenceHigherOrderStress,
                                           variableMatrix &dHigherOrderStressdDeformationGradient,
                                           variableMatrix &dHigherOrderStressdMicroDeformation );

    errorOut computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress );

    errorOut computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress,
                                                 variableMatrix &dDeviatoricStressdStress );

    errorOut computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen,  variableType &pressure );

    errorOut computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen, variableType &pressure,
                                                        variableVector &dpdStress, variableVector &dpdRCG );

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress );

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &invRCG, variableType &CS,
                                                          variableVector &deviatoricSecondOrderReferenceStress );

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG );

    errorOut computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress );
    
    errorOut computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress,
                                                 variableMatrix &dDeviatoricHigherOrderStressdHigherOrderStress);

    errorOut computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress, 
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure );

    errorOut computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress, 
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableMatrix &dpdM, variableMatrix &dpdC );

    errorOut computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableMatrix &dpdM, variableMatrix &dpdC,
                                                        variableMatrix &d2pdMdC );

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress );
    
    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG );
}

#endif
