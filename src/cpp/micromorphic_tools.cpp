/*!
 * Micromorphic Tools.cpp
 * ====================
 *
 * A collection of tools and utilities which seek to make the 
 * implementation of micromorphic continuum mechanics based constitutive 
 * theories easier.
 *
 */

#include<micromorphic_tools.h>

namespace micromorphicTools{

    errorOut computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
                         variableVector &Psi ){
        /*!
         * Computes the micromorphic quantity psi defined as:
         * \Psi_{IJ} = F_{i I} \Xi_{i J}
         *
         * :param const variableVector &deformationGradient: The deformation gradient.
         * :param const variableVector &microDeformation: The micro-deformation.
         * :param variableVector &Psi: The micro-displacement metric Psi.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( deformationGradient.size() != dim * dim ){
            return new errorNode( "computePsi", "The deformation gradient doesn't have the correct size" );
        }

        if ( microDeformation.size() != dim * dim ){
            return new errorNode( "computePsi", "The micro-deformation doesn't have the correct size" );
        }

        Psi = vectorTools::matrixMultiply( deformationGradient, microDeformation,
                                           dim, dim, dim, dim, 1, 0 );

        return NULL;
    }

    errorOut computePsi( const variableVector &deformationGradient, const variableVector &microDeformation,
                         variableVector &Psi, variableMatrix &dPsidF, variableMatrix &dPsidXi ){
        /*!
         * Computes the micromorphic quantity psi defined as:
         * \Psi_{IJ} = F_{i I} \Xi_{i J}
         *
         * along with the jacobians
         *
         * \frac{ \partial \Psi_{IJ} }{ \partial F_{k K} } = \delta_{I K} \Xi_{k J}
         * \frac{ \partial \Psi_{IJ} }{ \partial \Xi_{k K} } = F_{k I} \delta_{J K}
         *
         * :param const variableVector &deformationGradient: The deformation gradient.
         * :param const variableVector &microDeformation: The micro-deformation.
         * :param variableVector &Psi: The micro-displacement metric Psi
         * :param variableMatrix &dPsidF: The jacobian of Psi w.r.t. the deformation gradient.
         * :param variableMatrix &dPsidXi: The jacobian of Psi w.r.t. the micro-deformation.
         */

        //Assume 3d
        unsigned int dim = 3;

        errorOut error = computePsi( deformationGradient, microDeformation, Psi );

        if (error){
            errorOut result = new errorNode( "computePsi (jacobian)", "Error in the computation of Psi" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim, 0 );
        vectorTools::eye( eye );

        dPsidF  = variableMatrix( Psi.size(), variableVector( deformationGradient.size(), 0 ) );
        dPsidXi = variableMatrix( Psi.size(), variableVector( microDeformation.size(), 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        dPsidF[ dim * I + J ][ dim * k + K ] = eye[ dim * I + K ] * microDeformation[ dim * k + J ];
                        dPsidXi[ dim * I + J ][ dim * k + K ] = deformationGradient[ dim * k + I ] * eye[ dim * J + K ];
                    }
                }
            }
        }
        return NULL;
    }

    errorOut computeGamma( const variableVector &deformationGradient, const variableVector &gradXi,
                           variableVector &Gamma ){
        /*!
         * Compute the deformation metric Gamma:
         *
         * Gamma_{IJK} = F_{iI} \Xi_{iJ,K}
         *
         * :param const variableVector &deformationGradient: The deformation gradient.
         * :param const variableVector &gradXi: The gradient of the micro-deformation tensor
         *     w.r.t. the reference configuration.
         * :param variableVector &Gamma: The micromorphic deformation metric Gamma.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( deformationGradient.size() != dim * dim ){
            return new errorNode("computeGamma", "The deformation gradient isn't the right size");
        }

        if ( gradXi.size() != dim * dim * dim ){
            return new errorNode("computeGamma", "The micro-deformation gradient isn't the right size");
        }

        Gamma = variableVector( dim * dim * dim, 0 );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int i = 0; i < dim; i++ ){
                        Gamma[ dim * dim * I + dim * J + K ] += deformationGradient[ dim * i + I ] * gradXi[ dim * dim * i + dim * J + K ];
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeGamma( const variableVector &deformationGradient, const variableVector &gradXi,
                           variableVector &Gamma, variableMatrix &dGammadF, variableMatrix &dGammadGradXi ){
        /*!
         * Compute the deformation metric Gamma:
         *
         * Gamma_{IJK} = F_{iI} \Xi_{iJ,K}
         *
         * Also return the Jacobians
         * \frac{ \partial Gamma_{IJK} }{ \partial F_{lL} } = \delta_{IL} \Xi_{lJ,K}
         * \frac{ \partial Gamma_{IJK} }{ \partial \Xi_{lL,M} } = F_{lI} \delta_{JL} \delta_{KM}
         *
         * :param const variableVector &deformationGradient: The deformation gradient.
         * :param const variableVector &gradXi: The gradient of the micro-deformation tensor
         *     w.r.t. the reference configuration.
         * :param variableVector &Gamma: The micromorphic deformation metric Gamma.
         * :param variableMatrix &dGammadF: The gradient of Gamma w.r.t. the deformation gradient.
         * :param variableMatrix &dGammadGradXi: The gradient of Gamma w.r.t. the gradient of Xi in the reference 
         *     configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        errorOut error = computeGamma( deformationGradient, gradXi, Gamma );

        if ( error ){
            errorOut result = new errorNode("computeGamma (jacobian)", "Error in computation of Gamma");
            result->addNext(error);
            return result;
        }

        dGammadF      = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dGammadGradXi = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        constantVector eye( dim * dim, 0 );
        vectorTools::eye( eye );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        for ( unsigned int L = 0; L < dim; L++ ){
                            dGammadF[ dim * dim * I + dim * J + K ][ dim * l + L ] = eye[ dim * I + L ]
                                                                                   * gradXi[ dim * dim * l + dim * J + K ];
                            for ( unsigned int M = 0; M < dim; M++ ){
                                dGammadGradXi[ dim * dim * I + dim * J + K ][ dim * dim * l + dim * L + M ] = deformationGradient[ dim * l + I ]
                                                                                                            * eye[ dim * J + L ] 
                                                                                                            * eye[ dim * K + M ];
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeMicroStrain( const variableVector &Psi, variableVector &microStrain ){
        /*!
         * Compute the microstrain defined as:
         * \Epsilon_{IJ} = \Psi_{IJ} - eye_{IJ}
         *
         * :param const variableVector &Psi: The micro-deformation metric Psi.
         * :param variableVector &microStrain: The micro-strain.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( Psi.size() != dim * dim ){
            return new errorNode( "computeMicroStrain", "Psi is not of the correct size" );
        }

        constantVector eye( dim * dim, 0 );
        vectorTools::eye( eye );

        microStrain = Psi - eye;

        return NULL;
    }

    errorOut computeMicroStrain( const variableVector &Psi, variableVector &microStrain,
                                 variableMatrix &dMicroStraindPsi ){
        /*!
         * Compute the microstrain defined as:
         * \Epsilon_{IJ} = \Psi_{IJ} - eye_{IJ}
         *
         * and also compute the jacobian
         * \frac{ \partial \Epsilon_{IJ} }{ \partial \Psi_{KL} } = \delta_{IK} \delta_{JL}
         *
         * :param const variableVector &Psi: The micro-deformation metric Psi.
         * :param variableVector &microStrain: The micro-strain.
         * :param variableMatrix &dMicroStraindPsi: The jacobian of the micro-strain.
         */

        errorOut error = computeMicroStrain( Psi, microStrain );

        if (error){
            errorOut result = new errorNode( "computeMicroStrain (jacobian)", "Error in computation of micro-strain" );
            result->addNext( error );
            return result;
        }

        dMicroStraindPsi = vectorTools::eye< variableType >( Psi.size() );

        return NULL;
    }

    errorOut pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress ){
        /*!
         * Push forward the PK2-stress in the reference configuration to the 
         * configuration (the Cauchy stress) indicated by the deformation gradient.
         *
         * \cauchy_{ij} = (1 / J ) F_{i I} S_{I J} F_{j J}
         *
         * :param const variableVector &PK2Stress: The PK2 stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * :param variableVector &cauchyStress: The Cauchy stress in the current 
         *     configuration.
         */

        variableType detF;
        errorOut error = pushForwardReferenceMicroStress( PK2Stress, deformationGradient, 
                                                          detF, cauchyStress );

        if ( error ){
            errorOut result = new errorNode( "pushForwardPK2Stress",
                                             "Error in push-forward operation (micro-stress and PK2 are identical)" );
            result->addNext( error );
            return result;
        }
        return NULL;
    }

    errorOut pushForwardPK2Stress( const variableVector &PK2Stress,
                                   const variableVector &deformationGradient,
                                   variableVector &cauchyStress,
                                   variableMatrix &dCauchyStressdPK2Stress,
                                   variableMatrix &dCauchyStressdDeformationGradient ){
        /*!
         * Push forward the PK2 stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * \sigma_{ij} = (1 / J ) F_{iI} S_{IJ} F_{jJ}
         *
         * Also computes the jacobians:
         * \frac{ \partial \cauchy_{ij} }{\partial \Sigma_{KL} } = ( 1 / J ) F_{iK} F_{jL}
         * \frac{ \partial \cauchy_{ij} }{\partial F_{kK} } = ( \delta_{i k} \delta_{I K} S_{I J} F_{j J}
         *                                                  + F_{i I} S_{I J} \delta_{j k} \delta_{J K}
         *                                                  - \cauchy_{i j} dDetFdF_{kK} ) / J
         *
         * :param const variableVector &referenceMicroStress: The PK2 stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * :param variableVector &cauchyStress: The Cauchy stress in the current 
         *     configuration.
         * :param variableMatrix &CauchyStressdReferenceMicroStress: The jacobian of 
         *     the Cauchy w.r.t. the PK2 tress in the reference configuration.
         * :param variableMatrix &dCauchyStressdDeformationGradient: The jacobian of 
         *     the Cauchy stress w.r.t. the deformation gradient.
         */

        errorOut error = pushForwardReferenceMicroStress( PK2Stress, deformationGradient, cauchyStress,
                                                          dCauchyStressdPK2Stress,
                                                          dCauchyStressdDeformationGradient );

        if ( error ){
            errorOut result = new errorNode( "pushForwardPK2Stress (jacobian)",
                                             "Error in push-forward operation (micro-stress and PK2 are identical)" );
            result->addNext( error );
            return result;
        }
        return NULL;
    }


    errorOut pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableVector &microStress ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * s_{ij} = (1 / J ) F_{i I} \Sigma_{I J} F_{j J}
         *
         * :param const variableVector &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * :param variableVector &microStress: The micro-stress in the current 
         *     configuration.
         */

        variableType detF;
        return pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, 
                                                detF, microStress );

    }

    errorOut pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableType &detF, variableVector &microStress ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * s_{ij} = (1 / J ) F_{i I} \Sigma_{I J} F_{j J}
         *
         * :param const variableVector &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * :param const variableType &detF: The determinant of the deformation gradient.
         * :param variableVector &microStress: The micro-stress in the current 
         *     configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( referenceMicroStress.size() != dim * dim ){
            return new errorNode("pushForwardReferenceMicroStress", "The reference micro-stress has an incorrect size");
        }

        if ( deformationGradient.size() != dim * dim ){
            return new errorNode("pushForwardReferenceMicroStress", "The deformation gradient has an incorrect size");
        }

        microStress = variableVector( dim * dim, 0 );

        detF = vectorTools::determinant( deformationGradient, dim, dim );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int I = 0; I < dim; I++ ){
                    for ( unsigned int J = 0; J < dim; J++ ){
                        microStress[ dim * i + j ] += deformationGradient[ dim * i + I ]
                                                    * referenceMicroStress[ dim * I + J ]
                                                    * deformationGradient[ dim * j + J ] / detF;
                    }
                }
            }
        }

        return NULL;
    }

    errorOut pushForwardReferenceMicroStress( const variableVector &referenceMicroStress,
                                              const variableVector &deformationGradient,
                                              variableVector &microStress,
                                              variableMatrix &dMicroStressdReferenceMicroStress,
                                              variableMatrix &dMicroStressdDeformationGradient ){
        /*!
         * Push forward the micro-stress in the reference configuration to the 
         * configuration indicated by the deformation gradient.
         *
         * s_{ij} = (1 / J ) F_{iI} \Sigma_{IJ} F_{jJ}
         *
         * Also computes the jacobians:
         * \frac{ \partial s_{ij} }{\partial \Sigma_{KL} } = ( 1 / J ) F_{iK} F_{jL}
         * \frac{ \partial s_{ij} }{\partial F_{kK} } = ( \delta_{i k} \delta_{I K} \Sigma_{I J} F_{j J}
         *                                              + F_{i I} \Sigma_{I J} \delta_{j k} \delta_{J K}
         *                                              - s_{i j} dDetFdF_{kK} ) / J
         *
         * :param const variableVector &referenceMicroStress: The micro-stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient 
         *     mapping between configurations.
         * :param variableVector &microStress: The micro-stress in the current 
         *     configuration.
         * :param variableMatrix &dmicroStressdReferenceMicroStress: The jacobian of 
         *     the micro-stress w.r.t. the micro-stress in the reference configuration.
         * :param variableMatrix &dmicroSTressdDeformationGradient: The jacobian of 
         *     the micro-stress w.r.t. the deformation gradient.
         */

        //Assume 3d
        unsigned int dim = 3;

        variableType detF;
        errorOut error = pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient,
                                                          detF, microStress );

        if (error){
            errorOut result = new errorNode( "pushForwardReferenceMicroStress (jacobian)", "Error in computation of push forward of micro-stress" );
            result->addNext(error);
            return result;
        }

        //Assemble the jacobian of the determinant of the deformation gradient
        variableVector inverseDeformationGradient = vectorTools::inverse( deformationGradient, dim, dim );

        variableVector dDetFdF( dim * dim, 0 );

        for (unsigned int i = 0; i < dim; i++ ){
            for (unsigned int I = 0; I < dim; I++ ){
                dDetFdF[ dim * i + I ] = inverseDeformationGradient[ dim * I + i ] * detF;
            }
        }

        //Assemble the jacobians
        dMicroStressdReferenceMicroStress = variableMatrix( microStress.size(), variableVector( referenceMicroStress.size(), 0 ) );
        dMicroStressdDeformationGradient = variableMatrix( microStress.size(), variableVector( deformationGradient.size(), 0 ) );

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int K = 0; K < dim; K++ ){
                        dMicroStressdReferenceMicroStress[ dim * i + j ][ dim * k + K ] = deformationGradient[ dim * i + k ]
                                                                                        * deformationGradient[ dim * j + K ] / detF;
                        
                        for ( unsigned int I = 0; I < dim; I++ ){

                            dMicroStressdDeformationGradient[ dim * i + j ][ dim * k + K] += eye[ dim * i + k ] * referenceMicroStress[ dim * K + I ] * deformationGradient[ dim * j + I ]
                                                                                          + deformationGradient[ dim * i + I ] * referenceMicroStress[ dim * I + K ] * eye[ dim * j + k ];
                        }

                        dMicroStressdDeformationGradient[ dim * i + j ][ dim * k + K] -= microStress[ dim * i + j ] * dDetFdF[ dim * k + K ];
                        dMicroStressdDeformationGradient[ dim * i + j ][ dim * k + K] /= detF; 
                    }
                }
            }
        }

        return NULL;
    }

    errorOut pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress ){
        /*!
         * Compute the push-forward operation on the higher order stress.
         *
         * m_{ijk} = \frac{1}{J} F_{iI} F_{jJ} \Xi_{kK} M_{IJK}
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * :param const variableVector &microDeformation: The micro-deformation tensor.
         * :param variableVector &higherOrderStress: The higher order stress in the current configuration.
         */

        variableType detF;
        return pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient,
                                             microDeformation, detF, higherOrderStress );
    }

    errorOut pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableType &detF,
                                           variableVector &higherOrderStress ){
        /*!
         * Compute the push-forward operation on the higher order stress.
         *
         * m_{ijk} = \frac{1}{J} F_{iI} F_{jJ} \Xi_{kK} M_{IJK}
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * :param const variableVector &microDeformation: The micro-deformation tensor.
         * :param const variableType &detF: The determinant of the deformation gradient.
         * :param variableVector &higherOrderStress: The higher order stress in the current configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( referenceHigherOrderStress.size() != dim * dim * dim ){
            return new errorNode( "pushForwardHigherOrderStress", "The reference higher order stress doesn't have the correct size" );
        }

        if ( deformationGradient.size() != dim * dim ){
            return new errorNode( "pushForwardHigherOrderStress", "The deformation gradient doesn't have the correct size" );
        }

        if ( microDeformation.size() != dim * dim ){
            return new errorNode( "pushForwardHigherOrderStress", "The micro-deformation doesn't have the correct size" );
        }

        detF = vectorTools::determinant( deformationGradient, dim, dim );

        higherOrderStress = variableVector( dim * dim * dim, 0 );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int I = 0; I < dim; I++ ){
                        for ( unsigned int J = 0; J < dim; J++ ){
                            for ( unsigned int K = 0; K < dim; K++ ){
                                higherOrderStress[ dim * dim * i + dim * j + k ] += deformationGradient[ dim * i + I ]
                                                                                  * deformationGradient[ dim * j + J ]
                                                                                  * microDeformation[ dim * k + K ]
                                                                                  * referenceHigherOrderStress[ dim * dim * I + dim * J + K ];
                            }
                        }
                    }
                    higherOrderStress[ dim * dim * i + dim * j + k ] /= detF;
                }
            }
        }

        return NULL;
    }

    errorOut pushForwardHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                           const variableVector &deformationGradient,
                                           const variableVector &microDeformation,
                                           variableVector &higherOrderStress,
                                           variableMatrix &dHigherOrderStressdReferenceHigherOrderStress,
                                           variableMatrix &dHigherOrderStressdDeformationGradient,
                                           variableMatrix &dHigherOrderStressdMicroDeformation ){
        /*!
         * Compute the push-forward operation on the higher order stress.
         *
         * m_{ijk} = \frac{1}{J} F_{iI} F_{jJ} \Xi_{kK} M_{IJK}
         *
         * Also returns the Jacobians
         *
         * \frac{ \partial m_{ijk} }{ \partial M_{LMN} } = \frac{1}{J} F_{iL} F_{jM} \Xi_{kN}
         * \frac{ \partial m_{ijk} }{ \partial F_{lM} } = \left( \delta_{il} F_{jN} \Xi_{kO} M_{MNO}
         *                                                     + F_{iN} \delta_{jl} \Xi_{kO} M_{NMO}
         *                                                     - m_{ijk} dDetFdF_{lM} \right)/J
         * \frac{ \partial m_{ijk} }{ \partial \Xi_{lM} } = \frac{1}{J} F_{iN} F_{jO} \delta_{kl} M_{NOM}
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * :param const variableVector &deformationGradient: The deformation gradient which maps 
         *     between the reference and current configurations.
         * :param const variableVector &microDeformation: The micro-deformation tensor.
         * :param variableVector &higherOrderStress: The higher order stress in the current configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        variableType detF;

        errorOut error = pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient, microDeformation,
                                                       detF, higherOrderStress );
        
        if (error){
            errorOut result = new errorNode( "pushForwardHigherOrderStress (jacobian)", "Error in computation of push forward of the higher order stress" );
            result->addNext(error);
            return result;
        }

        //Assemble the jacobian of the determinant of the deformation gradient
        variableVector inverseDeformationGradient = vectorTools::inverse( deformationGradient, dim, dim );

        variableVector dDetFdF( dim * dim, 0 );

        for (unsigned int i = 0; i < dim; i++ ){
            for (unsigned int I = 0; I < dim; I++ ){
                dDetFdF[ dim * i + I ] = inverseDeformationGradient[ dim * I + i ] * detF;
            }
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dHigherOrderStressdReferenceHigherOrderStress = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
        dHigherOrderStressdDeformationGradient = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );
        dHigherOrderStressdMicroDeformation = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                dHigherOrderStressdReferenceHigherOrderStress[ dim * dim * i + dim * j + k ][ dim * dim * l + dim * M + N ] = deformationGradient[ dim * i + l ] * deformationGradient[ dim * j + M ] * microDeformation[ dim * k + N] / detF;
                                for ( unsigned int O = 0; O < dim; O++ ){
                                    dHigherOrderStressdDeformationGradient[ dim * dim * i + dim * j + k ][ dim * l + M ] += 
                                        eye[ dim * i + l ] * deformationGradient[ dim * j + N ] * microDeformation[ dim * k + O ] * referenceHigherOrderStress[ dim * dim * M + dim * N + O ]
                                      + deformationGradient[ dim * i + N ] * eye[ dim * j + l ] * microDeformation[ dim * k + O ] * referenceHigherOrderStress[ dim * dim * N + dim * M + O ];
                                    dHigherOrderStressdMicroDeformation[ dim * dim * i + dim * j + k ][ dim * l + M ] += deformationGradient[ dim * i + N ] * deformationGradient[ dim * j + O ] * eye[ dim * k + l ] * referenceHigherOrderStress[ dim * dim * N + dim * O + M ];
                                }
                            }
                            dHigherOrderStressdDeformationGradient[ dim * dim * i + dim * j + k ][ dim * l + M ] -= higherOrderStress[ dim * dim * i + dim * j + k ] * dDetFdF[ dim * l + M ];
                            dHigherOrderStressdDeformationGradient[ dim * dim * i + dim * j + k ][ dim * l + M ] /= detF;
                            dHigherOrderStressdMicroDeformation[ dim * dim * i + dim * j + k ][ dim * l + M ] /= detF;
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress ){
        /*!
         * Compute the deviatoric part of the higher order stress.
         *
         * dev ( m_{ijk} ) = m_{ijk} - ( 1 / 3 ) m_{llk} \delta_{ij}
         *
         * :param const variableVector &higherOrderStress: The higher order stress in the current configuration.
         * :param variableVector &deviatoricHigherOrderStress: The deviatoric part of the higher order stress.
         */

        //Assume 3d
        unsigned int dim = 3;

        if ( higherOrderStress.size() != dim * dim * dim ){
            return new errorNode( "computeDeviatoricHigherOrderStress", "The higher order stress has an incorrect size" );
        }

        deviatoricHigherOrderStress = higherOrderStress;

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        deviatoricHigherOrderStress[ dim * dim * i + dim * j + k ] -= higherOrderStress[ dim * dim * l + dim * l + k ] 
                                                                                    * eye[ dim * i + j ] / 3;
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeDeviatoricHigherOrderStress( const variableVector &higherOrderStress,
                                                 variableVector &deviatoricHigherOrderStress,
                                                 variableMatrix &dDeviatoricHigherOrderStressdHigherOrderStress){
        /*!
         * Compute the deviatoric part of the higher order stress.
         *
         * dev ( m_{ijk} ) = m_{ijk} - ( 1 / 3 ) m_{llk} \delta_{ij}
         *
         * Also compute the Jacobian
         * \frac{ \partial dev ( m_{ijk} ) }{ \partial m_{mno} } = \delta_{im} \delta_{jn} \delta_{ko} - ( 1 / 3 ) \delta_{mn} \delta_{ko} \delta_{ij}
         *
         * :param const variableVector &higherOrderStress: The higher order stress in the current configuration.
         * :param variableVector &deviatoricHigherOrderStress: The deviatoric part of the higher order stress.
         * :param variableMatrix &dDeviatoricHigherOrderStressdHigherOrderStress: The gradient of the deviatoric part of the 
         *     higher order stress w.r.t. the higher order stress.
         */

        //Assume 3d
        unsigned int dim = 3;

        errorOut error = computeDeviatoricHigherOrderStress( higherOrderStress, deviatoricHigherOrderStress );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricHigherOrderStress (jacobian)",
                                             "Error in the computation of the deviatoric part of the higher order stress" );
            result->addNext(error);
            return result;
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dDeviatoricHigherOrderStressdHigherOrderStress = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int m = 0; m < dim; m++ ){
                        for ( unsigned int n = 0; n < dim; n++ ){
                            for ( unsigned int o = 0; o < dim; o++ ){
                                dDeviatoricHigherOrderStressdHigherOrderStress[ dim * dim * i + dim * j + k ][ dim * dim * m + dim * n + o ] = eye[ dim * i + m ] * eye[ dim * j + n ] * eye[ dim * k + o ] - eye[ dim * m + n ] * eye[ dim * k + o ] * eye[ dim * i + j ] / 3;
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure ){
        /*!
         * Compute the pressure for a higher-order stress in the reference configuration.
         * $p_K = \frac{1}{3} C_{AB} M_{ABK}$
         *
         * where $C_{AB}$ is the right Cauchy-Green deformation tensor and 
         * M_{ABK} is the higher order stress tensor in the reference configuration.
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation
         *     tensor.
         * :param variableVector &referenceHigherOrderPressure: The higher order pressure.
         */

        //Assume 3D
        unsigned int dim = 3;

        if ( rightCauchyGreenDeformation.size() != dim * dim ){
            return new errorNode( "computeReferenceHigherOrderStressPressure",
                                  "The right Cauchy-Green deformation tensor must have nine terms." );
        }

        if ( referenceHigherOrderStress.size() != dim * dim * dim ){
            return new errorNode( "computeReferenceHigherOrderStressPressure",
                                  "The higher order stress tensor must have 27 terms." );
        }

        referenceHigherOrderPressure = variableVector( dim, 0 );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int A = 0; A < dim; A++ ){
                for ( unsigned int B = 0; B < dim; B++ ){
                    referenceHigherOrderPressure[K] += rightCauchyGreenDeformation[ dim * A + B ]
                                                     * referenceHigherOrderStress[ dim * dim * A + dim * B + K ];
                }
            }
        }

        referenceHigherOrderPressure /= 3;
        return NULL;
    }

    errorOut computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableMatrix &dpdM, variableMatrix &dpdC ){
        /*!
         * Compute the pressure for a higher-order stress in the reference configuration.
         * $p_K = \frac{1}{3} C_{AB} M_{ABK}$
         *
         * Also compute the Jacobians
         * $\frac{ \partial p_K }{ \partial M_{NOP} } = \frac{1}{3} C_{NO} \delta_{KP}$
         * $\frac{ \partial p_K }{ \partial C_{NO} } = \frac{1}{3} M_{NOK}$
         *
         * where $C_{AB}$ is the right Cauchy-Green deformation tensor and 
         * M_{ABK} is the higher order stress tensor in the reference configuration.
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the 
         *     reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation
         *     tensor.
         * :param variableVector &referenceHigherOrderPressure: The higher order pressure.
         * :param variableMatrix &dpdM: The Jacobian of the pressure w.r.t. the higher order stress.
         * :param variableMatrix &dpdC: The Jacobian of the pressure w.r.t. the right Cauchy-Green 
         *     deformation tensor.
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                                    referenceHigherOrderPressure );

        if ( error ){
            errorOut result = new errorNode( "computeReferenceHigherOrderStressPressure (jacobian)",
                                             "Error in computation of reference higher order pressure" );
            result->addNext( error );
            return result;
        }

        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dpdM = variableMatrix( referenceHigherOrderPressure.size(), variableVector( referenceHigherOrderStress.size(), 0 ) );
        dpdC = variableMatrix( referenceHigherOrderPressure.size(), variableVector( rightCauchyGreenDeformation.size(), 0 ) );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int N = 0; N < dim; N++ ){
                for ( unsigned int O = 0; O < dim; O++ ){
                    dpdC[ K ][ dim * N + O ] = referenceHigherOrderStress[ dim * dim * N + dim * O + K ];
                    for ( unsigned int P = 0; P < dim; P++ ){
                        dpdM[ K ][ dim * dim * N + dim * O + P ] = rightCauchyGreenDeformation[ dim * N + O ] * eye[ dim * K + P ];
                    }
                }
            }
        }

        dpdM /= 3;
        dpdC /= 3;

        return NULL;
    }

    errorOut computeReferenceHigherOrderStressPressure( const variableVector &referenceHigherOrderStress,
                                                        const variableVector &rightCauchyGreenDeformation,
                                                        variableVector &referenceHigherOrderPressure,
                                                        variableMatrix &dpdM, variableMatrix &dpdC,
                                                        variableMatrix &d2pdMdC ){
        /*!
         * Compute the pressure for a higher-order stress in the reference configuration.
         * $p_K = \frac{1}{3} C_{AB} M_{ABK}$
         *
         * Also compute the Jacobians
         * $\frac{ \partial p_K }{ \partial M_{NOP} } = \frac{1}{3} C_{NO} \delta_{KP}$
         * $\frac{ \partial p_K }{ \partial C_{NO} } = \frac{1}{3} M_{NOK}$
         * $\frac{ \partial^2 p_K}{ \partial M_{NOP} C_{QR} } = \frac{1}{3} \delta_{NQ} \delta_{OR} \delta_{KP}
         *
         * where $C_{AB}$ is the right Cauchy-Green deformation tensor and
         * M_{ABK} is the higher order stress tensor in the reference configuration.
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the
         *     reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation
         *     tensor.
         * :param variableVector &referenceHigherOrderPressure: The higher order pressure.
         * :param variableMatrix &dpdM: The Jacobian of the pressure w.r.t. the higher order stress.
         * :param variableMatrix &dpdC: The Jacobian of the pressure w.r.t. the right Cauchy-Green
         *     deformation tensor.
         * :param variableMatrix &d2pdMdC: The second order jacobian of the pressure w.r.t the 
         *     reference higher order stress and right Cauchy-Green deformation tensor. This Jacobian is organized
         *     [ K ][ NOPQR ]
         */

        //Assume 3D
        unsigned int dim = 3;

        errorOut error = computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation,
                                                                    referenceHigherOrderPressure, dpdM, dpdC );

        if ( error ){
            errorOut result = new errorNode( "computeReferenceHigherOrderStressPressure (second order jacobian)",
                                             "Error in computation of reference higher order pressure" );
            result->addNext( error );
            return result;
        }

        variableVector eye( dim * dim );
        vectorTools::eye( eye );

        d2pdMdC = variableMatrix( dim, variableVector( dim * dim * dim * dim * dim, 0 ) );

        for ( unsigned int K = 0; K < dim; K++ ){
            for ( unsigned int N = 0; N < dim; N++ ){
                for ( unsigned int O = 0; O < dim; O++ ){
                    for ( unsigned int P = 0; P < dim; P++ ){
                        for ( unsigned int Q = 0; Q < dim; Q++ ){
                            for ( unsigned int R = 0; R < dim; R++ ){
                                d2pdMdC[ K ][ dim * dim * dim * dim * N + dim * dim * dim * O + dim * dim * P + dim * Q + R ] = 
                                    eye[ dim * N + Q ] * eye[ dim * O + R ] * eye[ dim * K + P ] / 3;
                            }
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * dev ( M_{IJK} ) = M_{IJK} - ( 1 / 3 ) (C^{-1})_{IJ} C_{AB} M_{ABK}
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * :param variableVector &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        //Compute the pressure
        variableVector pressure;
        errorOut error = computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation, pressure );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricReferenceHigherOrderStress",
                                             "Error in computation of higher order pressure" );
            result->addNext( error );
            return result;
        }

        deviatoricReferenceHigherOrderStress = referenceHigherOrderStress;

        variableVector invRCG = vectorTools::inverse( rightCauchyGreenDeformation, dim, dim );
        
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    deviatoricReferenceHigherOrderStress[ dim * dim * I + dim * J + K ] -= invRCG[ dim * I + J ] * pressure[ K ];
                }
            }
        }

        return NULL;
    }

    errorOut computeDeviatoricReferenceHigherOrderStress( const variableVector &referenceHigherOrderStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress,
                                                          variableMatrix &dDeviatoricReferenceHigherOrderStressdRCG ){
        /*!
         * Compute the deviatoric part of the higher order stress in the reference configuration.
         *
         * dev ( M_{IJK} ) = M_{IJK} - \frac{1}{3} (C^{-1})_{IJ} C_{AB} M_{ABK} = M_{IJK} - ( C^{-1} )_{IJ} p_{K}
         *
         * Also compute Jacobians:
         * \frac{ \partial dev ( M_{IJK} ) }{ \partial M_{LMN} } = \delta_{IL} \delta_{JM} \delta_{KN} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial M_{LMN} }
         * \frac{ \partial dev ( M_{IJK} ) }{ \partial C_{LM} } = \left( (C^{-1})_{IL} (C^{-1})_{MJ} p_{K} - (C^{-1})_{IJ} \frac{ \partial p_K }{ \partial C_{LM} } \right)
         *
         * :param const variableVector &referenceHigherOrderStress: The higher order stress in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy-Green deformation tensor.
         * :param variableVector &deviatoricReferenceHigherOrderStress: The deviatoric part of the higher order tensor in the 
         *     reference configuration.
         */
        
        //Assume 3d
        unsigned int dim = 3;

        variableVector invRCG, pressure;
        variableMatrix dpdM, dpdC;

        //Compute the pressure
        errorOut error = computeReferenceHigherOrderStressPressure( referenceHigherOrderStress, rightCauchyGreenDeformation, pressure,
                                                                    dpdM, dpdC );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricReferenceHigherOrderStress (jacobian)",
                                             "Error in computation of higher order pressure" );
            result->addNext( error );
            return result;
        }

        deviatoricReferenceHigherOrderStress = referenceHigherOrderStress;

        invRCG = vectorTools::inverse( rightCauchyGreenDeformation, dim, dim );

        //Compute the deviatoric higher order stress        
        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    deviatoricReferenceHigherOrderStress[ dim * dim * I + dim * J + K ] -= invRCG[ dim * I + J ] * pressure[ K ];
                }
            }
        }

        //Compute the jacobians
        constantVector eye( dim * dim );
        vectorTools::eye( eye );

        dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress = variableMatrix( dim * dim * dim, variableVector( dim * dim * dim, 0 ) );
        dDeviatoricReferenceHigherOrderStressdRCG = variableMatrix( dim * dim * dim, variableVector( dim * dim, 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        for ( unsigned int M = 0; M < dim; M++ ){
                            for ( unsigned int N = 0; N < dim; N++ ){
                                dDeviatoricReferenceHigherOrderStressdReferenceHigherOrderStress[ dim * dim * I + dim * J + K ][ dim * dim * L + dim * M + N ]
                                    = eye[ dim * I + L ] * eye[ dim * J + M ] * eye[ dim * K + N ] - invRCG[ dim * I + J ] * dpdM[ K ][ dim * dim * L + dim * M + N ];
                            }

                            dDeviatoricReferenceHigherOrderStressdRCG[ dim * dim * I + dim * J + K ][ dim * L + M ]
                                = invRCG[ dim * I + L ] * invRCG[ dim * M + J ] * pressure[ K ] - invRCG[ dim * I + J ] * dpdC[ K ][ dim * L + M ];
                        }
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the current configuration.
         * \hat{s}_{ij} = s_{ij} - \frac{1}{3} s_{mm} \delta_{ij}
         *
         * :param const variableVector &secondOrderStress: The stress measure in the current configuration.
         * :param variableVector &deviatoricSecondOrderStress: The deviatoric part of the second order stress 
         *     measure.
         */

        //Assume 3d
        unsigned int dim = 3;

        variableVector eye( dim * dim );
        vectorTools::eye( eye );

        variableType trace = vectorTools::trace( secondOrderStress );

        deviatoricSecondOrderStress = secondOrderStress - trace * eye / 3.;

        return NULL;
    }

    errorOut computeDeviatoricSecondOrderStress( const variableVector &secondOrderStress,
                                                 variableVector &deviatoricSecondOrderStress,
                                                 variableMatrix &dDeviatoricStressdStress ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the current configuration.
         * \hat{s}_{ij} = s_{ij} - \frac{1}{3} s_{mm} \delta_{ij}
         *
         * Also return the jacobian
         * \frac{\partial \hat{s}_{ij}}{\partial s_{kl} } \delta_{ik} \delta_{jl} - \frac{1}{3} \delta_{ij} \delta_{kl}
         *
         * :param const variableVector &secondOrderStress: The stress measure in the current configuration.
         * :param variableVector &deviatoricSecondOrderStress: The deviatoric part of the second order stress 
         *     measure.
         * :param variableMatrix &dDeviatoricStressdStress: The jacobian of the deviatoric stress.
         */

        //Assume 3d
        unsigned int dim = 3;

        variableVector eye( dim * dim );
        vectorTools::eye( eye );

        errorOut error = computeDeviatoricSecondOrderStress( secondOrderStress, deviatoricSecondOrderStress );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricSecondOrderStress (jacobian)",
                                             "Error in the computation of the deviatoric stress" );
            result->addNext( error );
            return result;
        }

        dDeviatoricStressdStress = variableMatrix( deviatoricSecondOrderStress.size(), variableVector( secondOrderStress.size(), 0 ) );

        for ( unsigned int i = 0; i < dim; i++ ){
            for ( unsigned int j = 0; j < dim; j++ ){
                for ( unsigned int k = 0; k < dim; k++ ){
                    for ( unsigned int l = 0; l < dim; l++ ){
                        dDeviatoricStressdStress[ dim * i + j ][ dim * k + l ] = eye[ dim * i + k ] * eye[ dim * j + l ]
                                                                               - eye[ dim * i + j ] * eye[ dim * k + l ] / 3;
                    }
                }
            }
        }

        return NULL;
    }

    errorOut computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen,  variableType &pressure ){
        /*!
         * Compute the pressure part of a second order stress measure in the reference configuration.
         * p = \frac{1}{3} C_{IJ} S_{IJ}
         *
         * where C_{IJ} is the right Cauchy-Green deformation tensor and S_{IJ} is the stress measure.
         *
         * :param const variableVector &referenceStressMeasure: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreen: The right Cauchy-Green deformation tensor between the 
         *     current configuration and the reference configuration of the stress tensor.
         * :param variableType &pressure: The computed pressure.
         */

        if ( referenceStressMeasure.size() != rightCauchyGreen.size() ){
            return new errorNode( "computeReferenceSecondOrderStressPressure",
                                  "The stress measure and right Cauchy-Green deformation tensors aren't the same size" );
        }

        pressure = vectorTools::dot( referenceStressMeasure, rightCauchyGreen ) / 3;

        return NULL;
    }

    errorOut computeReferenceSecondOrderStressPressure( const variableVector &referenceStressMeasure,
                                                        const variableVector &rightCauchyGreen,  variableType &pressure,
                                                        variableVector &dpdStress, variableVector &dpdRCG ){
        /*!
         * Compute the pressure part of a second order stress measure in the reference configuration.
         * p = \frac{1}{3} C_{IJ} S_{IJ}
         *
         * where C_{IJ} is the right Cauchy-Green deformation tensor and S_{IJ} is the stress measure.
         *
         * :param const variableVector &referenceStressMeasure: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreen: The right Cauchy-Green deformation tensor between the 
         *     current configuration and the reference configuration of the stress tensor.
         * :param variableType &pressure: The computed pressure.
         * :param variableVector &dpdStress: The Jacobian of the pressure w.r.t. the stress.
         * :param variableVector &dpdRCG: The Jacobian of the pressure w.r.t. the right Cauchy-Green deformation tensor.
         */

        errorOut error = computeReferenceSecondOrderStressPressure( referenceStressMeasure, rightCauchyGreen, pressure );

        if ( error ){
            errorOut result = new errorNode( "computeReferenceSecondOrderStressPressure (jacobian)",
                                             "Error in computation of pressure in the reference configuration" );
            result->addNext( error );
            return result;
        }

        dpdStress = rightCauchyGreen / 3;
        dpdRCG    = referenceStressMeasure / 3;

        return NULL;
    }

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}
         *
         * :param const variableVector &secondOrderReferenceStress: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy Green Deformation tensor of the 
         *     deformation between configurations.
         * :param variableVector &deviatoricSecondOrderReferenceStrain: The deviatoric part of the second order 
         *     stress in the reference configuration.
         */

        variableVector invRCG;
        variableType pressure;
        return computeDeviatoricReferenceSecondOrderStress( secondOrderReferenceStress, rightCauchyGreenDeformation,
                                                            invRCG, pressure, deviatoricSecondOrderReferenceStress );

    }

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &invRCG, variableType &pressure,
                                                          variableVector &deviatoricSecondOrderReferenceStress ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}
         *
         * :param const variableVector &secondOrderReferenceStress: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy Green Deformation tensor of the 
         *     deformation between configurations.
         * :param variableVector &invRCG: The inverse of the right Cauchy-Green deformation tensor.
         * :param variableType &pressure: The pressure in the reference configuration.
         * :param variableVector &deviatoricSecondOrderReferenceStrain: The deviatoric part of the second order 
         *     stress in the reference configuration.
         */

        //Assume 3d
        unsigned int dim = 3;

        errorOut error = computeReferenceSecondOrderStressPressure( secondOrderReferenceStress, rightCauchyGreenDeformation, pressure );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricReferenceSecondOrderStress",
                                             "Error in computation of the reference pressure" );
            result->addNext( error );
            return result;
        }

        invRCG = vectorTools::inverse( rightCauchyGreenDeformation, dim, dim );

        deviatoricSecondOrderReferenceStress = secondOrderReferenceStress - pressure * invRCG;

        return NULL;
    }

    errorOut computeDeviatoricReferenceSecondOrderStress( const variableVector &secondOrderReferenceStress,
                                                          const variableVector &rightCauchyGreenDeformation,
                                                          variableVector &deviatoricSecondOrderReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdReferenceStress,
                                                          variableMatrix &dDeviatoricReferenceStressdRCG ){
        /*!
         * Compute the deviatoric part of a second order stress measure in the reference configuration.
         * \hat{S}_{IJ} = S_{IJ} - \frac{1}{3} C_{AB} S_{AB} (C^{-1})_{IJ}
         * 
         * Also compute the Jacobians.
         * \frac{ \partial \hat{S}_{IJ} }{ \partial S_{KL} } = \delta_{IK} \delta_{LJ} - \frac{1}{3} C_{KL} (C^{-1})_{IJ}
         * \frac{ \partial \hat{S}_{IJ} }{ \partial C_{KL} } = \frac{1}{3} \left( C_{AB} S_{AB} (C^{-1})_{IK} (C^{-1})_{LJ} - S_{KL} (C^{-1}_{IJ}) \right)
         *
         * :param const variableVector &secondOrderReferenceStress: The stress measure in the reference configuration.
         * :param const variableVector &rightCauchyGreenDeformation: The right Cauchy Green Deformation tensor of the 
         *     deformation between configurations.
         * :param variableVector &deviatoricSecondOrderReferenceStrain: The deviatoric part of the second order 
         *     stress in the reference configuration.
         * :param variableMatrix &dDeviatoricReferenceStressdReferenceStress: The jacobian w.r.t. the reference stress.
         * :param variableMatrix &dDeviatoricreferenceStressdRCG: The jacobian w.r.t. the right Cauchy Green deformation
         *     tensor.
         */

        //Assume 3d
        unsigned int dim = 3;
        constantVector eye ( dim * dim );
        vectorTools::eye( eye );

        variableVector invRCG = vectorTools::inverse( rightCauchyGreenDeformation, dim, dim );
        variableVector dpdStress, dpdRCG;
        variableType pressure;

        errorOut error = computeReferenceSecondOrderStressPressure( secondOrderReferenceStress, rightCauchyGreenDeformation, pressure,
                                                                    dpdStress, dpdRCG );

        if ( error ){
            errorOut result = new errorNode( "computeDeviatoricReferenceSecondOrderStress (jacobian)",
                                             "Error in computation of the reference pressure" );
            result->addNext( error );
            return result;
        }

        //Compute the deviatoric part of the reference stress
        deviatoricSecondOrderReferenceStress = secondOrderReferenceStress - pressure * invRCG;

        //Assemblet the Jacobian
        dDeviatoricReferenceStressdReferenceStress = variableMatrix( deviatoricSecondOrderReferenceStress.size(),
                                                                     variableVector( secondOrderReferenceStress.size(), 0 ) );

        dDeviatoricReferenceStressdRCG = variableMatrix( deviatoricSecondOrderReferenceStress.size(),
                                                         variableVector( rightCauchyGreenDeformation.size(), 0 ) );

        for ( unsigned int I = 0; I < dim; I++ ){
            for ( unsigned int J = 0; J < dim; J++ ){
                for ( unsigned int K = 0; K < dim; K++ ){
                    for ( unsigned int L = 0; L < dim; L++ ){
                        dDeviatoricReferenceStressdReferenceStress[ dim * I + J ][ dim * K + L ] = eye[ dim * I + K ] * eye[ dim * L + J ]
                            - dpdStress[ dim * K + L ] * invRCG[ dim * I + J ];

                        dDeviatoricReferenceStressdRCG[ dim * I + J ][ dim * K + L ] = ( pressure * invRCG[ dim * I + K ] * invRCG[ dim * L + J ] - dpdRCG[ dim * K + L ] * invRCG[ dim * I + J ] );
                    }
                }
            }
        }

        return NULL;
    }



}
