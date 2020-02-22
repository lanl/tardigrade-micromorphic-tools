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

}
