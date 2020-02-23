//Tests for constitutive_tools

#include<micromorphic_tools.h>
#include<sstream>
#include<fstream>
#include<iostream>

typedef micromorphicTools::constantType constantType;
typedef micromorphicTools::constantVector constantVector;
typedef micromorphicTools::constantMatrix constantMatrix;

typedef micromorphicTools::parameterType parameterType;
typedef micromorphicTools::parameterVector parameterVector;
typedef micromorphicTools::parameterMatrix parameterMatrix;

typedef micromorphicTools::variableType variableType;
typedef micromorphicTools::variableVector variableVector;
typedef micromorphicTools::variableMatrix variableMatrix;

typedef micromorphicTools::errorNode errorNode;
typedef micromorphicTools::errorOut errorOut;

struct cout_redirect{
    cout_redirect( std::streambuf * new_buffer)
        : old( std::cout.rdbuf( new_buffer ) )
    { }

    ~cout_redirect( ) {
        std::cout.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

struct cerr_redirect{
    cerr_redirect( std::streambuf * new_buffer)
        : old( std::cerr.rdbuf( new_buffer ) )
    { }

    ~cerr_redirect( ) {
        std::cerr.rdbuf( old );
    }

    private:
        std::streambuf * old;
};

int test_computePsi( std::ofstream &results ){
    /*!
     * Tests of the compute Psi function.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector F  = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
    variableVector Xi = { 9, 10, 11, 12, 13, 14, 15, 16, 17 };

    variableVector answer = { 126, 135, 144, 162, 174, 186, 198, 213, 228 };

    variableVector result;

    errorOut error = micromorphicTools::computePsi( F, Xi, result );

    if ( error ){
        error->print();
        results << "test_computePsi & False\n";
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_computePsi (test 1) & False\n";
        return 1;
    }

    //Test Jacobians
    variableVector resultJ;
    variableMatrix dPsidF, dPsidXi;

    error = micromorphicTools::computePsi( F, Xi, resultJ, dPsidF, dPsidXi );

    if ( error ){
        error->print();
        results << "test_computePsi & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_computePsi (test 2) & False\n";
        return 1;
    }

    //Test dPsidF
    
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < F.size(); i++ ){
        constantVector delta( F.size(), 0 );
        delta[i] = eps * fabs( F[i] ) + eps;

        error = micromorphicTools::computePsi( F + delta, Xi, resultJ );

        if ( error ){
            error->print();
            results << "test_computePsi & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dPsidF[j][i] ) ){
                results << "test_computePsi (test 3) & False\n";
                return 1;
            }
        }
    }

    for ( unsigned int i = 0; i < Xi.size(); i++ ){
        constantVector delta( Xi.size(), 0 );
        delta[i] = eps * fabs( Xi[i] ) + eps;

        error = micromorphicTools::computePsi( F, Xi + delta, resultJ );

        if ( error ){
            error->print();
            results << "test_computePsi & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dPsidXi[j][i] ) ){
                results << "test_computePsi (test 4) & False\n";
                return 1;
            }
        }
    }

    results << "test_computePsi & True\n";
    return 0;
}

int test_computeMicroStrain( std::ofstream &results ){
    /*!
     * Test the computation of the micro-strain
     *
     * :param std::ofstream &results: The output file
     */

    variableVector Psi = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    variableVector answer = { 0, 2, 3, 4, 4, 6, 7, 8, 8 };

    variableVector result;

    errorOut error = micromorphicTools::computeMicroStrain( Psi, result );

    if ( error ){
        error->print();
        results << "test_computeMicroStrain & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_computeMicroStrain (test 1) & False\n";
        return 1;
    }

    //Test Jacobians 
    
    variableVector resultJ;
    variableMatrix dMicroStraindPsi;
    error = micromorphicTools::computeMicroStrain( Psi, resultJ, dMicroStraindPsi );

    if ( error ){
        error->print();
        results << "test_computeMicroStrain & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_computeMicroStrain (test 2) & False\n";
        return 1;
    }

    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < Psi.size(); i++ ){
        constantVector delta( Psi.size(), 0 );
        delta[i] = eps * fabs( Psi[i] ) + eps;

        error = micromorphicTools::computeMicroStrain( Psi + delta, resultJ );

        if ( error ){
            error->print();
            results << "test_computeMicroStrain & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dMicroStraindPsi[j][i] ) ){
                results << "test_computeMicroStrain (test 3) & False\n";
                return 1;
            }
        }
    }

    results << "test_computeMicroStrain & True\n";
    return 0;
}

int test_pushForwardPK2Stress( std::ofstream &results ){
    /*!
     * Test the computation of the push-foward operation on the PK2 Stress.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector PK2Stress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    variableVector deformationGradient  = { -1.08831037, -0.66333427, -0.48239487,
                                            -0.904554  ,  1.28942848, -0.02156112,
                                            -0.08464824, -0.07730218,  0.86415668 };

    variableVector answer = { -10.75427056,   2.95576352,   4.7810659 ,
                                5.36943821,  -1.06947471,  -1.91553073,
                                7.58260457,  -1.61246489,  -2.82366599 };

    variableVector result;

    errorOut error = micromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient, result );

    if ( error ){
        error->print();
        results << "test_pushForwardPK2Stress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_pushForwardPK2Stress (test 1) & False\n";
        return 1;
    }

    //Test Jacobian
    variableVector resultJ;
    variableMatrix dCauchydPK2, dCauchydF;

    error = micromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient, resultJ,
                                                     dCauchydPK2, dCauchydF );

    if ( error ){
        error->print();
        results << "test_pushForwardPK2Stress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_pushForwardPK2Stress (test 2) & False\n";
        return 1;
    }

    //Test dCauchydPK2
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < PK2Stress.size(); i++ ){
        constantVector delta( PK2Stress.size(), 0 );
        delta[i] = eps * fabs( PK2Stress[i] ) + eps;

        error = micromorphicTools::pushForwardPK2Stress( PK2Stress + delta, deformationGradient, resultJ );

        if ( error ){
            error->print();
            results << "test_pushForwardPK2Stress & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dCauchydPK2[j][i] ) ){
                results << "test_pushForwardPK2Stress (test 3) & False\n";
                return 1;
            }
        }
    }

    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[i] = eps * fabs( deformationGradient[i] ) + eps;

        error = micromorphicTools::pushForwardPK2Stress( PK2Stress, deformationGradient + delta, resultJ );

        if ( error ){
            error->print();
            results << "test_pushForwardPK2Stress & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dCauchydF[j][i], 1e-5, 1e-5 ) ){
                results << "test_pushForwardPK2Stress (test 4) & False\n";
                return 1;
            }
        }
    }

    results << "test_pushForwardPK2Stress & True\n";
    return 0;
}

int test_pushForwardReferenceMicroStress( std::ofstream &results ){
    /*!
     * Test the computation of the push-foward operation on the reference micro-stress.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector referenceMicroStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    variableVector deformationGradient  = { -1.08831037, -0.66333427, -0.48239487,
                                            -0.904554  ,  1.28942848, -0.02156112,
                                            -0.08464824, -0.07730218,  0.86415668 };

    variableVector answer = { -10.75427056,   2.95576352,   4.7810659 ,
                                5.36943821,  -1.06947471,  -1.91553073,
                                7.58260457,  -1.61246489,  -2.82366599 };

    variableVector result;

    errorOut error = micromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, result );

    if ( error ){
        error->print();
        results << "test_pushForwardReferenceMicroStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_pushForwardReferenceMicroStress (test 1) & False\n";
        return 1;
    }

    //Test Jacobian
    variableVector resultJ;
    variableMatrix dsdS, dsdF;

    error = micromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient, resultJ,
                                                                dsdS, dsdF );

    if ( error ){
        error->print();
        results << "test_pushForwardReferenceMicroStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_pushForwardReferenceMicroStress (test 2) & False\n";
        return 1;
    }

    //Test dsdS
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < referenceMicroStress.size(); i++ ){
        constantVector delta( referenceMicroStress.size(), 0 );
        delta[i] = eps * fabs( referenceMicroStress[i] ) + eps;

        error = micromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress + delta, deformationGradient, resultJ );

        if ( error ){
            error->print();
            results << "test_pushForwardReferenceMicroStress & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dsdS[j][i] ) ){
                results << "test_pushForwardReferenceMicroStress (test 3) & False\n";
                return 1;
            }
        }
    }

    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[i] = eps * fabs( deformationGradient[i] ) + eps;

        error = micromorphicTools::pushForwardReferenceMicroStress( referenceMicroStress, deformationGradient + delta, resultJ );

        if ( error ){
            error->print();
            results << "test_pushForwardReferenceMicroStress & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dsdF[j][i], 1e-5, 1e-5 ) ){
                results << "test_pushForwardReferenceMicroStress (test 4) & False\n";
                return 1;
            }
        }
    }

    results << "test_pushForwardReferenceMicroStress & True\n";
    return 0;
}

int test_computeGamma( std::ofstream &results ){
    /*!
     * Test the computation of the deformation gradient Gamma
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector deformationGradient = { -1, -2, -3, -4, -5, -6, -7, -8, -9 };
    variableVector gradXi = { 1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector answer = { -174, -186, -198, -210, -222, -234, -246, -258, -270,
                              -204, -219, -234, -249, -264, -279, -294, -309, -324,
                              -234, -252, -270, -288, -306, -324, -342, -360, -378 };

    variableVector result;

    errorOut error = micromorphicTools::computeGamma( deformationGradient, gradXi, result );

    if ( error ){
        error->print();
        results << "test_computeGamma & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_computeGamma (test 1) & False\n";
        return 1;
    }

    //Test Jacobian
    variableVector resultJ;
    variableMatrix dGammadF, dGammadGradXi;

    error = micromorphicTools::computeGamma( deformationGradient, gradXi, resultJ, 
                                             dGammadF, dGammadGradXi );

    if ( error ){
        error->print();
        results << "test_computeGamma & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_computeGamma (test 2) & False\n";
        return 1;
    }

    //Test dGammadF
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[i] = eps * fabs( deformationGradient[i] ) + eps;

        error = micromorphicTools::computeGamma( deformationGradient + delta, gradXi, resultJ );

        if ( error ){
            error->print();
            results << "test_computeGamma & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dGammadF[j][i] ) ){
                results << "test_computeGamma (test 3) & False\n";
                return 1;
            }
        }
    }

    //Test dGammadGradXi
    for ( unsigned int i = 0; i < gradXi.size(); i++ ){
        constantVector delta( gradXi.size(), 0 );
        delta[i] = eps * fabs( gradXi[i] ) + eps;

        error = micromorphicTools::computeGamma( deformationGradient, gradXi + delta, resultJ );

        if ( error ){
            error->print();
            results << "test_computeGamma & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dGammadGradXi[j][i] ) ){
                results << "test_computeGamma (test 4) & False\n";
                return 1;
            }
        }
    }

    results << "test_computeGamma & True\n";
    return 0;
}

int test_pushForwardHigherOrderStress( std::ofstream &results){
    /*!
     * Tests for the push-forward operation of the higher order stress.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector referenceHigherOrderStress = {  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                                  10, 11, 12, 13, 14, 15, 16, 17, 18,
                                                  19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector deformationGradient = { 0.29524861, -0.25221581, -1.60534711,
                                          -0.08817703,  0.28447808,  0.06451703,
                                           0.39849454,  0.38681512,  0.87840084 };

    variableVector microDeformation = { -0.25781969, -0.39826899, -0.79493259,
                                         0.38104724, -0.00830511, -0.51985409,
                                        -0.36415661, -0.6871168 ,  0.54018665 };

    variableVector answer = { -370.21000924,  -44.9887908 , -120.76625915,   57.76488049,
                                 7.10106323,   18.73840733,  356.3462823 ,   44.06705478,
                               115.25802557,   49.68640146,    6.28202573,   15.89296064,
                                -7.62049271,   -0.98037632,   -2.41571071,  -46.58548828,
                                -6.04841312,  -14.69638769,  280.56462547,   36.38392302,
                                88.5657901 ,  -42.53702297,   -5.63795901,  -13.27041478,
                              -258.4236358 ,  -34.6543976 ,  -80.10152498 };

    variableVector result;

    errorOut error = micromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress, deformationGradient,
                                                                      microDeformation, result );

    if ( error ){
        error->print();
        results << "test_pushForwardHigherOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_pushForwardHigherOrderStress (test 1) & False\n";
        return 1;
    }

    //Test Jacobian
    variableVector resultJ;
    variableMatrix dHigherOrderStressdReferenceHigherOrderStress;
    variableMatrix dHigherOrderStressdDeformationGradient;
    variableMatrix dHigherOrderStressdMicroDeformation;

    error = micromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress,
                                                             deformationGradient, microDeformation,
                                                             resultJ,
                                                             dHigherOrderStressdReferenceHigherOrderStress,
                                                             dHigherOrderStressdDeformationGradient,
                                                             dHigherOrderStressdMicroDeformation );

    if ( error ){
        error->print();
        results << "test_pushForwardHigherOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_pushForwardHigherOrderStress (test 2) & False\n";
        return 1;
    }

    //Test dHigherOrderStressdReferenceHigherOrderStress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < referenceHigherOrderStress.size(); i++ ){
        constantVector delta( referenceHigherOrderStress.size(), 0 );
        delta[i] = eps * fabs( referenceHigherOrderStress[i] ) + eps;

        error = micromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress + delta,
                                                                 deformationGradient, microDeformation,
                                                                 resultJ );

        if ( error ){
            error->print();
            results << "test_pushForwardHigherOrderStress & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dHigherOrderStressdReferenceHigherOrderStress[j][i] ) ){
                results << "test_pushForwardHigherOrderStress (test 3) & False\n";
                return 1;
            }
        }
    }

    //Test dHigherOrderStressdDeformationGradient
    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[i] = eps * fabs( deformationGradient[i] ) + eps;

        error = micromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress,
                                                                 deformationGradient + delta, microDeformation,
                                                                 resultJ );

        if ( error ){
            error->print();
            results << "test_pushForwardHigherOrderStress & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dHigherOrderStressdDeformationGradient[j][i], 1e-5, 1e-5 ) ){
                results << "test_pushForwardHigherOrderStress (test 4) & False\n";
                return 1;
            }
        }
    }

    //Test dHigherOrderStressdMicroDeformation
    for ( unsigned int i = 0; i < microDeformation.size(); i++ ){
        constantVector delta( microDeformation.size(), 0 );
        delta[i] = eps * fabs( microDeformation[i] ) + eps;

        error = micromorphicTools::pushForwardHigherOrderStress( referenceHigherOrderStress,
                                                                 deformationGradient, microDeformation + delta,
                                                                 resultJ );

        if ( error ){
            error->print();
            results << "test_pushForwardHigherOrderStress & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dHigherOrderStressdMicroDeformation[j][i] ) ){
                results << "test_pushForwardHigherOrderStress (test 5) & False\n";
                return 1;
            }
        }
    }

    results << "test_pushForwardHigherOrderStress & True\n";
    return 0;
}

int test_computeDeviatoricHigherOrderStress( std::ofstream &results ){
    /*!
     * Test the computation of the deviatoric higher order stress.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector higherOrderStress = {  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                         10, 11, 12, 13, 14, 15, 16, 17, 18,
                                         19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector answer = { -12., -12., -12.,   4.,   5.,   6.,   7.,   8.,   9.,
                               10.,  11.,  12.,   0.,   0.,   0.,  16.,  17.,  18.,
                               19.,  20.,  21.,  22.,  23.,  24.,  12.,  12.,  12. };

    variableVector result;

    errorOut error = micromorphicTools::computeDeviatoricHigherOrderStress( higherOrderStress, result );

    if ( error ){
        results << "test_computeDeviatoricHigherOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_computeDeviatoricHigherOrderStress (test 1) & False\n";
        return 1;
    }

    //Test Jacobians
    
    variableVector resultJ;
    variableMatrix dDeviatoricHigherOrderStressdHigherOrderStress;

    error = micromorphicTools::computeDeviatoricHigherOrderStress( higherOrderStress, resultJ,
                                                                   dDeviatoricHigherOrderStressdHigherOrderStress );

    if ( error ){
        results << "test_computeDeviatoricHigherOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, resultJ ) ){
        results << "test_computeDeviatoricHigherOrderStress (test 2) & False\n";
        return 1;
    }

    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < higherOrderStress.size(); i++ ){
        constantVector delta( higherOrderStress.size(), 0 );
        delta[i] = eps * fabs( higherOrderStress[i] ) + eps;

        error = micromorphicTools::computeDeviatoricHigherOrderStress( higherOrderStress + delta, resultJ );

        if ( error ){
            error->print();
            results << "test_computeDeviatoricHigherOrderStress & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDeviatoricHigherOrderStressdHigherOrderStress[j][i] ) ){
                results << "test_pushForwardHigherOrderStress (test 3) & False\n";
                return 1;
            }
        }
    }


    results << "test_computeDeviatoricHigherOrderStress & True\n";
    return 0;
}

int test_computeDeviatoricReferenceHigherOrderStress( std::ofstream &results ){
    /*!
     * Test the computation of the deviatoric part of the reference higher order stress.
     *
     * :param std::ofstream &results: The output file
     */

    variableVector M = { 0.80732114,  0.79202055,  0.17990022,  0.97454675,  0.703207  ,
                        -0.58236697,  0.53324571, -0.93438873, -0.40650796,  0.14071918,
                         0.66933708, -0.67854069, -0.30317772, -0.93821882,  0.97270622,
                         0.00295302, -0.12441126,  0.30539971, -0.0580227 ,  0.89696105,
                         0.17567709, -0.9592962 ,  0.63535407,  0.95437804, -0.64531877,
                         0.69978907,  0.81327586 };

    variableVector C = { 0.3991656 ,  0.43459435,  0.15398811,
                        -0.20239202, -0.50763359,  0.04756988,
                        -0.0573016 , -0.95939895,  0.2693173 };

    variableVector answer = { 5.27633771,   1.68477127,  -3.21111472,  13.96704974,
                              3.29864744, -10.4408607 ,  -4.31689678,  -1.90327486,
                              3.27369895,  -2.4001685 ,   0.161758  ,   1.24944234,
                             -6.01118653,  -2.07847616,   5.30384773,   2.46397506,
                              0.3672135 ,  -1.56198261,  -8.15866529,  -0.72125951,
                              6.32230906, -18.52878201,  -2.87440491,  14.2858097 ,
                              4.98150383,   1.82382829,  -3.45626291 };

    variableVector result;

    errorOut error = micromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C, result );

    if ( error ){
        results << "test_computeDeviatoricReferenceHigherOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_computeDeviatoricReferenceHigherOrderStress (test 1) & False\n";
        return 1;
    }

    variableVector resultJ;
    variableMatrix dDevMdM, dDevMdC;

    error = micromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C, resultJ, dDevMdM, dDevMdC );

    if ( error ){
        results << "test_computeDeviatoricReferenceHigherOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_computeDeviatoricReferenceHigherOrderStress (test 2) & False\n";
        return 1;
    }

    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < M.size(); i++ ){
        constantVector delta( M.size(), 0 );
        delta[i] = eps * fabs( M[i] ) + eps;

        error = micromorphicTools::computeDeviatoricReferenceHigherOrderStress( M + delta, C, resultJ );

        if ( error ){
            error->print();
            results << "test_computeDeviatoricReferenceHigherOrderStress & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevMdM[j][i] ) ){
                results << "test_pushForwardHigherOrderStress (test 3) & False\n";
                return 1;
            }
        }
    }

    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        error = micromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C + delta, resultJ );

        if ( error ){
            error->print();
            results << "test_computeDeviatoricReferenceHigherOrderStress & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevMdC[j][i], 1e-4, 1e-5 ) ){
                results << "test_computeDeviatoricReferenceHigherOrderStress (test 4) & False\n";
                return 1;
            }
        }
    }

    results << "test_computeDeviatoricReferenceHigherOrderStress & True\n";
    return 0;
}

int main(){
    /*!
    The main loop which runs the tests defined in the 
    accompanying functions. Each function should output
    the function name followed by & followed by True or False 
    if the test passes or fails respectively.
    */

    //Open the results file
    std::ofstream results;
    results.open("results.tex");

    //Run the tests
    test_computePsi( results );
    test_computeGamma( results );
    test_computeMicroStrain( results );
    test_pushForwardPK2Stress( results );
    test_pushForwardReferenceMicroStress( results );
    test_pushForwardHigherOrderStress( results );
    test_computeDeviatoricHigherOrderStress( results );
    test_computeDeviatoricReferenceHigherOrderStress( results );

    //Close the results file
    results.close();

    return 0;
}
