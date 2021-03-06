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
    variableVector Chi = { 9, 10, 11, 12, 13, 14, 15, 16, 17 };

    variableVector answer = { 126, 135, 144, 162, 174, 186, 198, 213, 228 };

    variableVector result;

    errorOut error = micromorphicTools::computePsi( F, Chi, result );

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
    variableMatrix dPsidF, dPsidChi;

    error = micromorphicTools::computePsi( F, Chi, resultJ, dPsidF, dPsidChi );

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

        error = micromorphicTools::computePsi( F + delta, Chi, resultJ );

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

    for ( unsigned int i = 0; i < Chi.size(); i++ ){
        constantVector delta( Chi.size(), 0 );
        delta[i] = eps * fabs( Chi[i] ) + eps;

        error = micromorphicTools::computePsi( F, Chi + delta, resultJ );

        if ( error ){
            error->print();
            results << "test_computePsi & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dPsidChi[j][i] ) ){
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
    variableVector gradChi = { 1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector answer = { -174, -186, -198, -210, -222, -234, -246, -258, -270,
                              -204, -219, -234, -249, -264, -279, -294, -309, -324,
                              -234, -252, -270, -288, -306, -324, -342, -360, -378 };

    variableVector result;

    errorOut error = micromorphicTools::computeGamma( deformationGradient, gradChi, result );

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
    variableMatrix dGammadF, dGammadGradChi;

    error = micromorphicTools::computeGamma( deformationGradient, gradChi, resultJ, 
                                             dGammadF, dGammadGradChi );

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

        error = micromorphicTools::computeGamma( deformationGradient + delta, gradChi, resultJ );

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

    //Test dGammadGradChi
    for ( unsigned int i = 0; i < gradChi.size(); i++ ){
        constantVector delta( gradChi.size(), 0 );
        delta[i] = eps * fabs( gradChi[i] ) + eps;

        error = micromorphicTools::computeGamma( deformationGradient, gradChi + delta, resultJ );

        if ( error ){
            error->print();
            results << "test_computeGamma & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dGammadGradChi[j][i] ) ){
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

    variableVector pressure;
    variableMatrix dpdM, dpdC, d2pdMdC;
    micromorphicTools::computeReferenceHigherOrderStressPressure( M, C, pressure, dpdM, dpdC, d2pdMdC );

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

    error = micromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C, pressure, result );

    if ( error ){
        results << "test_computeDeviatoricReferenceHigherOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_computeDeviatoricReferenceHigherOrderStress (test 2) & False\n";
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
        results << "test_computeDeviatoricReferenceHigherOrderStress (test 3) & False\n";
        return 1;
    }

    variableVector resultJP;
    variableMatrix dDevMdMP, dDevMdCP;

    error = micromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C, pressure, dpdM, dpdC, resultJP, dDevMdMP, dDevMdCP );

    if ( error ){
        results << "test_computeDeviatoricReferenceHigherOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJP, answer ) ){
        results << "test_computeDeviatoricReferenceHigherOrderStress (test 4) & False\n";
        return 1;
    }

    variableVector resultJ2;
    variableMatrix dDevMdMJ2, dDevMdCJ2, d2DevMdMdC;

    error = micromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C, resultJ2, dDevMdMJ2, dDevMdCJ2, d2DevMdMdC );

    if ( error ){
        results << "test_computeDeviatoricReferenceHigherOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_computeDeviatoricReferenceHigherOrderStress (test 5) & False\n";
        return 1;
    }

    variableVector resultJ2P;
    variableMatrix dDevMdMJ2P, dDevMdCJ2P, d2DevMdMdCP;

    error = micromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C, pressure, dpdM, dpdC, d2pdMdC,
                                                                            resultJ2P, dDevMdMJ2P, dDevMdCJ2P, d2DevMdMdCP );

    if ( error ){
        results << "test_computeDeviatoricReferenceHigherOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ2P, answer ) ){
        results << "test_computeDeviatoricReferenceHigherOrderStress (test 6) & False\n";
        return 1;
    }


    //Test dDevMdM
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
                results << "test_computeDeviatoricReferenceHigherOrderStress (test 7) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevMdMP[j][i] ) ){
                results << "test_computeDeviatoricReferenceHigherOrderStress (test 8) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevMdMJ2[j][i] ) ){
                results << "test_computeDeviatoricReferenceHigherOrderStress (test 9) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevMdMJ2P[j][i] ) ){
                results << "test_computeDeviatoricReferenceHigherOrderStress (test 10) & False\n";
                return 1;
            }
        }
    }

    //Test dDevMdC
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
                results << "test_computeDeviatoricReferenceHigherOrderStress (test 11) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevMdCP[j][i], 1e-4, 1e-5 ) ){
                results << "test_computeDeviatoricReferenceHigherOrderStress (test 12) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevMdCJ2[j][i], 1e-4, 1e-5 ) ){
                results << "test_computeDeviatoricReferenceHigherOrderStress (test 13) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevMdCJ2P[j][i], 1e-4, 1e-5 ) ){
                results << "test_computeDeviatoricReferenceHigherOrderStress (test 14) & False\n";
                return 1;
            }
        }
    }

    //Test d2DevMdMdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        error = micromorphicTools::computeDeviatoricReferenceHigherOrderStress( M, C + delta, resultJ2, dDevMdMJ2, dDevMdCJ2 );

        if ( error ){
            error->print();
            results << "test_computeDeviatoricReferenceHigherOrderStress & False\n";
            return 1;
        }

        constantMatrix gradCol = ( dDevMdMJ2 - dDevMdM ) / delta[i];

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        for ( unsigned int n = 0; n < 3; n++ ){
                            for ( unsigned int o = 0; o < 3; o++ ){
                                if ( !vectorTools::fuzzyEquals( gradCol[ 9 * j + 3 * k + l ][ 9 * m + 3 * n + o ],
                                    d2DevMdMdC[ 9 * j + 3 * k + l ][ 81 * m + 27 * n + 9 * o + 3 * (int)( i / 3 ) + i % 3 ], 1e-4, 1e-4 ) ){

                                    results << "test_computeDeviatoricReferenceHigherOrderStress (test 15) & False\n";
                                    return 1;
                                }
                                if ( !vectorTools::fuzzyEquals( gradCol[ 9 * j + 3 * k + l ][ 9 * m + 3 * n + o ],
                                    d2DevMdMdCP[ 9 * j + 3 * k + l ][ 81 * m + 27 * n + 9 * o + 3 * (int)( i / 3 ) + i % 3 ], 1e-4, 1e-4 ) ){

                                    results << "test_computeDeviatoricReferenceHigherOrderStress (test 16) & False\n";
                                    return 1;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    results << "test_computeDeviatoricReferenceHigherOrderStress & True\n";
    return 0;
}

int test_computeDeviatoricSecondOrderStress( std::ofstream &results ){
    /*!
     * Test the computation the deviatoric part of a second order 
     * stress tensor in the current configuration.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector stress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    variableVector answer = { -4.,  2.,  3.,  4.,  0.,  6.,  7.,  8.,  4. };

    variableVector result;

    errorOut error = micromorphicTools::computeDeviatoricSecondOrderStress( stress, result );

    if ( error ){
        results << "test_computeDeviatoricSecondOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_computeDeviatoricSecondOrderStress (test 1) & False\n";
        return 0;
    }

    variableVector resultJ;
    variableMatrix jacobian;

    error = micromorphicTools::computeDeviatoricSecondOrderStress( stress, resultJ, jacobian );

    if ( error ){
        results << "test_computeDeviatoricSecondOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, resultJ ) ){
        results << "test_computeDeviatoricSecondOrderStress (test 2) & False\n";
        return 0;
    }

    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < stress.size(); i++ ){
        constantVector delta( stress.size(), 0 );
        delta[i] = eps * fabs( stress[i] ) + eps;

        error = micromorphicTools::computeDeviatoricSecondOrderStress( stress + delta, resultJ );

        if ( error ){
            results << "test_computeDeviatoricSecondOrderStress & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], jacobian[j][i] ) ){
                results << "test_computeDeviatoricSecondOrderStress (test 3) & False\n";
                return 1;
            }
        }
    }

    results << "test_computeDeviatoricSecondOrderStress & True\n";
    return 0;
}

int test_computeReferenceSecondOrderStressPressure( std::ofstream &results ){
    /*!
     * Test the computation of the pressure of a second order stress in the 
     * reference configuration.
     *
     * :param std::ofstream &rsults: The output file.
     */

    variableVector S = { 0.77189588, -0.84417528,  0.95929231,
                        -0.50465708,  0.50576944,  0.05335127,
                         0.81510751,  0.76814059, -0.82146208 };

    variableVector C = { 0.03468919, -0.31275742, -0.57541261,
                        -0.27865312, -0.45844965,  0.52325004,
                        -0.0439162 , -0.80201065, -0.44921044 };

    variableType answer = -0.20245462701026676;

    variableType result;

    errorOut error = micromorphicTools::computeReferenceSecondOrderStressPressure( S, C, result );

    if ( error ){
        error->print();
        results << "test_computeReferenceSecondOrderStressPressure & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_computeReferenceSecondOrderStressPressure (test 1) & False\n";
        return 1;
    }

    //Test the Jacobians

    variableType resultJ;
    variableVector dpdS, dpdC;

    error = micromorphicTools::computeReferenceSecondOrderStressPressure( S, C, resultJ, dpdS, dpdC );

    if ( error ){
        error->print();
        results << "test_computeReferenceSecondOrderStressPressure & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_computeReferenceSecondOrderStressPressure (test 2) & False\n";
        return 1;
    }

    variableType resultJ2;
    variableVector dpdSJ2, dpdCJ2;
    variableMatrix d2pdSdC;

    error = micromorphicTools::computeReferenceSecondOrderStressPressure( S, C, resultJ2, dpdSJ2, dpdCJ2, d2pdSdC );

    if ( error ){
        error->print();
        results << "test_computeReferenceSecondOrderStressPressure & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ2, answer ) ){
        results << "test_computeReferenceSecondOrderStressPressure (test 3) & False\n";
        return 1;
    }


    //Test dpdS
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < S.size(); i++ ){
        constantVector delta( S.size(), 0 );
        delta[i] = eps * fabs( S[i] ) + eps;

        error = micromorphicTools::computeReferenceSecondOrderStressPressure( S + delta, C, resultJ );

        if ( error ){
            error->print();
            results << "test_computeReferenceSecondOrderStressPressure & False\n";
            return 1;
        }

        constantType gradCol = ( resultJ - result ) / delta[i];

        if ( !vectorTools::fuzzyEquals( gradCol, dpdS[i] ) ){
            results << "test_computeReferenceSecondOrderStressPressure (test 4) & False\n";
            return 1;
        }

        if ( !vectorTools::fuzzyEquals( gradCol, dpdSJ2[i] ) ){
            results << "test_computeReferenceSecondOrderStressPressure (test 5) & False\n";
            return 1;
        }

    }

    //Test dpdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        error = micromorphicTools::computeReferenceSecondOrderStressPressure( S, C + delta, resultJ );

        if ( error ){
            error->print();
            results << "test_computeReferenceSecondOrderStressPressure & False\n";
            return 1;
        }

        constantType gradCol = ( resultJ - result ) / delta[i];

        if ( !vectorTools::fuzzyEquals( gradCol, dpdC[i] ) ){
            results << "test_computeReferenceSecondOrderStressPressure (test 6) & False\n";
            return 1;
        }

        if ( !vectorTools::fuzzyEquals( gradCol, dpdCJ2[i] ) ){
            results << "test_computeReferenceSecondOrderStressPressure (test 7) & False\n";
            return 1;
        }
    }

    //Test d2pdSdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        error = micromorphicTools::computeReferenceSecondOrderStressPressure( S, C + delta, resultJ, dpdSJ2, dpdCJ2 );

        if ( error ){
            error->print();
            results << "test_computeReferenceSecondOrderStressPressure & False\n";
            return 1;
        }

        constantVector gradCol = ( dpdSJ2 - dpdS ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], d2pdSdC[j][i] ) ){
                results << "test_computeReferenceSecondOrderStressPressure (test 8) & False\n";
                return 1;
            }
        }
    }

    results << "test_computeReferenceSecondOrderStressPressure & True\n";
    return 0;
}

int test_computeDeviatoricReferenceSecondOrderStress( std::ofstream &results ){
    /*!
     * Test the computation of the deviatoric part of the second higher order stress.
     *
     * :param std::ofstream &results: The output file
     */

    variableVector S = { 0.77189588, -0.84417528,  0.95929231,
                        -0.50465708,  0.50576944,  0.05335127,
                         0.81510751,  0.76814059, -0.82146208 };

    variableVector C = { 0.03468919, -0.31275742, -0.57541261,
                        -0.27865312, -0.45844965,  0.52325004,
                        -0.0439162 , -0.80201065, -0.44921044 };

    variableVector answer = { -1.81433044, -2.17117481,  2.726381  ,
                               0.10781401,  0.6746562 , -0.53446584,
                              -0.0255477 ,  0.59634539, -0.39543207 };

    variableType pressure;
    variableVector dpdS, dpdC;
    variableMatrix d2pdSdC;

    errorOut error = micromorphicTools::computeReferenceSecondOrderStressPressure( S, C, pressure, dpdS, dpdC, d2pdSdC );

    if ( error ){
        results << "test_computeDeviatoricReferenceSecondOrderStress & False\n";
        return 1;
    }

    variableVector result;

    error = micromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C, result );

    if ( error ){
        results << "test_computeDeviatoricReferenceSecondOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_computeDeviatoricReferenceSecondOrderStress (test 1) & False\n";
        return 1;
    }

    error = micromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C, pressure, result );

    if ( error ){
        results << "test_computeDeviatoricReferenceSecondOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_computeDeviatoricReferenceSecondOrderStress (test 2) & False\n";
        return 1;
    }

    //Test the jacobians
    variableVector resultJ;
    variableMatrix dDevSdS, dDevSdC;

    error = micromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C, resultJ, dDevSdS, dDevSdC );

    if ( error ){
        results << "test_computeDeviatoricReferenceSecondOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_computeDeviatoricReferenceSecondOrderStress (test 3) & False\n";
        return 1;
    }

    variableVector resultJP;
    variableMatrix dDevSdSJP, dDevSdCJP;

    error = micromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C, pressure, dpdS, dpdC, 
                                                                            resultJP, dDevSdSJP, dDevSdCJP );

    if ( error ){
        results << "test_computeDeviatoricReferenceSecondOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJP, answer ) ){
        results << "test_computeDeviatoricReferenceSecondOrderStress (test 4) & False\n";
        return 1;
    }

    variableVector resultJ2;
    variableMatrix dDevSdSJ2, dDevSdCJ2, d2DevSdSdC;

    error = micromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C, resultJ2, dDevSdSJ2, dDevSdCJ2, d2DevSdSdC );

    if ( error ){
        results << "test_computeDeviatoricReferenceSecondOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ2, answer ) ){
        results << "test_computeDeviatoricReferenceSecondOrderStress (test 5) & False\n";
        return 1;
    }

    variableVector resultJ2P;
    variableMatrix dDevSdSJ2P, dDevSdCJ2P, d2DevSdSdCP;

    error = micromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C, pressure, dpdS, dpdC, d2pdSdC,
                                                                            resultJ2P, dDevSdSJ2P, dDevSdCJ2P, d2DevSdSdCP );

    if ( error ){
        results << "test_computeDeviatoricReferenceSecondOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ2P, answer ) ){
        results << "test_computeDeviatoricReferenceSecondOrderStress (test 6) & False\n";
        return 1;
    }

    //Test of dpdS
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < S.size(); i++ ){
        constantVector delta( S.size(), 0 );
        delta[i] = eps * fabs( S[i] ) + eps;

        error = micromorphicTools::computeDeviatoricReferenceSecondOrderStress( S + delta, C, resultJ );

        if ( error ){
            error->print();
            results << "test_computeDeviatoricReferenceSecondOrderStress & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevSdS[j][i] ) ){
                results << "test_computeDeviatoricReferenceSecondOrderStress (test 7) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevSdSJP[j][i] ) ){
                results << "test_computeDeviatoricReferenceSecondOrderStress (test 8) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevSdSJ2[j][i] ) ){
                results << "test_computeDeviatoricReferenceSecondOrderStress (test 9) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevSdSJ2P[j][i] ) ){
                results << "test_computeDeviatoricReferenceSecondOrderStress (test 10) & False\n";
                return 1;
            }
        }
    }

    //Test of dpdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        error = micromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C + delta, resultJ );

        if ( error ){
            error->print();
            results << "test_computeDeviatoricReferenceSecondOrderStress & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevSdC[j][i], 1e-4, 1e-5 ) ){
                results << "test_computeDeviatoricReferenceSecondOrderStress (test 11) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevSdCJP[j][i], 1e-4, 1e-5 ) ){
                results << "test_computeDeviatoricReferenceSecondOrderStress (test 12) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevSdCJ2[j][i], 1e-4, 1e-5 ) ){
                results << "test_computeDeviatoricReferenceSecondOrderStress (test 13) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){

            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevSdCJ2P[j][i], 1e-4, 1e-5 ) ){
                results << "test_computeDeviatoricReferenceSecondOrderStress (test 14) & False\n";
                return 1;
            }
        }
    }

    //Test of d2pdSdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        error = micromorphicTools::computeDeviatoricReferenceSecondOrderStress( S, C + delta, resultJ, dDevSdSJ2, dDevSdCJ2 );

        if ( error ){
            error->print();
            results << "test_computeDeviatoricReferenceSecondOrderStress & False\n";
            return 1;
        }

        constantMatrix gradCol = ( dDevSdSJ2 - dDevSdS ) / delta[i];

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        if ( !vectorTools::fuzzyEquals( gradCol[ 3 * j + k ][ 3 * l + m ],
                                                        d2DevSdSdC[ 3 * j + k ][ 27 * l + 9 * m + 3 * ( int )( i / 3 ) + i % 3 ],
                                                        1e-4, 1e-5 ) ){
                            results << "test_computeDeviatoricReferenceSecondOrderStress (test 15) & False\n";
                            return 1;
                        }
                        if ( !vectorTools::fuzzyEquals( gradCol[ 3 * j + k ][ 3 * l + m ],
                                                        d2DevSdSdCP[ 3 * j + k ][ 27 * l + 9 * m + 3 * ( int )( i / 3 ) + i % 3 ],
                                                        1e-4, 1e-5 ) ){
                            results << "test_computeDeviatoricReferenceSecondOrderStress (test 16) & False\n";
                            return 1;
                        }
                    }
                }
            }
        }
    }

    results << "test_computeDeviatoricReferenceSecondOrderStress & True\n";
    return 0;
}

int test_computeSecondOrderReferenceStressDecomposition( std::ofstream &results ){
    /*!
     * Test the computation of the deviatoric - volumetric (pressure) stress 
     * decomposition.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector S = { 0.77189588, -0.84417528,  0.95929231,
                        -0.50465708,  0.50576944,  0.05335127,
                         0.81510751,  0.76814059, -0.82146208 };

    variableVector C = { 0.03468919, -0.31275742, -0.57541261,
                        -0.27865312, -0.45844965,  0.52325004,
                        -0.0439162 , -0.80201065, -0.44921044 };

    variableVector deviatoricAnswer = { -1.81433044, -2.17117481,  2.726381  ,
                                         0.10781401,  0.6746562 , -0.53446584,
                                        -0.0255477 ,  0.59634539, -0.39543207 };

    variableType pressureAnswer = -0.20245462701026676;

    variableType pressureResult;
    variableVector deviatoricResult;

    errorOut error = micromorphicTools::computeSecondOrderReferenceStressDecomposition( S, C, deviatoricResult, pressureResult );

    if ( error ){
        error->print();
        results << "test_computeSecondOrderReferenceStressDecomposition & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( pressureResult, pressureAnswer ) ){
        results << "test_computeSecondOrderReferenceStressDecomposition (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( deviatoricResult, deviatoricAnswer ) ){
        results << "test_computeSecondOrderReferenceStressDecomposition (test 2) & False\n";
        return 1;
    }

    //Test the Jacobians
    variableType pressureResultJ;
    variableVector deviatoricResultJ;

    variableVector dPressuredStress, dPressuredRCG;
    variableMatrix dDevStressdStress, dDevStressdRCG;

    error = micromorphicTools::computeSecondOrderReferenceStressDecomposition( S, C, deviatoricResultJ, pressureResultJ,
                                                                               dDevStressdStress, dDevStressdRCG,
                                                                               dPressuredStress, dPressuredRCG );

    if ( error ){
        error->print();
        results << "test_computeSecondOrderReferenceStressDecomposition & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( pressureResultJ, pressureAnswer ) ){
        results << "test_computeSecondOrderReferenceStressDecomposition (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( deviatoricResultJ, deviatoricAnswer ) ){
        results << "test_computeSecondOrderReferenceStressDecomposition (test 4) & False\n";
        return 1;
    }

    variableType pressureResultJ2;
    variableVector deviatoricResultJ2;

    variableVector dPressuredStressJ2, dPressuredRCGJ2;
    variableMatrix dDevStressdStressJ2, dDevStressdRCGJ2;

    variableMatrix d2DevStressdStressdRCGJ2, d2PressuredStressdRCGJ2;

    error = micromorphicTools::computeSecondOrderReferenceStressDecomposition( S, C, deviatoricResultJ2, pressureResultJ2,
                                                                               dDevStressdStressJ2, dDevStressdRCGJ2,
                                                                               dPressuredStressJ2, dPressuredRCGJ2,
                                                                               d2DevStressdStressdRCGJ2, d2PressuredStressdRCGJ2 );

    if ( error ){
        error->print();
        results << "test_computeSecondOrderReferenceStressDecomposition & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( pressureResultJ2, pressureAnswer ) ){
        results << "test_computeSecondOrderReferenceStressDecomposition (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( deviatoricResultJ2, deviatoricAnswer ) ){
        results << "test_computeSecondOrderReferenceStressDecomposition (test 6) & False\n";
        return 1;
    }

    //Test deriatives w.r.t. the stress.
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < S.size(); i++ ){
        constantVector delta( S.size(), 0 );
        delta[i] = eps * fabs( S[i] ) + eps;

        error = micromorphicTools::computeSecondOrderReferenceStressDecomposition( S + delta, C, deviatoricResultJ, pressureResultJ );

        if ( error ){
            error->print();
            results << "test_computeSecondOrderReferenceStressDecomposition & False\n";
            return 1;
        }

        constantVector gradCol = ( deviatoricResultJ - deviatoricResult ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevStressdStress[j][i] ) ){
                results << "test_computeSecondOrderReferenceStressDecomposition (test 7) & False\n";
                return 1;
            }
            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevStressdStressJ2[j][i] ) ){
                results << "test_computeSecondOrderReferenceStressDecomposition (test 8) & False\n";
                return 1;
            }
        }

        constantType gradScalar = ( pressureResultJ - pressureResult ) / delta[i];

        if ( !vectorTools::fuzzyEquals( gradScalar, dPressuredStress[i] ) ){
            results << "test_computeSecondOrderReferenceStressDecomposition (test 9) & False\n";
            return 1;
        }

        if ( !vectorTools::fuzzyEquals( gradScalar, dPressuredStressJ2[i] ) ){
            results << "test_computeSecondOrderReferenceStressDecomposition (test 10) & False\n";
            return 1;
        }
    }

    //Test deriatives w.r.t. the right Cauchy-Green deformation tensor
    variableMatrix dDevStressdStressP, dDevStressdRCGP;
    variableVector dPressuredStressP, dPressuredRCGP;
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        error = micromorphicTools::computeSecondOrderReferenceStressDecomposition( S, C + delta, deviatoricResultJ, pressureResultJ );

        if ( error ){
            error->print();
            results << "test_computeSecondOrderReferenceStressDecomposition & False\n";
            return 1;
        }

        constantVector gradCol = ( deviatoricResultJ - deviatoricResult ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevStressdRCG[j][i], 1e-4 ) ){
                results << "test_computeSecondOrderReferenceStressDecomposition (test 11) & False\n";
                return 1;
            }
            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevStressdRCGJ2[j][i], 1e-4 ) ){
                results << "test_computeSecondOrderReferenceStressDecomposition (test 12) & False\n";
                return 1;
            }
        }

        constantType gradScalar = ( pressureResultJ - pressureResult ) / delta[i];

        if ( !vectorTools::fuzzyEquals( gradScalar, dPressuredRCG[i], 1e-4 ) ){
            results << "test_computeSecondOrderReferenceStressDecomposition (test 13) & False\n";
            return 1;
        }

        if ( !vectorTools::fuzzyEquals( gradScalar, dPressuredRCGJ2[i], 1e-4 ) ){
            results << "test_computeSecondOrderReferenceStressDecomposition (test 14) & False\n";
            return 1;
        }

        error = micromorphicTools::computeSecondOrderReferenceStressDecomposition( S, C + delta, deviatoricResultJ, pressureResultJ,
                                                                                   dDevStressdStressP, dDevStressdRCGP, 
                                                                                   dPressuredStressP, dPressuredRCGP );

        if ( error ){
            error->print();
            results << "test_computeSecondOrderReferenceStressDecomposition & False\n";
            return 1;
        }

        constantMatrix gradMat = ( dDevStressdStressP - dDevStressdStress ) / delta[i];

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        if ( !vectorTools::fuzzyEquals( gradMat[ 3 * j + k ][ 3 * l + m ],
                                                        d2DevStressdStressdRCGJ2[ 3 * j + k ][ 27 * l + 9 * m + 3 * ( int )( i / 3 ) + i % 3 ],
                                                        1e-4, 1e-5 ) ){
                            results << "test_computeSecondOrderReferenceStressDecomposition (test 15) & False\n";
                            return 1;
                        }
                    }
                }
            }
        }

        gradCol = ( dPressuredStressP - dPressuredStress ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], d2PressuredStressdRCGJ2[j][i], 1e-4 ) ){
                results << "test_computeSecondOrderReferenceStressDecomposition (test 16) & False\n";
                return 1;
            }
        }
    }

    results << "test_computeSecondOrderReferenceStressDecomposition & True\n";
    return 0;
}

int test_computeReferenceHigherOrderStressPressure( std::ofstream &results ){
    /*!
     * Test the computation of the higher order stress pressure.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector M = {  1,  2,  3,  4,  5,  6,  7,  8,  9,
                         10, 11, 12, 13, 14, 15, 16, 17, 18,
                         19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector C = { 0.34852835, 0.47540122, 1.11252634,
                         0.47540122, 1.49184663, 1.57435946,
                         1.11252634, 1.57435946, 3.68235756 };

    variableVector answer = { 69.06947837, 73.01858056, 76.96768276 };

    variableVector result;

    errorOut error = micromorphicTools::computeReferenceHigherOrderStressPressure( M, C, result );

    if ( error ){
        error->print();
        results << "test_computeReferenceHigherOrderStressPressure & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_computeReferenceHigherOrderStressPressure (test 1) & False\n";
        return 1;
    }

    //Test the Jacobians
    
    variableVector resultJ;
    variableMatrix dpdM, dpdC;
    error = micromorphicTools::computeReferenceHigherOrderStressPressure( M, C, resultJ, dpdM, dpdC );

    if ( error ){
        error->print();
        results << "test_computeReferenceHigherOrderStressPressure & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, resultJ ) ){
        results << "test_computeReferenceHigherOrderStressPressure (test 2) & False\n";
        return 1;
    }

    variableVector resultJ2;
    variableMatrix dpdMJ2, dpdCJ2, d2pdMdC;
    error = micromorphicTools::computeReferenceHigherOrderStressPressure( M, C, resultJ2, dpdMJ2, dpdCJ2, d2pdMdC );

    if ( error ){
        error->print();
        results << "test_computeReferenceHigherOrderStressPressure & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, resultJ2 ) ){
        results << "test_computeReferenceHigherOrderStressPressure (test 3) & False\n";
        return 1;
    }

    //Test dpdM
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < M.size(); i++ ){
        constantVector delta( M.size(), 0 );
        delta[i] = eps * fabs( M[i] ) + eps;

        error = micromorphicTools::computeReferenceHigherOrderStressPressure( M + delta, C, resultJ );

        if ( error ){
            error->print();
            results << "test_computeReferenceHigherOrderStressPressure & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dpdM[j][i] ) ){
                results << "test_computeReferenceHigherOrderStressPressure (test 4) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dpdMJ2[j][i] ) ){
                results << "test_computeReferenceHigherOrderStressPressure (test 5) & False\n";
                return 1;
            }
        }
    }

    //Test dpdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        error = micromorphicTools::computeReferenceHigherOrderStressPressure( M, C + delta, resultJ );

        if ( error ){
            error->print();
            results << "test_computeReferenceHigherOrderStressPressure & False\n";
            return 1;
        }

        constantVector gradCol = ( resultJ - result ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dpdC[j][i] ) ){
                results << "test_computeReferenceHigherOrderStressPressure (test 6) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dpdCJ2[j][i] ) ){
                results << "test_computeReferenceHigherOrderStressPressure (test 7) & False\n";
                return 1;
            }
        }
    }

    //Test d2pdMdC
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        error = micromorphicTools::computeReferenceHigherOrderStressPressure( M, C + delta, resultJ, dpdMJ2, dpdCJ2 );

        if ( error ){
            error->print();
            results << "test_computeReferenceHigherOrderStressPressure & False\n";
            return 1;
        }

        constantMatrix gradCol = ( dpdMJ2 - dpdM ) / delta[i];

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        if ( !vectorTools::fuzzyEquals( gradCol[ j ][ 9 * k + 3 * l + m ], 
                                                        d2pdMdC[ j ][ 81 * k + 27 * l + 9 * m + 3 * (int)( i / 3 ) + i % 3 ] ) ){
                            results << "test_computeReferenceHigherOrderStressPressure (test 7) & False\n";
                            return 1;
                        }
                    }
                }
            }
        }
    }

    results << "test_computeReferenceHigherOrderStressPressure & True\n";
    return 0;
}

int test_computeHigherOrderReferenceStressDecomposition( std::ofstream &results ){
    /*!
     * Test the computation of the decomposition of the higher order stress
     * into deviatoric and volumetric (pressure) parts.
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

    variableVector deviatoricAnswer = { 5.27633771,   1.68477127,  -3.21111472,  13.96704974,
                                        3.29864744, -10.4408607 ,  -4.31689678,  -1.90327486,
                                        3.27369895,  -2.4001685 ,   0.161758  ,   1.24944234,
                                       -6.01118653,  -2.07847616,   5.30384773,   2.46397506,
                                        0.3672135 ,  -1.56198261,  -8.15866529,  -0.72125951,
                                        6.32230906, -18.52878201,  -2.87440491,  14.2858097 ,
                                        4.98150383,   1.82382829,  -3.45626291 };

    variableVector pressureAnswer = { 0.56778037,  0.11342234, -0.43082224 };

    variableVector deviatoricResult, pressureResult;
    errorOut error = micromorphicTools::computeHigherOrderReferenceStressDecomposition( M, C, deviatoricResult, pressureResult );

    if ( error ){
        error->print();
        results << "test_computeHigherOrderReferenceStressDecomposition & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( deviatoricResult, deviatoricAnswer ) ){
        results << "test_computeHigherOrderReferenceStressDecomposition (test 1) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( pressureResult, pressureAnswer ) ){
        results << "test_computeHigherOrderReferenceStressDecomposition (test 2) & False\n";
        return 1;
    }

    //Test the Jacobians
    variableVector pressureResultJ;
    variableVector deviatoricResultJ;

    variableMatrix dPressuredStress, dPressuredRCG;
    variableMatrix dDevStressdStress, dDevStressdRCG;

    error = micromorphicTools::computeHigherOrderReferenceStressDecomposition( M, C, deviatoricResultJ, pressureResultJ,
                                                                               dDevStressdStress, dDevStressdRCG,
                                                                               dPressuredStress, dPressuredRCG );

    if ( error ){
        error->print();
        results << "test_computeHigherOrderReferenceStressDecomposition & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( pressureResultJ, pressureAnswer ) ){
        results << "test_computeHigherOrderReferenceStressDecomposition (test 3) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( deviatoricResultJ, deviatoricAnswer ) ){
        results << "test_computeHigherOrderReferenceStressDecomposition (test 4) & False\n";
        return 1;
    }


    variableVector pressureResultJ2;
    variableVector deviatoricResultJ2;

    variableMatrix dPressuredStressJ2, dPressuredRCGJ2;
    variableMatrix dDevStressdStressJ2, dDevStressdRCGJ2;

    variableMatrix d2DevStressdStressdRCGJ2, d2PressuredStressdRCGJ2;

    error = micromorphicTools::computeHigherOrderReferenceStressDecomposition( M, C, deviatoricResultJ2, pressureResultJ2,
                                                                               dDevStressdStressJ2, dDevStressdRCGJ2,
                                                                               dPressuredStressJ2, dPressuredRCGJ2,
                                                                               d2DevStressdStressdRCGJ2, d2PressuredStressdRCGJ2 );

    if ( error ){
        error->print();
        results << "test_computeHigherOrderReferenceStressDecomposition & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( pressureResultJ2, pressureAnswer ) ){
        results << "test_computeHigherOrderReferenceStressDecomposition (test 5) & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( deviatoricResultJ2, deviatoricAnswer ) ){
        results << "test_computeHigherOrderReferenceStressDecomposition (test 6) & False\n";
        return 1;
    }

    //Test deriatives w.r.t. the stress.
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < M.size(); i++ ){
        constantVector delta( M.size(), 0 );
        delta[i] = eps * fabs( M[i] ) + eps;

        error = micromorphicTools::computeHigherOrderReferenceStressDecomposition( M + delta, C, deviatoricResultJ, pressureResultJ );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderReferenceStressDecomposition & False\n";
            return 1;
        }

        constantVector gradCol = ( deviatoricResultJ - deviatoricResult ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevStressdStress[j][i] ) ){
                results << "test_computeHigherOrderReferenceStressDecomposition (test 7) & False\n";
                return 1;
            }
            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevStressdStressJ2[j][i] ) ){
                results << "test_computeHigherOrderReferenceStressDecomposition (test 8) & False\n";
                return 1;
            }
        }

        gradCol = ( pressureResultJ - pressureResult ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dPressuredStress[j][i] ) ){
                results << "test_computeHigherOrderReferenceStressDecomposition (test 9) & False\n";
                return 1;
            }
    
            if ( !vectorTools::fuzzyEquals( gradCol[j], dPressuredStressJ2[j][i] ) ){
                results << "test_computeHigherOrderReferenceStressDecomposition (test 10) & False\n";
                return 1;
            }
        }
    }

    //Test deriatives w.r.t. the right Cauchy-Green deformation tensor
    variableMatrix dDevStressdStressP, dDevStressdRCGP;
    variableMatrix dPressuredStressP, dPressuredRCGP;
    for ( unsigned int i = 0; i < C.size(); i++ ){
        constantVector delta( C.size(), 0 );
        delta[i] = eps * fabs( C[i] ) + eps;

        error = micromorphicTools::computeHigherOrderReferenceStressDecomposition( M, C + delta, deviatoricResultJ, pressureResultJ );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderReferenceStressDecomposition & False\n";
            return 1;
        }

        constantVector gradCol = ( deviatoricResultJ - deviatoricResult ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevStressdRCG[j][i], 1e-4 ) ){
                results << "test_computeHigherOrderReferenceStressDecomposition (test 11) & False\n";
                return 1;
            }
            if ( !vectorTools::fuzzyEquals( gradCol[j], dDevStressdRCGJ2[j][i], 1e-4 ) ){
                results << "test_computeHigherOrderReferenceStressDecomposition (test 12) & False\n";
                return 1;
            }
        }

        gradCol = ( pressureResultJ - pressureResult ) / delta[i];

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[j], dPressuredRCG[j][i], 1e-4 ) ){
                results << "test_computeHigherOrderReferenceStressDecomposition (test 13) & False\n";
                return 1;
            }
    
            if ( !vectorTools::fuzzyEquals( gradCol[j], dPressuredRCGJ2[j][i], 1e-4 ) ){
                results << "test_computeHigherOrderReferenceStressDecomposition (test 14) & False\n";
                return 1;
            }
        }

        error = micromorphicTools::computeHigherOrderReferenceStressDecomposition( M, C + delta, deviatoricResultJ, pressureResultJ,
                                                                                   dDevStressdStressP, dDevStressdRCGP, 
                                                                                   dPressuredStressP, dPressuredRCGP );

        if ( error ){
            error->print();
            results << "test_computeSecondOrderReferenceStressDecomposition & False\n";
            return 1;
        }

        constantMatrix gradMat = ( dDevStressdStressP - dDevStressdStress ) / delta[i];

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        for ( unsigned int n = 0; n < 3; n++ ){
                            for ( unsigned int o = 0; o < 3; o++ ){
                                if ( !vectorTools::fuzzyEquals( gradMat[ 9 * j + 3 * k + l ][ 9 * m + 3 * n + o ],
                                                                d2DevStressdStressdRCGJ2[ 9 * j + 3 * k + l ][ 81 * m + 27 * n + + 9 * o + 3 * ( int )( i / 3 ) + i % 3 ],
                                                                1e-4, 1e-5 ) ){
                                    results << "test_computeSecondOrderReferenceStressDecomposition (test 15) & False\n";
                                    return 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        gradMat = ( dPressuredStressP - dPressuredStress ) / delta[i];

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        if ( !vectorTools::fuzzyEquals( gradMat[ j ][ 9 * k + 3 * l + m ],
                                                        d2PressuredStressdRCGJ2[ j ][ 81 * k + 27 * l + 9 * m + 3 * ( int )( i / 3 ) + i % 3 ],
                                                        1e-4 ) ){
                            results << "test_computeSecondOrderReferenceStressDecomposition (test 16) & False\n";
                            return 1;
                        }
                    }
                }
            }
        }
    }
    
    results << "test_computeHigherOrderReferenceStressDecomposition & True\n";
    return 0;
}

int test_computeHigherOrderStressNorm( std::ofstream &results ){
    /*!
     * Test the computation of the special higher order stress norm.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector M = { 0.80732114,  0.79202055,  0.17990022,  0.97454675,  0.703207  ,
                        -0.58236697,  0.53324571, -0.93438873, -0.40650796,  0.14071918,
                         0.66933708, -0.67854069, -0.30317772, -0.93821882,  0.97270622,
                         0.00295302, -0.12441126,  0.30539971, -0.0580227 ,  0.89696105,
                         0.17567709, -0.9592962 ,  0.63535407,  0.95437804, -0.64531877,
                         0.69978907,  0.81327586 };

    variableVector answer = { 1.82692071, 2.24422424, 1.90780645 };

    variableVector result;

    errorOut error = micromorphicTools::computeHigherOrderStressNorm( M, result );

    if ( error ){
        results << "test_computeHigherOrderStressNorm & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_computeHigherOrderStressNorm (test 1) & False\n";
        return 1;
    }

    //Test the Jacobians
    variableVector resultJ;
    variableMatrix dNormMdMJ;

    error = micromorphicTools::computeHigherOrderStressNorm( M, resultJ, dNormMdMJ );

    if ( error ){
        results << "test_computeHigherOrderStressNorm & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, resultJ ) ){
        results << "test_computeHigherOrderStressNorm (test 2) & False\n";
        return 1;
    }

    variableVector resultJ2;
    variableMatrix dNormMdMJ2, d2NormMdM2J2;

    error = micromorphicTools::computeHigherOrderStressNorm( M, resultJ2, dNormMdMJ2, d2NormMdM2J2 );

    if ( error ){
        results << "test_computeHigherOrderStressNorm & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, resultJ2 ) ){
        results << "test_computeHigherOrderStressNorm (test 3) & False\n";
        return 1;
    }

    //Test dNormMdM
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < M.size(); i++ ){
        constantVector delta( M.size(), 0 );
        delta[i] = eps * fabs( M[i] ) + eps;

        variableVector resultP, resultM;

        error = micromorphicTools::computeHigherOrderStressNorm( M + delta, resultP );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderStressNorm & False\n";
            return 1;
        }

        error = micromorphicTools::computeHigherOrderStressNorm( M - delta, resultM );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderStressNorm & False\n";
            return 1;
        }

        variableVector gradCol = ( resultP - resultM ) / ( 2 * delta[i] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dNormMdMJ[ j ][ i ] ) ){
                results << "test_computeHigherOrderStressNorm (test 4) & False\n";
                return 1;
            }
        }

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dNormMdMJ2[ j ][ i ] ) ){
                results << "test_computeHigherOrderStressNorm (test 5) & False\n";
                return 1;
            }
        }

        variableMatrix derP, derM;

        error = micromorphicTools::computeHigherOrderStressNorm( M + delta, resultP, derP );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderStressNorm & False\n";
            return 1;
        }

        error = micromorphicTools::computeHigherOrderStressNorm( M - delta, resultM, derM );

        if ( error ){
            error->print();
            results << "test_computeHigherOrderStressNorm & False\n";
            return 1;
        }

        variableMatrix gradMat = ( derP - derM ) / ( 2 * delta[i] );

        unsigned int n = ( int )( i / 9 );
        unsigned int o = ( int )( ( i - 9 * n ) / 3 );
        unsigned int p = ( i - 9 * n - 3 * o ) % 3;

        for ( unsigned int j = 0; j < 3; j++ ){
            for ( unsigned int k = 0; k < 3; k++ ){
                for ( unsigned int l = 0; l < 3; l++ ){
                    for ( unsigned int m = 0; m < 3; m++ ){
                        if ( !vectorTools::fuzzyEquals( gradMat[ j ][ 9 * k + 3 * l + m ],
                                                        d2NormMdM2J2[ j ][ 243 * k + 81 * l + 27 * m + 9 * n + 3 * o + p ] ) ){
                            results << "test_computeHigherOrderStressNorm (test 6) & False\n";
                            return 1;
                        }
                    }
                }
            }
        }
    }

    results << "test_computeHigherOrderStressNorm & True\n";
    return 0;
}

int test_assembleDeformationGradient( std::ofstream &results ){
    /*!
     * Test the assembly of the deformation gradient from the gradient of the displacement.
     *
     * :param std::ofstream &results: The output file.
     */

    variableMatrix displacementGradient = { { 1, 2, 3 },
                                            { 4, 5, 6 },
                                            { 7, 8, 9 } };

    variableVector answer = { 2, 2, 3, 4, 6, 6, 7, 8, 10 };

    variableVector result;

    errorOut error = micromorphicTools::assembleDeformationGradient( displacementGradient, result );

    if ( error ){
        error->print();
        results << "test_assembleDeformationGradient & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_assembleDeformationGradient (test 1) & False\n";
        return 1;
    }

    //Tests of the Jacobians
    variableVector resultJ;
    variableMatrix dFdGradU;
    error = micromorphicTools::assembleDeformationGradient( displacementGradient, resultJ, dFdGradU );

    if ( error ){
        error->print();
        results << "test_assembleDeformationGradient & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, resultJ ) ){
        results << "test_assembleDeformationGradient (test 2) & False\n";
        return 1;
    }

    //Test the Jacobian w.r.t. the displacement gradient
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < 9; i++ ){
        constantMatrix delta = constantMatrix( 3, constantVector( 3, 0 ) );
        delta[ ( int )( i / 3) ][ i % 3 ] = eps * fabs( displacementGradient[ ( int )( i / 3 ) ][ i % 3 ] ) + eps;

        variableVector P, M;

        error = micromorphicTools::assembleDeformationGradient( displacementGradient + delta, P );

        if ( error ){
            error->print();
            results << "test_assembleDeformationGradient & False\n";
            return 1;
        }

        error = micromorphicTools::assembleDeformationGradient( displacementGradient - delta, M );

        if ( error ){
            error->print();
            results << "test_assembleDeformationGradient & False\n";
            return 1;
        }

        variableVector gradCol = ( P - M ) / ( 2 * delta[ ( int )( i / 3 ) ][ i % 3 ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dFdGradU[ j ][ i ] ) ){
                results << "test_assembledDeformationGradient (test 3) & False\n";
                return 1;
            }
        }
    }

    results << "test_assembleDeformationGradient & True\n";
    return 0;
}

int test_assembleMicroDeformation( std::ofstream &results ){
    /*!
     * Test the assembly of the micro deformation from the micro displacement
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector microDisplacement = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };

    variableVector answer = { 2, 2, 3, 4, 6, 6, 7, 8, 10 };

    variableVector result;

    errorOut error = micromorphicTools::assembleMicroDeformation( microDisplacement, result );

    if ( error ){
        error->print();
        results << "test_assembleMicroDeformation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_assembleMicroDeformation (test 1) & False\n";
        return 1;
    }

    //Test the Jacobians
    variableVector resultJ;
    variableMatrix dChidPhi;

    error = micromorphicTools::assembleMicroDeformation( microDisplacement, resultJ, dChidPhi );

    if ( error ){
        error->print();
        results << "test_assembleMicroDeformation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, resultJ ) ){
        results << "test_assembleMicroDeformation (test 2) & False\n";
        return 1;
    }

    //Test the jacobians w.r.t. phi
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < microDisplacement.size(); i++ ){
        constantVector delta( microDisplacement.size(), 0 );
        delta[i] = eps * fabs( microDisplacement[ i ] ) + eps;

        variableVector P, M;

        error = micromorphicTools::assembleMicroDeformation( microDisplacement + delta, P );

        if ( error ){
            error->print();
            results << "test_assembleMicroDeformation & False\n";
            return 1;
        }

        error = micromorphicTools::assembleMicroDeformation( microDisplacement - delta, M );

        if ( error ){
            error->print();
            results << "test_assembleMicroDeformation & False\n";
            return 1;
        }

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dChidPhi[ j ][ i ] ) ){
                results << "test_assembleMicroDeformation (test 3) & False\n";
                return 1;
            }
        }
    }

    results << "test_assembleMicroDeformation & True\n";
    return 0;
}

int test_assembleGradientMicroDeformation( std::ofstream &results ){
    /*!
     * Test the assembly of the micro deformation gradient
     *
     * :param std::ofstream &results: The output file.
     */

    variableMatrix gradientMicroDisplacement = { {  1,  2,  3 },
                                                 {  4,  5,  6 },
                                                 {  7,  8,  9 },
                                                 { 10, 11, 12 },
                                                 { 13, 14, 15 },
                                                 { 16, 17, 18 },
                                                 { 19, 20, 21 },
                                                 { 22, 23, 24 },
                                                 { 25, 26, 27 } };

    variableVector answer = { 1,  2,  3,  4,  5,  6,  7,  8,  9,
                             10, 11, 12, 13, 14, 15, 16, 17, 18,
                             19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector result;

    errorOut error = micromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, result );

    if ( error ){
        error->print();
        results << "test_assembleGradientMicroDeformation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_assembleGradientMicroDeformation (test 1) & False\n";
        return 1;
    }

    //Test the Jacobians
    variableVector resultJ;
    variableMatrix dGradChidGradPhi;

    error = micromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement, resultJ, dGradChidGradPhi );

    if ( error ){
        error->print();
        results << "test_assembleGradientMicroDeformation & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, resultJ ) ){
        results << "test_assembleGradientMicroDeformation (test 2) & False\n";
        return 1;
    }

    //Test the Jacobian w.r.t. the micro displacement gradient
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < 27; i++ ){
        constantMatrix delta( 9, constantVector( 3, 0 ) );
        delta[ ( int )( i / 3 ) ][ i % 3 ] = eps * fabs( gradientMicroDisplacement[ ( int )( i / 3 ) ][ i % 3 ] ) + eps;

        variableVector P, M;

        error = micromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement + delta, P );

        if ( error ){
            error->print();
            results << "test_assembleGradientMicroDeformation & False\n";
            return 1;
        }

        error = micromorphicTools::assembleGradientMicroDeformation( gradientMicroDisplacement - delta, M );

        if ( error ){
            error->print();
            results << "test_assembleGradientMicroDeformation & False\n";
            return 1;
        }

        variableVector gradCol = ( P - M ) / ( 2 * delta[ ( int )( i / 3 ) ][ i % 3 ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dGradChidGradPhi[ j ][ i ] ) ){
                results << "test_assembleGradientMicroDeformation (test 3) & False\n";
                return 1;
            }
        }
    }

    results << "test_assembleGradientMicroDeformation & True\n";
    return 0;
}

int test_pullBackCauchyStress( std::ofstream &results ){
    /*!
     * Test the pull back operation on the Cauchy stress
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector cauchyStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    variableVector deformationGradient  = { -1.08831037, -0.66333427, -0.48239487,
                                            -0.904554  ,  1.28942848, -0.02156112,
                                            -0.08464824, -0.07730218,  0.86415668 };

    variableVector answer = { -12.84536712,  -0.74608931,  13.40350133,
                                0.10367879,   0.04818336,   0.22121131,
                               20.04248815,   1.49351009, -18.33987341 };

    variableVector result, resultPF;

    errorOut error = micromorphicTools::pullBackCauchyStress( cauchyStress, deformationGradient, result );

    if ( error ){
        error->print();
        results << "test_pullBackCauchyStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_pullBackCauchyStress (test 1) & False\n";
        return 1;
    }

    error = micromorphicTools::pushForwardPK2Stress( result, deformationGradient, resultPF );

    if ( !vectorTools::fuzzyEquals( resultPF, cauchyStress ) ){
        results << "test_pullBackCauchyStress (test 2) & False\n";
        return 1;
    }

    //Test the Jacobians
    variableVector resultJ;
    variableMatrix dPK2dCauchy, dPK2dF;

    error = micromorphicTools::pullBackCauchyStress( cauchyStress, deformationGradient, resultJ, dPK2dCauchy, dPK2dF );

    if ( error ){
        error->print();
        results << "test_pullBackCauchyStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_pullBackCauchyStress (test 3) & False\n";
        return 1;
    }


    //Test the Jacobian w.r.t. the Cauchy stress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < cauchyStress.size(); i++ ){
        constantVector delta( cauchyStress.size(), 0 );
        delta[i] = eps * fabs( cauchyStress[ i ] ) + eps;

        variableVector P, M;

        error = micromorphicTools::pullBackCauchyStress( cauchyStress + delta, deformationGradient, P );

        if ( error ){
            error->print();
            results << "test_pullBackCauchyStress & False\n";
            return 1;
        }

        error = micromorphicTools::pullBackCauchyStress( cauchyStress - delta, deformationGradient, M );

        if ( error ){
            error->print();
            results << "test_pullBackCauchyStress & False\n";
            return 1;
        }

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPK2dCauchy[ j ][ i ] ) ){
                results << "test_pullBackCauchyStress (test 4) & False\n";
                return 1;
            }
        }
    }

    //Test the Jacobian w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[i] = eps * fabs( deformationGradient[ i ] ) + eps;

        variableVector P, M;

        error = micromorphicTools::pullBackCauchyStress( cauchyStress, deformationGradient + delta, P );

        if ( error ){
            error->print();
            results << "test_pullBackCauchyStress & False\n";
            return 1;
        }

        error = micromorphicTools::pullBackCauchyStress( cauchyStress, deformationGradient - delta, M );

        if ( error ){
            error->print();
            results << "test_pullBackCauchyStress & False\n";
            return 1;
        }

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dPK2dF[ j ][ i ] ) ){
                results << "test_pullBackCauchyStress (test 5) & False\n";
                return 1;
            }
        }
    }

    results << "test_pullBackCauchyStress & True\n";
    return 0;
}

int test_pullBackMicroStress( std::ofstream &results ){
    /*!
     * Test the pull back operation on the micro stress
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector microStress = { 1, 2, 3, 4, 5, 6, 7, 8, 9 };
    variableVector deformationGradient  = { -1.08831037, -0.66333427, -0.48239487,
                                            -0.904554  ,  1.28942848, -0.02156112,
                                            -0.08464824, -0.07730218,  0.86415668 };

    variableVector answer = { -12.84536712,  -0.74608931,  13.40350133,
                                0.10367879,   0.04818336,   0.22121131,
                               20.04248815,   1.49351009, -18.33987341 };

    variableVector result, resultPF;

    errorOut error = micromorphicTools::pullBackMicroStress( microStress, deformationGradient, result );

    if ( error ){
        error->print();
        results << "test_pullBackMicroStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( result, answer ) ){
        results << "test_pullBackMicroStress (test 1) & False\n";
        return 1;
    }

    error = micromorphicTools::pushForwardReferenceMicroStress( result, deformationGradient, resultPF );

    if ( !vectorTools::fuzzyEquals( resultPF, microStress ) ){
        results << "test_pullBackMicroStress (test 2) & False\n";
        return 1;
    }

    //Test the Jacobians
    variableVector resultJ;
    variableMatrix dSigmads, dSigmadF;

    error = micromorphicTools::pullBackMicroStress( microStress, deformationGradient, resultJ, dSigmads, dSigmadF );

    if ( error ){
        error->print();
        results << "test_pullBackMicroStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( resultJ, answer ) ){
        results << "test_pullBackMicroStress (test 3) & False\n";
        return 1;
    }


    //Test the Jacobian w.r.t. the micro stress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < microStress.size(); i++ ){
        constantVector delta( microStress.size(), 0 );
        delta[i] = eps * fabs( microStress[ i ] ) + eps;

        variableVector P, M;

        error = micromorphicTools::pullBackMicroStress( microStress + delta, deformationGradient, P );

        if ( error ){
            error->print();
            results << "test_pullBackMicroStress & False\n";
            return 1;
        }

        error = micromorphicTools::pullBackMicroStress( microStress - delta, deformationGradient, M );

        if ( error ){
            error->print();
            results << "test_pullBackMicroStress & False\n";
            return 1;
        }

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dSigmads[ j ][ i ] ) ){
                results << "test_pullBackMicroStress (test 4) & False\n";
                return 1;
            }
        }
    }

    //Test the Jacobian w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[i] = eps * fabs( deformationGradient[ i ] ) + eps;

        variableVector P, M;

        error = micromorphicTools::pullBackMicroStress( microStress, deformationGradient + delta, P );

        if ( error ){
            error->print();
            results << "test_pullBackMicroStress & False\n";
            return 1;
        }

        error = micromorphicTools::pullBackMicroStress( microStress, deformationGradient - delta, M );

        if ( error ){
            error->print();
            results << "test_pullBackMicroStress & False\n";
            return 1;
        }

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dSigmadF[ j ][ i ] ) ){
                results << "test_pullBackMicroStress (test 5) & False\n";
                return 1;
            }
        }
    }

    results << "test_pullBackMicroStress & True\n";
    return 0;
}

int test_pullBackHigherOrderStress( std::ofstream &results ){
    /*!
     * Test the push forward operation on the higher order stress tensor.
     *
     * :param std::ofstream &results: The output file.
     */

    variableVector higherOrderStress = {  1,  2,  3,  4,  5,  6,  7,  8,  9,
                                         10, 11, 12, 13, 14, 15, 16, 17, 18,
                                         19, 20, 21, 22, 23, 24, 25, 26, 27 };

    variableVector deformationGradient = { 0.29524861, -0.25221581, -1.60534711,
                                          -0.08817703,  0.28447808,  0.06451703,
                                           0.39849454,  0.38681512,  0.87840084 };

    variableVector microDeformation = { -0.25781969, -0.39826899, -0.79493259,
                                         0.38104724, -0.00830511, -0.51985409,
                                        -0.36415661, -0.6871168 ,  0.54018665 };

    variableVector answer = {   12.16637846,  -17.55206802,   -2.44609546,   46.49442604,
                               -66.93990198,   -9.1624221 ,  -11.22600619,   16.12630462,
                                 2.16292349,   36.56411317,  -52.34132435,   -6.7949869 ,
                               133.6526827 , -190.51226706,  -23.73373276,  -30.65330414,
                                43.46860157,    5.13634497,   -6.18698921,    8.71843323,
                                 0.96159717,  -20.55921994,   28.62926331,    2.72986588,
                                 4.14350171,   -5.66535471,   -0.40778148 };

    variableVector result;

    errorOut error = micromorphicTools::pullBackHigherOrderStress( higherOrderStress, deformationGradient,
                                                                   microDeformation, result );

    if ( error ){
        error->print();
        results << "test_pullBackHigherOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, result ) ){
        results << "test_pullBackHigherOrderStress (test 1) & False\n";
        return 1;
    }

    variableVector resultPF;

    error = micromorphicTools::pushForwardHigherOrderStress( result, deformationGradient, microDeformation, resultPF );

    if ( error ){
        error->print();
        results << "test_pullBackHigherOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( higherOrderStress, resultPF ) ){
        results << "test_pullBackHigherOrderStress (test 2) & False\n";
        return 1;
    }

    //Test the Jacobian
    variableVector resultJ;
    variableMatrix dMdm, dMdF, dMdChi;

    error = micromorphicTools::pullBackHigherOrderStress( higherOrderStress, deformationGradient,
                                                          microDeformation, resultJ, dMdm, dMdF, dMdChi );

    if ( error ){
        error->print();
        results << "test_pullBackHigherOrderStress & False\n";
        return 1;
    }

    if ( !vectorTools::fuzzyEquals( answer, resultJ ) ){
        results << "test_pullBackHigherOrderStress (test 3) & False\n";
        return 1;
    }

    //Test the jacobian w.r.t. the higher order stress
    constantType eps = 1e-6;
    for ( unsigned int i = 0; i < higherOrderStress.size(); i++ ){
        constantVector delta( higherOrderStress.size(), 0 );
        delta[ i ] = eps * fabs( higherOrderStress[ i ] ) + eps;

        variableVector P, M;

        error = micromorphicTools::pullBackHigherOrderStress( higherOrderStress + delta, deformationGradient,
                                                              microDeformation, P );
    
        if ( error ){
            error->print();
            results << "test_pullBackHigherOrderStress & False\n";
            return 1;
        }

        error = micromorphicTools::pullBackHigherOrderStress( higherOrderStress - delta, deformationGradient,
                                                              microDeformation, M );
    
        if ( error ){
            error->print();
            results << "test_pullBackHigherOrderStress & False\n";
            return 1;
        }

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMdm[ j ][ i ] ) ){
                results << "test_pullBackHigherOrderStress (test 4) & False\n";
                return 1;
            }
        }
    }

    //Test the jacobian w.r.t. the deformation gradient
    for ( unsigned int i = 0; i < deformationGradient.size(); i++ ){
        constantVector delta( deformationGradient.size(), 0 );
        delta[ i ] = eps * fabs( deformationGradient[ i ] ) + eps;

        variableVector P, M;

        error = micromorphicTools::pullBackHigherOrderStress( higherOrderStress, deformationGradient + delta,
                                                              microDeformation, P );
    
        if ( error ){
            error->print();
            results << "test_pullBackHigherOrderStress & False\n";
            return 1;
        }

        error = micromorphicTools::pullBackHigherOrderStress( higherOrderStress, deformationGradient - delta,
                                                              microDeformation, M );
    
        if ( error ){
            error->print();
            results << "test_pullBackHigherOrderStress & False\n";
            return 1;
        }

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMdF[ j ][ i ] ) ){
                results << "test_pullBackHigherOrderStress (test 5) & False\n";
                return 1;
            }
        }
    }

    //Test the jacobian w.r.t. the micro deformation
    for ( unsigned int i = 0; i < microDeformation.size(); i++ ){
        constantVector delta( microDeformation.size(), 0 );
        delta[ i ] = eps * fabs( microDeformation[ i ] ) + eps;

        variableVector P, M;

        error = micromorphicTools::pullBackHigherOrderStress( higherOrderStress, deformationGradient,
                                                              microDeformation + delta, P );
    
        if ( error ){
            error->print();
            results << "test_pullBackHigherOrderStress & False\n";
            return 1;
        }

        error = micromorphicTools::pullBackHigherOrderStress( higherOrderStress, deformationGradient,
                                                              microDeformation - delta, M );
    
        if ( error ){
            error->print();
            results << "test_pullBackHigherOrderStress & False\n";
            return 1;
        }

        variableVector gradCol = ( P - M ) / ( 2 * delta[ i ] );

        for ( unsigned int j = 0; j < gradCol.size(); j++ ){
            if ( !vectorTools::fuzzyEquals( gradCol[ j ], dMdChi[ j ][ i ] ) ){
                results << "test_pullBackHigherOrderStress (test 6) & False\n";
                return 1;
            }
        }
    }

    results << "test_pullBackHigherOrderStress & True\n";
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
    test_pullBackCauchyStress( results );
    test_pullBackMicroStress( results );
    test_pullBackHigherOrderStress( results );
    test_computeDeviatoricHigherOrderStress( results );
    test_computeDeviatoricReferenceHigherOrderStress( results );
    test_computeReferenceSecondOrderStressPressure( results );
    test_computeDeviatoricSecondOrderStress( results );
    test_computeDeviatoricReferenceSecondOrderStress( results );
    test_computeReferenceHigherOrderStressPressure( results );
    test_computeSecondOrderReferenceStressDecomposition( results );
    test_computeHigherOrderReferenceStressDecomposition( results );
    test_computeHigherOrderStressNorm( results );
    test_assembleDeformationGradient( results );
    test_assembleMicroDeformation( results );
    test_assembleGradientMicroDeformation( results );

    //Close the results file
    results.close();

    return 0;
}
