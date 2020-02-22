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

    results << "test_pushForwardReferenceMicroStrain & True\n";
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
    test_pushForwardReferenceMicroStress( results );

    //Close the results file
    results.close();

    return 0;
}
