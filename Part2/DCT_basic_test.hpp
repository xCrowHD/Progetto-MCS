#ifndef DCT_BASIC_TEST
#define DCT_BASIC_TEST

#include<Eigen/Dense>
#include"DCT.hpp"
#include<iostream>

class DCTBasicTest
{

public:
    Eigen::VectorXd vector;
    Eigen::MatrixXd matrix;

public:

    DCTBasicTest(): vector(8), matrix(8, 8)
    {

        vector << 231, 32, 233, 161, 24, 71, 140, 245;

        matrix << 231,  32, 233, 161,  24,  71, 140, 245,
         247,  40, 248, 245, 124, 204,  36, 107,
         234, 202, 245, 167,   9, 217, 239, 173,
         193, 190, 100, 167,  43, 180,   8,  70,
          11,  24, 210, 177,  81, 243,   8, 112,
          97, 195, 203,  47, 125, 114, 165, 181,
         193,  70, 174, 167,  41,  30, 127, 245,
          87, 149,  57, 192,  65, 129, 178, 228;
    }

    ~DCTBasicTest() = default;

    void test_vector(const Eigen::VectorXd& v)
    {
        //std::cout.setf(std::ios::scientific);
        //std::cout.precision(2);
        DCT dct;
        Eigen::VectorXd a = dct.run_DCT(v);
        std::cout << "DCT1 del vettore:\n" << a.transpose() << std::endl; //Possibili differenze date da come viene calcolato il seno se ben minime
        //std::cout << "In tempo:\n" << dct.time << std::endl;
    }


    void test_matrix(const Eigen::MatrixXd& m)
    {
        //std::cout.setf(std::ios::scientific);
        //std::cout.precision(2);   
        DCT dct;
        Eigen::MatrixXd matrix_dct = dct.run_DCT2(m);
        std::cout << "DCT2 della matrice:\n" << matrix_dct << '\n'; //Possibili differenze date da come viene calcolato il seno se ben minime
        //std::cout << "In tempo:\n" << dct.time << std::endl;
    }

};

#endif