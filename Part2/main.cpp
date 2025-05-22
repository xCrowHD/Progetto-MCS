#include<iostream>
#include<Eigen/Dense>
#include"DCT.hpp"

int main(int argc, char const *argv[])
{
    Eigen::VectorXd test_vector(8);
    test_vector << 231, 32, 233, 161, 24, 71, 140, 245;

    DCT dct;
    Eigen::VectorXd a = dct.run_DCT(test_vector);
    std::cout << a.transpose() << std::endl;


    return 0;
}
