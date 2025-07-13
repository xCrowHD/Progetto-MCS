#include<iostream>
#include"DCT_basic_test.hpp"

int main(int argc, char const *argv[])
{
    std::cout.setf(std::ios::scientific);
    std::cout.precision(2);

    dct_basic_test test;
    test.test_vector(test.vector);
    test.test_matrix(test.matrix);
    
    
    return 0;
}
