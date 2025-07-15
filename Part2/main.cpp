#include<iostream>
#include"DCT_basic_test.hpp"
#include<opencv2/core.hpp>
#include"DCT_FDCT.hpp"
#include <cmath>
#include <random>

void Basic_DCT_testing()
{
    DCTBasicTest test;
    test.test_vector(test.vector);
    test.test_matrix(test.matrix);
}

int main(int argc, char const *argv[])
{

    DCTFDCT dct_fdct;
    dct_fdct.test_matrix();

    return 0;
}
