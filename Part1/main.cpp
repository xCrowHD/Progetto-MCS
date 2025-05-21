#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "linear_resolver.hpp"

void _3x3Test()
{
    // Matrice DENSE
    Eigen::MatrixXd A_dense(3,3);
    A_dense << 10, -1, 2,
               -1, 11, -1,
                2, -1, 10;

    double tol = 1e-6;

    std::cout << "\n=== Test con MATRICE DENSE ===\n";
    linear_resolver<Eigen::MatrixXd> solver_dense(A_dense, tol);
    solver_dense.run_resolvers();

    // Matrice SPARSE equivalente
    Eigen::SparseMatrix<double> A_sparse(3, 3);
    A_sparse.insert(0, 0) = 10;
    A_sparse.insert(0, 1) = -1;
    A_sparse.insert(0, 2) = 2;

    A_sparse.insert(1, 0) = -1;
    A_sparse.insert(1, 1) = 11;
    A_sparse.insert(1, 2) = -1;

    A_sparse.insert(2, 0) = 2;
    A_sparse.insert(2, 1) = -1;
    A_sparse.insert(2, 2) = 10;

    A_sparse.makeCompressed();

    std::cout << "\n=== Test con MATRICE SPARSA ===\n";
    linear_resolver<Eigen::SparseMatrix<double>> solver_sparse(A_sparse, tol);
    solver_sparse.run_resolvers();
}

int main() {
    
    //_3x3Test();

    Eigen::SparseMatrix<double> spa1_sparse;

    if (!Eigen::loadMarket(spa1_sparse, "./dati/spa1.mtx")) 
    {
        std::cerr << "Errore nel caricamento della matrice sparsa." << std::endl;
    }

    Eigen::MatrixXd spa1_dense = Eigen::MatrixXd(spa1_sparse);

    double tol = 1e-6;

    std::cout << "\n=== Test con MATRICE DENSE ===\n";
    linear_resolver<Eigen::MatrixXd> dense_solver(spa1_dense, tol);
    dense_solver.run_resolvers();

    std::cout << "\n=== Test con MATRICE SPARSE ===\n";
    linear_resolver<Eigen::SparseMatrix<double>> sparse_solver(spa1_sparse, tol);
    sparse_solver.run_resolvers();   

    return 0;
}

