#include <boost/test/unit_test.hpp>

#include "mdot/include/babel_type.hpp"
#include "mdot/include/routine/eig_jacobi.hpp"

#include <iostream>

BOOST_AUTO_TEST_CASE(test_routine_mul_routine) {
  const size_t N = 3;
  const size_t M = 3;

// std::vector<std::vector<double>> A = {
//         {4, -2, 0, 0},
//         {-2, 4, -2, 0},
//         {0, -2, 3, -1},
//         {0, 0, -1, 2}
//     };

std::vector<std::vector<double>> A = {
        {4, 1},
        {1, 3}
    };

  // std::vector<dnum_t> matrix = {1, 2, 3, 4, 5, 6, 7, 8, 9};
  // std::vector<dnum_t> psi = {1, 1, 1};

  // size_t maxiter = 3;
  // const dnum_t tolerance = 1e-5;
  // dnum_t eigenvalue;
  // dnum_t eigenvector[3];



    DavidsonJacobi solver(100, 1e-6);
    std::vector<double> eigenvalues;
    std::vector<std::vector<double>> eigenvectors;

    solver.solve(A, 2, eigenvalues, eigenvectors);

    std::cout << "Eigenvalues:\n";
    for (double val : eigenvalues) {
        std::cout << val << std::endl;
    }

    std::cout << "Eigenvectors:\n";
    for (const auto& vec : eigenvectors) {
        for (double val : vec) {
            std::cout << val << " ";
        }
        std::cout << std::endl;
    }

  // mdot::lanczos_ev(matrix.data(), psi.data(), N, maxiter, tolerance, eigenvalue,
  //                  eigenvector);
  // BOOST_CHECK(eigenvalue == -4.0866302262564709);
}