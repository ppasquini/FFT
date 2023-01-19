#include <complex>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using cVector = Eigen::VectorXcd;
using cMatrix = Eigen::MatrixXcd;
using Complex = std::complex<double>;
using cSparseMatrix = Eigen::SparseMatrix<Complex>;
constexpr double pi = 3.14159265358979323846;