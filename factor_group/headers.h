#ifndef HEADERS_H
#define HEADERS_H
#include "../eigen-git-mirror/Eigen/Dense"
#include "../eigen-git-mirror/Eigen/Eigenvalues"
#include <vector>

#define PREC 1e-5

typedef Eigen::Matrix3d Matrix3d;
typedef Eigen::Vector3d Vector3d;
typedef Eigen::MatrixXi MatrixXi;

typedef Eigen::EigenSolver<Matrix3d> eigensolver;
Eigen::IOFormat CleanFmt(7, 0, ", ", "\n", "[", "]");

class sym_op_class
{
public:
    Matrix3d coord_matrix;
    Vector3d translation;
    std::string label;
};

class basis_class
{
public:
    std::string atom_type;
    Vector3d coordinates;
};

class crystal
{
public:
    Matrix3d lattice;
    std::vector<basis_class> basis;
};

#endif
