#ifndef ROSE499_TYPES_HPP
#define ROSE499_TYPES_HPP

#include <Eigen/Core>
#include <iostream>

namespace SimulatorTypes
{
    typedef long double ValueType;
    typedef Eigen::Matrix<ValueType, Eigen::Dynamic, Eigen::Dynamic> MatrixXT;
    typedef Eigen::Matrix<ValueType, Eigen::Dynamic, 1> VectorXT;
    typedef Eigen::Matrix<ValueType, 2, 1> Vector2T;
};

#endif
