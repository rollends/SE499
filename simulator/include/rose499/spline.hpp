#ifndef ROSE499_SPLINE_HPP
#define ROSE499_SPLINE_HPP

#include <cstdint>
#include <Eigen/Core>
#include "rose499/types.hpp"

/** 2D, 5th order spline data structure.
 *
 * Encapsulates evaluation of spline. This
 * class stores a 5th order spline.
 */
struct Spline
{
    typedef SimulatorTypes::ValueType ValueType;
    constexpr static int PolyOrder = 5;

    Spline(Eigen::MatrixXd points);
    Eigen::Array<ValueType, 2, 1> operator() (ValueType parameter, uint32_t derivative = 0) const;
    Eigen::Matrix<ValueType, 2, 2> frame(ValueType parameter) const;

    class SplineException{ };
    class InvalidParameterException : public SplineException { };
    class IncorrectSplineFormatException : public SplineException { };

private:
    const int mSplineCount;
    Eigen::Matrix<ValueType, Eigen::Dynamic, PolyOrder + 1> mPoly;
    Eigen::Matrix<ValueType, Eigen::Dynamic, PolyOrder + 1> mDPoly;
    Eigen::Matrix<ValueType, Eigen::Dynamic, PolyOrder + 1> mDDPoly;
};

#endif
