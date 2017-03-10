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

    /** Default Constructor for Spline
     *
     * This should only be used for testing purposes. It produces a single
     * 5th order polynomial.
     */
    Spline();
    Spline(Eigen::Matrix<Spline::ValueType, 2, Eigen::Dynamic> points);
    Eigen::Matrix<ValueType, 2, 1> operator() (ValueType parameter, uint32_t derivative = 0) const;
    Eigen::Matrix<ValueType, 2, 2> frame(ValueType parameter, uint32_t derivative = 0) const;
    ValueType speed(ValueType parameter) const;

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
