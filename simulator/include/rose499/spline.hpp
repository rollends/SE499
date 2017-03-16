#ifndef ROSE499_SPLINE_HPP
#define ROSE499_SPLINE_HPP

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <cstdint>
#include <Eigen/Core>
#include "rose499/types.hpp"

namespace geom = boost::geometry;

/** 2D, 5th order spline data structure.
 *
 * Encapsulates evaluation of spline. This
 * class stores a 5th order spline.
 */
struct Spline
{
    using ValueType = SimulatorTypes::ValueType;
    using Point = geom::model::point<double, 2, geom::cs::cartesian>;
    using Line = geom::model::segment< Point >;
    using ApproximateSpline = std::vector< Line >;

    constexpr static int PolyOrder = 5;
    constexpr static int CoeffCount = PolyOrder + 1;

    /** Default Constructor for Spline
     *
     * This should only be used for testing purposes. It produces a single
     * 5th order polynomial.
     */
    Spline();
    Spline(Eigen::Matrix<Spline::ValueType, 2, Eigen::Dynamic> points);
    Eigen::Matrix<ValueType, 2, 1> operator() (ValueType parameter, uint32_t derivative = 0) const;
    Eigen::Matrix<ValueType, 2, 2> frame(ValueType parameter, uint32_t derivative = 0) const;
    ValueType nearestPoint(Eigen::Matrix<ValueType, 2, 1>, ValueType estimate) const;
    ValueType speed(ValueType parameter) const;
    int splineIndexUsed(ValueType parameter) const;
    ApproximateSpline const & approximation() const;

    Eigen::Matrix<ValueType, Eigen::Dynamic, CoeffCount> poly() const;

    class SplineException{ };
    class InvalidParameterException : public SplineException { };
    class IncorrectSplineFormatException : public SplineException { };

private:
    void approximateSelf();

    int mSplineCount;
    Eigen::Matrix<ValueType, Eigen::Dynamic, CoeffCount> mPoly;
    Eigen::Matrix<ValueType, Eigen::Dynamic, CoeffCount> mDPoly;
    Eigen::Matrix<ValueType, Eigen::Dynamic, CoeffCount> mDDPoly;
    ApproximateSpline mApproximation;
};

#endif
