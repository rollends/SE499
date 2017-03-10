#include <cmath>
#include <Eigen/LU>
#include "rose499/spline.hpp"

using namespace Eigen;

namespace
{
    constexpr int ResampleCount = 9;

    using ValueType = Spline::ValueType;
    using MatrixXT = SimulatorTypes::MatrixXT;
    using VectorXT = SimulatorTypes::VectorXT;

    Matrix<ValueType, 2, Dynamic>
    generateintermediates(Matrix<ValueType, 2, Dynamic> waypoints )
    {
        auto segmentCount = waypoints.cols() - 1;
        Matrix<ValueType, 2, Dynamic> set(2, segmentCount * ResampleCount);

        for(int i = 0; i < segmentCount; ++i)
        {
            set.block(0, ResampleCount * i, 1, ResampleCount) =
                VectorXd::LinSpaced(11, waypoints(0, i), waypoints(0, i + 1)).block<1, ResampleCount>(0, 1);
            set.block(1, ResampleCount * i, 1, ResampleCount) =
                VectorXd::LinSpaced(11, waypoints(1, i), waypoints(1, i + 1)).block<1, ResampleCount>(0, 1);
        }
    }

    //
    Matrix<ValueType, Dynamic, Spline::PolyOrder + 1>
    polyfit(Matrix<ValueType, 2, Dynamic> waypoints)
    {
        using PolySpace = Array<Spline::ValueType, Spline::PolyOrder + 1, 1>;

        auto samplePoints = generateintermediates(waypoints);
        auto const sampleCount = samplePoints.cols();
        auto const waypointCount = waypoints.cols();
        auto const segmentCount = waypointCount - 1;
        auto const pointCount = sampleCount + waypointCount;

        // Preallocate the necessary matrices
        MatrixXT A(2 * sampleCount, 2 * (Spline::PolyOrder + 1) * segmentCount);
        MatrixXT C(2 * waypointCount + 2 * 2 * (waypointCount - 2), 2 * (Spline::PolyOrder + 1) * segmentCount);
        MatrixXT S(A.cols() + C.rows(), A.cols() + C.rows());
        VectorXT b(A.rows(), 1);
        VectorXT d(C.rows(), 1);
        VectorXT y(A.cols() + d.rows(), 1);

        S.fill(0);
        A.fill(0);
        C.fill(0);

        // Generate least squares matrix.
        PolySpace powers = PolySpace::LinSpaced(0, Spline::PolyOrder);
        int segmentInd = 0;
        for(int i = 0; i < pointCount; ++i)
        {
            // Skip waypoints
            if( i % (ResampleCount+1) == 0 )
                continue;

            ValueType lambda = i * (1.0 / (pointCount - 1));
            PolySpace lambdaPow = PolySpace::Constant(lambda);
            A.block<1, Spline::PolyOrder + 1>(2*i, (0 + segmentInd * 2) * (Spline::PolyOrder+1) ) = pow(lambdaPow, powers);
            A.block<1, Spline::PolyOrder + 1>(2*i, (1 + segmentInd * 2) * (Spline::PolyOrder+1) ) = pow(lambdaPow, powers);
        }

        // Generate equality constraint matrix.

        // Apply Least Squares with Equality constraints procedure described in
        // https://stanford.edu/class/ee103/lectures/constrained-least-squares/constrained-least-squares_slides.pdf
        // Involves solving linear equation of the form,
        //          S*x = y
        // Where x concatenates the polynomial coefficients and Lagrange multipliers.
        // Where y concatenates the transformed constraint values. (A'*b ; d)
        // For equations,
        //          min( |A*c - b| ) (least squares)
        //          C*c = d (equality)
        S.block(0, 0, A.cols(), A.cols()) = 2 * A.transpose() * A;
        S.block(0, A.cols(), C.cols(), C.rows()) = C.transpose();
        S.block(A.cols(), 0, C.rows(), C.cols()) = C;
        y.block(0, 0, A.cols(), 1) = 2 * A.transpose() * b;
        y.block(A.cols(), 0, d.rows(), 1) = d;

        auto c = S.fullPivLu().solve(y).head(A.cols());

        Matrix<ValueType, Dynamic, Spline::PolyOrder + 1> result(2, (Spline::PolyOrder + 1));

        // Enforce going through way points, construct
        return result;
    }
}

Spline::Spline(Matrix<Spline::ValueType, 2, Dynamic> points)
  : mSplineCount(points.cols() - 1),
    mPoly(2 * mSplineCount, PolyOrder + 1),
    mDPoly(2 * mSplineCount, PolyOrder + 1),
    mDDPoly(2 * mSplineCount, PolyOrder + 1)
{
    mPoly = ::polyfit(points);
}

Spline::Spline()
  : mSplineCount(1),
    mPoly(2 * mSplineCount, PolyOrder + 1),
    mDPoly(2 * mSplineCount, PolyOrder + 1),
    mDDPoly(2 * mSplineCount, PolyOrder + 1)
{
    mPoly.fill(0);
    mDPoly.fill(0);
    mDDPoly.fill(0);

    mPoly.block<2, 1>(0, 0);
}

Matrix<Spline::ValueType, 2, 2> Spline::frame(ValueType parameter, uint32_t derivative) const
{
    Matrix<ValueType, 2, 2> basis;
    Matrix<ValueType, 2, 1> tangent = (*this)(parameter, 1);
    ValueType speed = tangent.norm();

    if(derivative >= 0)
    {
        basis.col(0) = tangent / tangent.norm();
    }

    if(derivative >= 1)
    {
        Matrix<ValueType, 2, 1> accel = (*this)(parameter, 2);
        basis.col(0) = accel / speed - accel.dot(tangent) / (speed * speed) * basis.col(0);
    }

    if(derivative >= 2)
        throw InvalidParameterException();

    // Rotate the first column by 90 degrees (right handed frame) to form the
    // full basis for R^2
    basis.col(1) = basis.col(0).reverse();
    basis(0, 1) = -basis(0, 1);
    return basis;
}

Spline::ValueType Spline::speed(ValueType parameter) const
{
    return (*this)(parameter, 1).norm();
}

Matrix<Spline::ValueType, 2, 1> Spline::operator() (ValueType parameter, uint32_t derivative) const
{
    using Vector = Matrix<Spline::ValueType, 2, 1>;
    using PolySpace = Array<Spline::ValueType, PolyOrder+1, 1>;

    // Parameter out of bounds
    if( parameter < 0 || parameter > 1 )
        throw InvalidParameterException();

    // Identically zero derivative
    if( derivative > PolyOrder )
        return Vector::Zero();

    // Calculate x^7, x^6, ..., x^0
    PolySpace value = PolySpace::Constant(parameter);
    PolySpace powers = PolySpace::LinSpaced(0, PolyOrder);
    value = pow(value, powers).reverse();

    // Choose spline based on parameter value.
    int indSpline = std::min(mSplineCount, (int)std::floor(mSplineCount * parameter));

    // And evaluate!
    switch(derivative)
    {
    case 0:
        return mPoly.block<2, PolyOrder + 1>(2 * indSpline, 0) * value.matrix();

    case 1:
        return mDPoly.block<2, PolyOrder + 1>(2 * indSpline, 0) * value.matrix();

    case 2:
        return mDDPoly.block<2, PolyOrder + 1>(2 * indSpline, 0) * value.matrix();

    default:
        throw InvalidParameterException();
    };
}
