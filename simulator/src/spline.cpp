#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "rose499/spline.hpp"

using namespace Eigen;

constexpr int Spline::PolyOrder;
constexpr int Spline::CoeffCount;

namespace
{
    constexpr int ResampleCount = 9;

    using ValueType = Spline::ValueType;
    using MatrixXT = SimulatorTypes::MatrixXT;
    using VectorXT = SimulatorTypes::VectorXT;
    using PolySpace = Array<ValueType, 1, Spline::CoeffCount>;

    PolySpace polydiff(PolySpace poly)
    {
        PolySpace dpoly = PolySpace::Constant(0);
        dpoly.tail(Spline::CoeffCount - 1) = poly.head(Spline::CoeffCount - 1);
        dpoly *= PolySpace::LinSpaced(0, Spline::PolyOrder);
        return dpoly;
    }

    PolySpace polydiffshift(PolySpace poly)
    {
        PolySpace dpoly = PolySpace::Constant(0);
        dpoly.head(Spline::CoeffCount - 1) = poly.tail(Spline::CoeffCount - 1);
        dpoly.head(Spline::CoeffCount - 1) *= Array<ValueType, Spline::CoeffCount - 1, 1>::LinSpaced(1, Spline::PolyOrder);
        return dpoly;
    }

    /**
     * Splines a 5th order polynomial through the waypoints that ensures C2
     * conditions are met.
     *
     * If the system is normally underdetermined so polyfit normally calculates the
     * right pseduo-inverse and performs a full pivot LU decomposition in order to
     * solve for the coefficients.
     */
    Matrix<ValueType, Dynamic, Spline::CoeffCount>
    polyfit(Matrix<ValueType, 2, Dynamic> waypoints)
    {
        constexpr auto CoeffCount = Spline::CoeffCount;

        auto const waypointCount = waypoints.cols();
        auto const segmentCount = waypointCount - 1;

        // Preallocate the necessary matrices
        MatrixXT A(8 * (segmentCount-2) + 8 + 4, 2*CoeffCount*segmentCount);
        VectorXT b(A.rows(), 1);

        A.setZero();
        b.setZero();

        PolySpace powers = PolySpace::LinSpaced(0, Spline::PolyOrder);
        auto row = 0;
        for(auto indSeg = 0; indSeg < segmentCount; ++indSeg)
        {
            auto col = 2*CoeffCount*indSeg;

            ValueType leftLambda = (indSeg * 1.0) / segmentCount;
            ValueType rightLambda = ((indSeg + 1.0) * 1.0) / segmentCount;

            PolySpace leftTerms = PolySpace::Constant(leftLambda);
            PolySpace rightTerms = PolySpace::Constant(rightLambda);
            leftTerms = pow(leftTerms, powers);
            rightTerms = pow(rightTerms, powers);

            PolySpace dLeftTerms = polydiff(leftTerms);
            PolySpace ddLeftTerms = polydiff(dLeftTerms);

            PolySpace dRightTerms = polydiff(rightTerms);
            PolySpace ddRightTerms = polydiff(dRightTerms);

            // Left Equality (C0)
            A.block<1, CoeffCount>(row + 0, col) = leftTerms;
            A.block<1, CoeffCount>(row + 1, col + CoeffCount) = leftTerms;
            b.block<2, 1>(row, 0) = waypoints.block<2, 1>(0, indSeg);
            row += 2;

            // Right Equality (C0)
            A.block<1, CoeffCount>(row + 0, col) = rightTerms;
            A.block<1, CoeffCount>(row + 1, col + CoeffCount) = rightTerms;
            b.block<2, 1>(row, 0) = waypoints.block<2, 1>(0, indSeg + 1);
            row += 2;

            if( indSeg < segmentCount - 1 )
            {
                // Right Differential Continuity (C1)
                A.block<1, CoeffCount>(row + 0, col) = dRightTerms;
                A.block<1, CoeffCount>(row + 0, col + 2 * CoeffCount) = -dRightTerms;
                A.block<1, CoeffCount>(row + 1, col + CoeffCount) = dRightTerms;
                A.block<1, CoeffCount>(row + 1, col + 3 * CoeffCount) = -dRightTerms;
                row += 2;

                // Right Differential Continuity (C2)
                A.block<1, CoeffCount>(row + 0, col) = ddRightTerms;
                A.block<1, CoeffCount>(row + 0, col + 2 * CoeffCount) = -ddRightTerms;
                A.block<1, CoeffCount>(row + 1, col + CoeffCount) = ddRightTerms;
                A.block<1, CoeffCount>(row + 1, col + 3 * CoeffCount) = -ddRightTerms;
                row += 2;
            }
        }

        MatrixXT c(A.cols(), 1);

        if( A.cols() == A.rows() )          // Exactly determined system (hopefully lol)
            c = A.fullPivLu().solve(b);
        else if( A.cols() < A.rows() )      // Use Left Inverse
            c = (A.transpose() * A).fullPivLu().solve(A.transpose() * b);
        else                                // Use Right Inverse
            c = A.transpose() * (A * A.transpose()).fullPivLu().solve(b);

        // Map the resulting coefficient vector into our matrix form of the polynomial
        // spline!
        Matrix<ValueType, Dynamic, CoeffCount> result(2 * segmentCount, CoeffCount);
        result = Map< Matrix<ValueType, Dynamic, CoeffCount, RowMajor> >(c.data(), result.rows(), result.cols());

        return result;
    }


    Spline::ApproximateSpline approximate(Matrix<ValueType, 2, Dynamic> poly)
    {
        Spline::ApproximateSpline spline;
        return spline;
    }
}

Spline::Spline(Matrix<Spline::ValueType, 2, Dynamic> points)
  : mSplineCount(points.cols() - 1),
    mPoly(2 * mSplineCount, CoeffCount),
    mDPoly(2 * mSplineCount, CoeffCount),
    mDDPoly(2 * mSplineCount, CoeffCount)
{
    mPoly = ::polyfit(points);
    for( auto i = 0; i < mPoly.rows(); ++i )
    {
        mDPoly.row(i) = ::polydiffshift(mPoly.row(i));
    }
    for( auto i = 0; i < mPoly.rows(); ++i )
    {
        mDDPoly.row(i) = ::polydiffshift(mDPoly.row(i));
    }
    //mApproximation = ::approximate(mPoly);
}

Spline::Spline()
  : mSplineCount(1),
    mPoly(2 * mSplineCount, CoeffCount),
    mDPoly(2 * mSplineCount, CoeffCount),
    mDDPoly(2 * mSplineCount, CoeffCount)
{
    mPoly.setZero();
    mDPoly.setZero();
    mDDPoly.setZero();
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
        basis.col(0) = ( accel - accel.dot(basis.col(0)) * basis.col(0) ) / speed;

        // If we have a 0 curvature, generate the trivial basis (for {0})
        if( !basis.col(0).isZero() )
            basis.col(0) /= basis.col(0).norm();
    }

    if(derivative >= 2)
        throw InvalidParameterException();

    // Rotate the first column by 90 degrees (right handed frame) to form the
    // full basis for R^2
    basis.col(1) = basis.col(0).reverse();
    basis(0, 1) = -basis(0, 1);
    return basis;
}

Spline::ValueType Spline::nearestPoint(Matrix<ValueType, 2, 1> point, ValueType lambdaStar) const
{
    auto error =
        [point, this](ValueType l) -> ValueType {
            return (point - this->operator()(l)).squaredNorm();
        };
    auto gradient =
        [point, this](ValueType l) -> ValueType {
            Matrix<ValueType, 2, 1> verr = point - this->operator()(l);
            return -2*verr.dot(this->operator()(l, 1));
        };

    double alpha = 1.0 / (10 * mSplineCount);
    alpha = std::min(alpha, alpha / speed(lambdaStar));

    do
    {
        ValueType estimate = lambdaStar - alpha * gradient(lambdaStar);
        if( error(estimate) <= error(lambdaStar) )
        {
            lambdaStar = estimate;
            alpha *= 1.2;
        }
        else
            alpha *= 0.8;
    } while( std::abs(alpha * gradient(lambdaStar)) > 1e-12 );

    return lambdaStar;
}

int Spline::splineIndexUsed(ValueType parameter) const
{
    return std::max(0, std::min(mSplineCount - 1, (int)std::floor(mSplineCount * parameter)));
}

Spline::ValueType Spline::speed(ValueType parameter) const
{
    return (*this)(parameter, 1).norm();
}

Matrix<Spline::ValueType, 2, 1> Spline::operator() (ValueType parameter, uint32_t derivative) const
{
    using Vector = Matrix<Spline::ValueType, 2, 1>;
    using PolySpace = Array<Spline::ValueType, CoeffCount, 1>;

    // Identically zero derivative
    if( derivative > PolyOrder )
        return Vector::Zero();

    // Calculate x^0, x^1, ..., x^5
    PolySpace value = PolySpace::Constant(parameter);
    PolySpace powers = PolySpace::LinSpaced(0, PolyOrder);
    value = pow(value, powers);

    // Choose spline based on parameter value.
    int indSpline = splineIndexUsed(parameter);

    // And evaluate!
    switch(derivative)
    {
    case 0:
        return mPoly.block<2, CoeffCount>(2 * indSpline, 0) * value.matrix();

    case 1:
        return mDPoly.block<2, CoeffCount>(2 * indSpline, 0) * value.matrix();

    case 2:
        return mDDPoly.block<2, CoeffCount>(2 * indSpline, 0) * value.matrix();

    default:
        throw InvalidParameterException();
    };
}

Matrix<Spline::ValueType, Dynamic, Spline::CoeffCount> Spline::poly() const { return mPoly; }
Spline::ApproximateSpline const & Spline::approximation() const { return mApproximation; }
