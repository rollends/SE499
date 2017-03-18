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

    PolySpace polydiff_power(PolySpace poly, int order)
    {
        PolySpace prevPoly = poly;
        PolySpace dpoly = PolySpace::Constant(0);
        PolySpace powerCoeff = PolySpace::Constant(1);
        for(auto i = 0; i < order; i++)
        {
            Array<ValueType, 1, Dynamic> power(1, Spline::CoeffCount - i);
            power.setLinSpaced(0, Spline::PolyOrder - i);
            powerCoeff.tail(Spline::CoeffCount - i) = powerCoeff.tail(Spline::CoeffCount - i) * power;
        }
        dpoly.tail(Spline::CoeffCount - order) = poly.head(Spline::CoeffCount - order);
        dpoly = dpoly * powerCoeff;
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
     * The system is normally underdetermined so polyfit normally calculates the
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
        MatrixXT A(4*(segmentCount-1) + 4*(segmentCount), 2*CoeffCount*segmentCount);
        VectorXT b(A.rows(), 1);

        A.setZero();
        b.setZero();

        PolySpace powers = PolySpace::LinSpaced(0, Spline::PolyOrder);
        auto row = 0;
        for(auto indSeg = 0; indSeg < segmentCount; ++indSeg)
        {
            auto col = 2 * CoeffCount * indSeg;

            ValueType leftLambda = (indSeg * 1.0) / segmentCount;
            ValueType rightLambda = ((indSeg + 1.0) * 1.0) / segmentCount;

            PolySpace leftTerms = PolySpace::Constant(leftLambda);
            PolySpace rightTerms = PolySpace::Constant(rightLambda);
            leftTerms = pow(leftTerms, powers);
            rightTerms = pow(rightTerms, powers);

            PolySpace dLeftTerms = polydiff_power(leftTerms, 1);
            PolySpace ddLeftTerms = polydiff_power(leftTerms, 2);

            PolySpace dRightTerms = polydiff_power(rightTerms, 1);
            PolySpace ddRightTerms = polydiff_power(rightTerms, 2);

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

        // Moore-Penrose (Right) Pseudo-Inverse
        c = A.transpose() * (A * A.transpose()).fullPivLu().solve(b);

        // Map the resulting coefficient vector into our matrix form of the polynomial
        // spline!
        Matrix<ValueType, Dynamic, CoeffCount> result(2 * segmentCount, CoeffCount);
        result = Map< Matrix<ValueType, Dynamic, CoeffCount, RowMajor> >(c.data(), result.rows(), result.cols());

        return result;
    }

    /**
     * Splines a 5th order polynomial through the waypoints that ensures C2
     * conditions are met as well as a direction requirement at lambda=0
     *
     * The system is normally underdetermined so polyfit normally calculates the
     * right pseduo-inverse and performs a full pivot LU decomposition in order to
     * solve for the coefficients.
     */
    Matrix<ValueType, Dynamic, Spline::CoeffCount>
    polyfit(Matrix<ValueType, 2, Dynamic> waypoints, double direction)
    {
        constexpr auto CoeffCount = Spline::CoeffCount;

        auto const waypointCount = waypoints.cols();
        auto const segmentCount = waypointCount - 1;

        // Preallocate the necessary matrices
        constexpr size_t DirectionConstraintCount = 2;
        MatrixXT A(8 * (segmentCount-2) + 8 + 4 + DirectionConstraintCount, 2*CoeffCount*segmentCount);
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

            PolySpace dLeftTerms = polydiff_power(leftTerms, 1);
            PolySpace ddLeftTerms = polydiff_power(leftTerms, 2);

            PolySpace dRightTerms = polydiff_power(rightTerms, 1);
            PolySpace ddRightTerms = polydiff_power(rightTerms, 2);

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

            if( indSeg == 0 )
            {
                A.bottomRows(DirectionConstraintCount).block<1, CoeffCount>(0, col) = dLeftTerms;
                A.bottomRows(DirectionConstraintCount).block<1, CoeffCount>(1, col + CoeffCount) = dLeftTerms;

                Matrix<ValueType, 2, 1> vec;
                vec(0) = std::cos(direction);
                vec(1) = std::sin(direction);

                b.tail(DirectionConstraintCount) = vec;
            }

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

        // Moore-Penrose (Right) Pseudo-Inverse
        c = A.transpose() * (A * A.transpose()).fullPivLu().solve(b);

        // Map the resulting coefficient vector into our matrix form of the polynomial
        // spline!
        Matrix<ValueType, Dynamic, CoeffCount> result(2 * segmentCount, CoeffCount);
        result = Map< Matrix<ValueType, Dynamic, CoeffCount, RowMajor> >(c.data(), result.rows(), result.cols());

        return result;
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
    approximateSelf();
}

Spline::Spline(Eigen::Matrix<Spline::ValueType, 2, Eigen::Dynamic> points, double direction)
  : mSplineCount(points.cols() - 1),
    mPoly(2 * mSplineCount, CoeffCount),
    mDPoly(2 * mSplineCount, CoeffCount),
    mDDPoly(2 * mSplineCount, CoeffCount)
{
    mPoly = ::polyfit(points, direction);
    for( auto i = 0; i < mPoly.rows(); ++i )
    {
        mDPoly.row(i) = ::polydiffshift(mPoly.row(i));
    }
    for( auto i = 0; i < mPoly.rows(); ++i )
    {
        mDDPoly.row(i) = ::polydiffshift(mDPoly.row(i));
    }
    approximateSelf();
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

    double lambda = lambdaStar;
    double alpha = 1.0;
    alpha = std::min(alpha, alpha / speed(lambda));

    do
    {
        ValueType estimate = lambda - alpha * gradient(lambda);
        if( error(estimate) < error(lambda) )
        {
            lambda = estimate;
            alpha *= 1.2;
        }
        else
            alpha *= 0.8;
    } while( std::abs(alpha * gradient(lambda)) > 1e-12 );

    if( lambda < 0 )
    {
        // perturb initial condition.
        //return nearestPoint(point, lambdaStar + 1.0 / (2*mSplineCount));
    }

    return std::max(std::min(lambda, 1.0), 0.0);
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

void Spline::approximateSelf()
{
    Spline& self = *this;

    const double maxDelta = 1.0 / (2.0 * mSplineCount);
    double startLambda = 0;
    double endLambda = maxDelta;

    mApproximation.clear();
    do
    {
        Spline::Line currentSegment;
        Matrix<ValueType, 2, 1> start = self(startLambda);
        Matrix<ValueType, 2, 1> end = self(endLambda);

        double normedIntegral = 0.0;
        do
        {
            normedIntegral = 0.0;
            endLambda = startLambda / 2.0 + endLambda / 2.0;
            end = self(endLambda);

            Matrix<ValueType, 2, 1> midpoint = self(startLambda / 2.0 + endLambda / 2.0);

            // Integrate numerically
            for(int i = 0; (i < 1000) && (normedIntegral <= 0.1); ++i)
            {
                Matrix<ValueType, 3, 1> dError = Matrix<ValueType, 3, 1>::Constant(0);
                Matrix<ValueType, 3, 1> dLine = Matrix<ValueType, 3, 1>::Constant(0);
                dError.head(2) = self(startLambda + i * (endLambda - startLambda) / 1000)
                               - (start + i * (end - start) / 1000);
                dLine.head(2)  = (end - start) / 1000;
                normedIntegral += dError.cross(dLine).norm();
            }
        } while( normedIntegral > 0.1 );

        geom::set<0, 0>(currentSegment, start[0]);
        geom::set<0, 1>(currentSegment, start[1]);
        geom::set<1, 0>(currentSegment, end[0]);
        geom::set<1, 1>(currentSegment, end[1]);

        startLambda = endLambda;
        endLambda = std::min(endLambda + maxDelta, 1.0);

        mApproximation.push_back(currentSegment);
    } while( startLambda < 1.0 );
}

Matrix<Spline::ValueType, Dynamic, Spline::CoeffCount> Spline::poly() const { return mPoly; }
Matrix<Spline::ValueType, Dynamic, Spline::CoeffCount> Spline::dpoly() const{ return mDPoly; }
Matrix<Spline::ValueType, Dynamic, Spline::CoeffCount> Spline::ddpoly() const{ return mDDPoly; }
Spline::ApproximateSpline const & Spline::approximation() const { return mApproximation; }
