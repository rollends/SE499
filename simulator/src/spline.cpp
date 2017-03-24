#include <cmath>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <limits>
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

/*
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
*/
    /**
        Builds an equality matrix that only provides end point guarantees and
        C^2 continuity everywhere.
    **/
    void
    polyBuildSoftEqualityMatrix(    MatrixXT& Aeq,
                                    VectorXT& beq,
                                    Matrix<ValueType, 2, Dynamic> const & waypoints,
                                    double direction = std::numeric_limits<double>::infinity() )
    {
        constexpr size_t DirectionConstraintCount = 1;
        constexpr auto CoeffCount = Spline::CoeffCount;

        auto const waypointCount = waypoints.cols();
        auto const segmentCount = waypointCount - 1;

        bool const supportDirection = std::isfinite(direction);

        Aeq.resize( 4*(waypointCount-2)                                         // Differential Constraints
                        + 2*waypointCount                                       // Equality Constraints
                        + (supportDirection ? DirectionConstraintCount : 0)     // Direction Constraint
                    , 2*CoeffCount*segmentCount );
        beq.resize(Aeq.rows(), 1);

        Aeq.setZero();
        beq.setZero();

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


            if( indSeg < segmentCount - 1 )
            {
                // Right Equality (C0)
                Aeq.block<1, CoeffCount>(row + 0, col) = rightTerms;
                Aeq.block<1, CoeffCount>(row + 0, col + 2 * CoeffCount) = -rightTerms;
                Aeq.block<1, CoeffCount>(row + 1, col + CoeffCount) = rightTerms;
                Aeq.block<1, CoeffCount>(row + 1, col + 3 * CoeffCount) = -rightTerms;
                row += 2;

                // Right Differential Continuity (C1)
                Aeq.block<1, CoeffCount>(row + 0, col) = dRightTerms;
                Aeq.block<1, CoeffCount>(row + 0, col + 2 * CoeffCount) = -dRightTerms;
                Aeq.block<1, CoeffCount>(row + 1, col + CoeffCount) = dRightTerms;
                Aeq.block<1, CoeffCount>(row + 1, col + 3 * CoeffCount) = -dRightTerms;
                row += 2;

                // Right Differential Continuity (C2)
                Aeq.block<1, CoeffCount>(row + 0, col) = ddRightTerms;
                Aeq.block<1, CoeffCount>(row + 0, col + 2 * CoeffCount) = -ddRightTerms;
                Aeq.block<1, CoeffCount>(row + 1, col + CoeffCount) = ddRightTerms;
                Aeq.block<1, CoeffCount>(row + 1, col + 3 * CoeffCount) = -ddRightTerms;
                row += 2;
            }
            else
            {
                // Hard Right Equality
                Aeq.block<1, CoeffCount>(row + 0, col) = rightTerms;
                Aeq.block<1, CoeffCount>(row + 1, col + CoeffCount) = rightTerms;
                beq.block<2, 1>(row, 0) = waypoints.block<2, 1>(0, indSeg + 1);
                row += 2;
            }

            if( indSeg == 0 )
            {
                // Left Equality (C0)
                Aeq.block<1, CoeffCount>(row + 0, col) = leftTerms;
                Aeq.block<1, CoeffCount>(row + 1, col + CoeffCount) = leftTerms;
                beq.block<2, 1>(row, 0) = waypoints.block<2, 1>(0, indSeg);
                row += 2;

                if( supportDirection )
                {
                    while(direction > M_PI) direction -= 2*M_PI;
                    while(direction < -M_PI) direction += 2*M_PI;

                    if( std::abs(direction) <= std::atan2(1.0, 1.0)
                     || std::abs(direction - M_PI) <= std::atan2(1.0, 1.0)
                     || std::abs(direction + M_PI) <= std::atan2(1.0, 1.0) )
                    {
                        // Small tangent.
                        Aeq.bottomRows(DirectionConstraintCount).block<1, CoeffCount>(0, col) = std::tan(direction)*dLeftTerms;
                        Aeq.bottomRows(DirectionConstraintCount).block<1, CoeffCount>(0, col + CoeffCount) = -dLeftTerms;
                    }
                    else
                    {
                        // Large tangent
                        Aeq.bottomRows(DirectionConstraintCount).block<1, CoeffCount>(0, col) = -dLeftTerms;
                        Aeq.bottomRows(DirectionConstraintCount).block<1, CoeffCount>(0, col + CoeffCount) = dLeftTerms / std::tan(direction);
                    }
                }
            }
        }
    }

    /**
        Builds the equality matrix that is used in least squares to achieve
        C0, C1, and C2 equality at waypoints.

    **/
    void
    polyBuildEqualityMatrix( MatrixXT& Aeq,
                             VectorXT& beq,
                             Matrix<ValueType, 2, Dynamic> const & waypoints,
                             double direction = std::numeric_limits<double>::infinity() )
    {
        constexpr size_t DirectionConstraintCount = 2;
        constexpr auto CoeffCount = Spline::CoeffCount;

        auto const waypointCount = waypoints.cols();
        auto const segmentCount = waypointCount - 1;

        bool const supportDirection = std::isfinite(direction);

        Aeq.resize( 4*(segmentCount-1)                                          // Differential Constraints
                        + 4*(segmentCount)                                      // Equality Constraints
                        + (supportDirection ? DirectionConstraintCount : 0)     // Direction Constraint
                    , 2*CoeffCount*segmentCount );
        beq.resize(Aeq.rows(), 1);

        Aeq.setZero();
        beq.setZero();

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
            Aeq.block<1, CoeffCount>(row + 0, col) = leftTerms;
            Aeq.block<1, CoeffCount>(row + 1, col + CoeffCount) = leftTerms;
            beq.block<2, 1>(row, 0) = waypoints.block<2, 1>(0, indSeg);
            row += 2;

            // Right Equality (C0)
            Aeq.block<1, CoeffCount>(row + 0, col) = rightTerms;
            Aeq.block<1, CoeffCount>(row + 1, col + CoeffCount) = rightTerms;
            beq.block<2, 1>(row, 0) = waypoints.block<2, 1>(0, indSeg + 1);
            row += 2;

            if( (indSeg == 0) && supportDirection )
            {
                Aeq.bottomRows(DirectionConstraintCount).block<1, CoeffCount>(0, col) = dLeftTerms;
                Aeq.bottomRows(DirectionConstraintCount).block<1, CoeffCount>(1, col + CoeffCount) = dLeftTerms;

                Matrix<ValueType, 2, 1> vec;
                vec(0) = std::cos(direction);
                vec(1) = std::sin(direction);

                beq.tail(DirectionConstraintCount) = vec;
            }

            if( indSeg < segmentCount - 1 )
            {
                // Right Differential Continuity (C1)
                Aeq.block<1, CoeffCount>(row + 0, col) = dRightTerms;
                Aeq.block<1, CoeffCount>(row + 0, col + 2 * CoeffCount) = -dRightTerms;
                Aeq.block<1, CoeffCount>(row + 1, col + CoeffCount) = dRightTerms;
                Aeq.block<1, CoeffCount>(row + 1, col + 3 * CoeffCount) = -dRightTerms;
                row += 2;

                // Right Differential Continuity (C2)
                Aeq.block<1, CoeffCount>(row + 0, col) = ddRightTerms;
                Aeq.block<1, CoeffCount>(row + 0, col + 2 * CoeffCount) = -ddRightTerms;
                Aeq.block<1, CoeffCount>(row + 1, col + CoeffCount) = ddRightTerms;
                Aeq.block<1, CoeffCount>(row + 1, col + 3 * CoeffCount) = -ddRightTerms;
                row += 2;
            }
        }
    }


    void
    polyBuildLSQMatrix( MatrixXT& A,
                        VectorXT& b,
                        Matrix<ValueType, 2, Dynamic> const & waypoints )
    {
        constexpr auto CoeffCount = Spline::CoeffCount;
        constexpr auto IntermediateSampleCount = 9;
        const auto waypointCount = waypoints.cols();
        const auto segmentCount = waypointCount - 1;
        const auto sampleCount = segmentCount * IntermediateSampleCount;
        using SampleRow = Matrix<ValueType, 1, IntermediateSampleCount>;

        A.resize(2*sampleCount, 2*CoeffCount*segmentCount);
        b.resize(A.rows());

        A.setZero();
        b.setZero();

        PolySpace powers = PolySpace::LinSpaced(0, Spline::PolyOrder);
        size_t row = 0;
        size_t colX = 0, colY = CoeffCount;
        for(auto si = 0; si < segmentCount; ++si)
        {
            SampleRow sampleX = Matrix<ValueType, 1, IntermediateSampleCount+2>::LinSpaced(waypoints(0, si), waypoints(0, si+1)).middleCols(1, IntermediateSampleCount);
            SampleRow sampleY = Matrix<ValueType, 1, IntermediateSampleCount+2>::LinSpaced(waypoints(1, si), waypoints(1, si+1)).middleCols(1, IntermediateSampleCount);
            SampleRow sampleLambda = Matrix<ValueType, 1, IntermediateSampleCount+2>::LinSpaced(si * 1.0 / segmentCount, (si+1) * 1.0 / segmentCount).middleCols(1, IntermediateSampleCount);

            for(auto i = 0; i < sampleLambda.cols(); ++i)
            {
                PolySpace terms = PolySpace::Constant(sampleLambda(i));
                terms = pow(terms, powers);

                // Equality at the lambda
                A.block<1, CoeffCount>(row + 0, colX) = terms;
                A.block<1, CoeffCount>(row + 1, colY) = terms;
                b(row + 0) = sampleX(i);
                b(row + 1) = sampleY(i);

                row += 2;
            }

            colX += 2 * CoeffCount;
            colY += 2 * CoeffCount;
        }
    }

    /**
     * Splines a 5th order polynomial through the waypoints that ensures C2
     * conditions are met as well as a direction requirement at lambda=0
     *
     * The system is normally underdetermined so polyfit normally calculates the
     * right pseduo-inverse and performs a full pivot LU decomposition in order to
     * solve for the coefficients.
     *
     * Implements https://stanford.edu/class/ee103/lectures/constrained-least-squares/constrained-least-squares_slides.pdf
     */
    Matrix<ValueType, Dynamic, Spline::CoeffCount>
    polyfit(Matrix<ValueType, 2, Dynamic> waypoints, double direction = std::numeric_limits<double>::infinity() )
    {
        constexpr auto CoeffCount = Spline::CoeffCount;

        auto const waypointCount = waypoints.cols();
        auto const segmentCount = waypointCount - 1;

        MatrixXT Aeq, A;
        VectorXT beq, b;

        polyBuildSoftEqualityMatrix(Aeq, beq, waypoints, direction);
        polyBuildLSQMatrix(A, b, waypoints);

        // Moore-Penrose (Right) Pseudo-Inverse
        //c = Aeq.transpose() * (Aeq * Aeq.transpose()).fullPivLu().solve(beq);

        // Build full solution matrix F * (x, z)' = (2A'b, beq)'
        MatrixXT F( A.cols() + Aeq.rows(),
                    A.cols() + Aeq.rows() );
        VectorXT v( F.rows() );

        F.setZero();
        v.setZero();

        F.block(0, 0, A.cols(), A.cols()) = 2 * A.transpose() * A;
        F.block(0, A.cols(), Aeq.cols(), Aeq.rows()) = Aeq.transpose();
        F.block(A.cols(), 0, Aeq.rows(), Aeq.cols()) = Aeq;

        v.head(A.cols()) = 2 * A.transpose() * b;
        v.tail(beq.rows()) = beq;

        VectorXT d = F.fullPivHouseholderQr().solve(v);
        VectorXT c = d.head(Aeq.cols());
        /*
        MatrixXT c(Aeq.cols(), 1);

        // Moore-Penrose (Right) Pseudo-Inverse
        c = Aeq.transpose() * (Aeq * Aeq.transpose()).fullPivLu().solve(beq);
        */
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
        basis.col(0) = tangent / speed;

        // Rotate the first column by 90 degrees (right handed frame) to form the
        // full basis for R^2
        basis.col(1) = basis.col(0).reverse();
        basis(0, 1) = -basis(0, 1);
    }

    if(derivative == 1)
    {
        basis.col(0) = -basis.col(1);
        basis.col(1) = tangent / speed;
        return basis;
    }

    if(derivative >= 1)
    {
        Matrix<ValueType, 2, 1> accel = (*this)(parameter, 2);
        Matrix<ValueType, 2, 2> tframe = Matrix<ValueType, 2, 2>::Constant(0);
        //tframe.col(0) = -basis.col(0) / tangent.squaredNorm();
        //tframe.col(0) += Matrix<ValueType, 2, 1>::Constant(1.0 / speed);
        //basis.col(0) = tframe * accel;

        basis.col(0) = ( accel - accel.dot(basis.col(0)) * basis.col(0) ) / speed;
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

    ValueType lambda = lambdaStar;
    ValueType alpha = static_cast<ValueType>(1.0);
    alpha = std::min(alpha, alpha / speed(lambda));

    ValueType improvement = 1;
    int i = 500;
    while(std::abs(gradient(lambda)) > 1e-6 && (i > 0))
    {
        ValueType estimate = lambda - alpha * gradient(lambda);
        improvement = error(estimate) - error(lambda);
        if( improvement < 0 )
        {
            lambda = estimate;
            alpha *= 1.2;
        }
        else
            alpha *= 0.7;
        --i;
    }

    return std::max(std::min(lambda, static_cast<ValueType>(1.0)), static_cast<ValueType>(0.0));
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
