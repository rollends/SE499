#include "rose499/sfcontrol.hpp"

using namespace Eigen;

DriveController::ValueType SerretFrenetController::genTurnControl(DriveController::StateType x, double t )
{
    using Matrix3T = Matrix<ValueType, 3, 1>;
    using Matrix33T = Matrix<ValueType, 3, 3>;
    using Matrix2T = Matrix<ValueType, 2, 1>;
    using Matrix23T = Matrix<ValueType, 2, 3>;
    using Matrix22T = Matrix<ValueType, 2, 2>;

    // Assume constant speed control from the default controller
    // TODO: This is a really bad dependency that should be accounted for
    auto speed = genSpeedControl(x, t);

    // Calculate System Dynamics (i.e. no control)
    Matrix3T f;
    Matrix2T h;
    Matrix23T dh;
    Matrix33T df;
    Matrix3T g;

    f << speed * std::cos(x[2]),
         speed * std::sin(x[2]),
         0;
    g << 0,
         0,
         1;
    h = Map<Matrix2T>(x.data());
    df << 0, 0, -speed * std::sin(x[2]),
          0, 0, speed * std::cos(x[2]),
          0, 0, 0;
    dh << 1, 0, 0,
          0, 1, 0;

    // Nearest Point
    double lambda = path().nearestPoint(h, operatingPoint());
    operatingPoint(lambda);

    if( std::abs(path().speed(lambda)) < 1e-6 )
        throw InvalidPathException();

    // Formulate Linearized State
    Matrix2T xi;
    Matrix2T s = path()(lambda);
    Matrix2T ds = path()(lambda, 1);
    Matrix2T dds = path()(lambda, 2);
    Matrix22T frame = path().frame(lambda);
    Matrix22T dFrame = path().frame(lambda, 1);

    auto e1 = frame.block<2, 1>(0, 0);
    auto e2 = frame.block<2, 1>(0, 1);

    mXi(0) = e2.dot(h - s);
    mXi(1) = e2.dot(dh * f);

    constexpr ValueType damping = 2.0;
    const ValueType g1 = -std::sqrt(10); // choosing g1 to bound xi1 by 0.1
    const ValueType g2 = -2 * damping * std::pow(10, 0.25);

    auto dwdh = e1.transpose() / path().speed(lambda);
    auto de2dx = dFrame.block<2, 1>(0, 1) * dwdh * dh;
    auto lff = (de2dx * f).dot( dh * f );
    //auto lff = (speed*speed - mXi(1)*mXi(1)) * std::abs(dds(0)*ds(1) - dds(1)*ds(0)) / std::pow(path().speed(lambda), 4.0);
    //auto lgf = (   e2    ).dot( dh * df * g );
    auto lgf = speed * sqrt(speed*speed - mXi(1)*mXi(1));
    //auto ulin = -mXi(0) - sqrt(2) * mXi(1);
    auto ulin = g1 * mXi(0) + g2 * mXi(1);
    auto ucontrol = (ulin - lff) / lgf;

    assert(!isnan(lff) && !isinf(lff));
    assert(!isnan(lgf) && !isinf(lgf));
    assert(!isnan(ulin) && !isinf(ulin));
    assert(!isnan(ucontrol) && !isinf(ucontrol));

    return ucontrol;
}

SerretFrenetController::SerretFrenetController(DriveSystem& sys, Eigen::Matrix<ValueType, 2, 1> goal, ValueType goalRadius)
  : DriveController(sys, goal, goalRadius) { }


bool SerretFrenetController::hasDiverged() const
{
    return std::abs(linearizedState()[1]) >= 0.7;
}

std::ostream& SerretFrenetController::printSpecificHeaders(std::ostream& s) const
{
    return s << ", lambdaStar, spline_ind, sigma_1, sigma_2, xi_1, xi_2";
}

std::ostream& SerretFrenetController::printSpecificData(std::ostream& s) const
{
    auto state = linearizedState();
    auto lambda = operatingPoint();
    auto sigma = path()(lambda);
    return s << ", " << lambda
             << ", " << path().splineIndexUsed(lambda)
             << ", " << sigma[0]
             << "," << sigma[1]
             << ", " << state[0]
             << "," << state[1];
}

Matrix<DriveController::ValueType, 2, 1> const & SerretFrenetController::linearizedState() const { return mXi; }
