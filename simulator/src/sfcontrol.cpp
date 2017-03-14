#include "rose499/sfcontrol.hpp"

using namespace Eigen;
#include <iostream>
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
    double lambda = path().nearestPoint(h, mOperatingLambda);
    mOperatingLambda = lambda;

    // Formulate Linearized State
    Matrix2T xi;
    Matrix2T s = path()(lambda);
    Matrix22T frame = path().frame(lambda);
    Matrix22T dFrame = path().frame(lambda, 1);

    auto e1 = frame.block<2, 1>(0, 0);
    auto e2 = frame.block<2, 1>(0, 1);

    mXi(0) = e2.dot(h - s);
    mXi(1) = e2.dot(dh * f);

    auto dwdh = e1.transpose() / path().speed(lambda);
    auto de2dx = dFrame.block<2, 1>(0, 1) * dwdh * dh;

    auto lff = (de2dx * f).dot( dh * f );
    auto lgf = (   e2    ).dot( dh * df * g );
    auto ulin = -mXi(0) - sqrt(3) * mXi(1);
    auto ucontrol = (5 * ulin - lff) / lgf;

    assert(!isnan(lff) && !isinf(lff));
    assert(!isnan(lgf) && !isinf(lgf));
    assert(!isnan(ulin) && !isinf(ulin));
    assert(!isnan(ucontrol) && !isinf(ucontrol));

    return ucontrol;
}

SerretFrenetController::SerretFrenetController(DriveSystem& sys)
  : mOperatingLambda(0), DriveController(sys) { }


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

Spline::ValueType SerretFrenetController::operatingPoint() const { return mOperatingLambda; }
Matrix<DriveController::ValueType, 2, 1> const & SerretFrenetController::linearizedState() const { return mXi; }
