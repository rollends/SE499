#include "rose499/sylvester.hpp"

using namespace Eigen;

DriveController::ValueType SylvesterController::genTurnControl(DriveController::StateType x, double t )
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
    operatingPoint(path().nearestPoint(h, operatingPoint()));

    // Nearest Polynomial
    int indPoly = path().splineIndexUsed(operatingPoint());
    Matrix<ValueType, 2, Spline::CoeffCount> poly = path().poly().block<2, Spline::CoeffCount>(2 * indPoly, 0);

    auto J = quinticJacobian(   poly(0, 5),
                                poly(0, 4),
                                poly(0, 3),
                                poly(0, 2),
                                poly(0, 1),
                                poly(0, 0),
                                poly(1, 5),
                                poly(1, 4),
                                poly(1, 3),
                                poly(1, 2),
                                poly(1, 1),
                                poly(1, 0),
                                x[0],
                                x[1] );
    auto H = quinticHessian(    poly(0, 5),
                                poly(0, 4),
                                poly(0, 3),
                                poly(0, 2),
                                poly(0, 1),
                                poly(0, 0),
                                poly(1, 5),
                                poly(1, 4),
                                poly(1, 3),
                                poly(1, 2),
                                poly(1, 1),
                                poly(1, 0),
                                x[0],
                                x[1] ) ;
    // Formulate Linearized State
    mXi(0) = quinticLevelSet(   poly(0, 5),
                                poly(0, 4),
                                poly(0, 3),
                                poly(0, 2),
                                poly(0, 1),
                                poly(0, 0),
                                poly(1, 5),
                                poly(1, 4),
                                poly(1, 3),
                                poly(1, 2),
                                poly(1, 1),
                                poly(1, 0),
                                x[0],
                                x[1] );
    mXi(1) = J.dot(dh * f);

    long double lff = (dh * f).transpose() * H * (dh * f);
    long double lgf = J.dot(dh * df.col(2));

    constexpr ValueType damping = 2.0;
    const ValueType g1 = -std::sqrt(10); // choosing g1 to bound xi1 by 0.1
    const ValueType g2 = -2 * damping * std::pow(10, 0.25);

    auto ucontrol = g1 * mXi(0)/ lgf + g2 * mXi(1) / lgf - (lff / lgf);

    assert(!isnan(lff) && !isinf(lff));
    assert(!isnan(lgf) && !isinf(lgf));
    assert(!isnan(ucontrol) && !isinf(ucontrol));

    return ucontrol;
}

SylvesterController::SylvesterController(DriveSystem& sys, Eigen::Matrix<ValueType, 2, 1> goal, ValueType goalRadius)
  : DriveController(sys, goal, goalRadius) { }


std::ostream& SylvesterController::printSpecificHeaders(std::ostream& s) const
{
    return s << ", lambdaStar, spline_ind, sigma_1, sigma_2, xi_1, xi_2";
}

std::ostream& SylvesterController::printSpecificData(std::ostream& s) const
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

Matrix<DriveController::ValueType, 2, 1> const & SylvesterController::linearizedState() const { return mXi; }
