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
    double lambda;

    // Formulate Linearized State
    Matrix2T xi;
    Matrix2T s = mPath(lambda);
    Matrix22T frame = mPath.frame(lambda);
    Matrix22T dFrame = mPath.frame(lambda, 1);

    auto e1 = frame.block<2, 1>(0, 0);
    auto e2 = frame.block<2, 1>(0, 1);

    mXi(0) = e2.dot(h - s);
    mXi(1) = e2.dot(dh * f);

    auto dwdh = e1.transpose() / mPath.speed(lambda);
    auto de2dx = dFrame.block<2, 1>(0, 1) * dwdh * dh;

    auto lff = (de2dx * f).dot( dh * f );
    auto lgf = (   e2    ).dot( dh * df * g );
    auto ulin = 0;

    return (ulin - lff) / lgf;
}
