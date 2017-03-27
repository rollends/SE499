#include "rose499/trackcontrol.hpp"

using namespace Eigen;

constexpr double TrackingController::LeadDistance;
constexpr double TrackingController::LambdaRate;

DriveController::ValueType TrackingController::genTurnControl(DriveController::StateType x, double t )
{
    using Matrix3T = Matrix<ValueType, 3, 1>;
    using Matrix33T = Matrix<ValueType, 3, 3>;
    using Matrix2T = Matrix<ValueType, 2, 1>;
    using Matrix23T = Matrix<ValueType, 2, 3>;
    using Matrix22T = Matrix<ValueType, 2, 2>;

    t = t - lastPlanTime;

    ValueType ucontrol = 0;
    Matrix<ValueType, 2, 2> R;
    R << 0, -1,
         1, 0;

    // Calculate System Dynamics (i.e. no control)
    Matrix2T tau, h, sigma, dsigma;
    Matrix22T K = -0.5*Matrix22T::Identity();

    tau << std::cos(x[2]), std::sin(x[2]);
    h = Map<Matrix2T>(x.data());
    sigma = path()(t / LambdaRate);
    dsigma = path()(t / LambdaRate, 1);

    Matrix2T error = h - sigma;

    ucontrol = (R*tau).dot(K*error + dsigma / LambdaRate) / LeadDistance;

    return ucontrol;
}

DriveController::ValueType TrackingController::genSpeedControl(DriveController::StateType x, double t )
{
    using Matrix3T = Matrix<ValueType, 3, 1>;
    using Matrix33T = Matrix<ValueType, 3, 3>;
    using Matrix2T = Matrix<ValueType, 2, 1>;
    using Matrix23T = Matrix<ValueType, 2, 3>;
    using Matrix22T = Matrix<ValueType, 2, 2>;
    ValueType vcontrol = 0;

    t = t - lastPlanTime;

    // Calculate System Dynamics (i.e. no control)
    Matrix2T tau, h, sigma, dsigma;
    Matrix22T K = -0.5*Matrix22T::Identity();

    tau << std::cos(x[2]), std::sin(x[2]);
    h = Map<Matrix2T>(x.data());
    sigma = path()(t / LambdaRate);
    dsigma = path()(t / LambdaRate, 1);

    Matrix2T error = h - sigma;

    vcontrol = tau.dot(K*error + dsigma / LambdaRate);

    return vcontrol;
}

TrackingController::TrackingController(DriveSystem& sys, Eigen::Matrix<ValueType, 2, 1> goal, ValueType goalRadius)
  : DriveController(sys, goal, goalRadius) { }


std::ostream& TrackingController::printSpecificHeaders(std::ostream& s) const
{
    return s << ", lambdaStar, spline_ind, sigma_1, sigma_2";
}

std::ostream& TrackingController::printSpecificData(std::ostream& s) const
{
    auto lambda = system().time() - lastPlanTime;
    auto sigma = path()(lambda / LambdaRate);
    return s << ", " << lambda / LambdaRate
             << ", " << path().splineIndexUsed(lambda / LambdaRate)
             << ", " << sigma[0]
             << "," << sigma[1];
}
