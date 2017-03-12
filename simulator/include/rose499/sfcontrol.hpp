#ifndef ROSE499_SFCONTROL_HPP
#define ROSE499_SFCONTROL_HPP

#include <Eigen/Core>
#include "rose499/diffdrive.hpp"
#include "rose499/spline.hpp"

struct SerretFrenetController : public DriveController
{
    SerretFrenetController(DriveSystem&);
    virtual DriveController::ValueType genTurnControl(DriveController::StateType x, double t) override;

    Spline::ValueType operatingPoint() const;
    Eigen::Matrix<DriveController::ValueType, 2, 1> linearizedState() const;

private:
    Spline::ValueType mOperatingLambda;
    Eigen::Matrix<DriveController::ValueType, 2, 1> mXi;
};

#endif
