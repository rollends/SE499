#ifndef ROSE499_SFCONTROL_HPP
#define ROSE499_SFCONTROL_HPP

#include <Eigen/Core>
#include "rose499/diffdrive.hpp"
#include "rose499/spline.hpp"

struct SerretFrenetController : public DriveController
{
    virtual DriveController::ValueType genTurnControl(DriveController::StateType x, double t) override;

private:
    Spline mPath;
    Eigen::Matrix<DriveController::ValueType, 2, 1> mXi;
};

#endif
