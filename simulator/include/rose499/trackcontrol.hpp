#ifndef ROSE499_TRACKCONTROL_HPP
#define ROSE499_TRACKCONTROL_HPP

#include <Eigen/Core>
#include "rose499/diffdrive.hpp"
#include "rose499/spline.hpp"

struct TrackingController : public DriveController
{
    static constexpr double LeadDistance = 1;
    static constexpr double LambdaRate = 600.0;

    TrackingController(DriveSystem&, Eigen::Matrix<ValueType, 2, 1> goal, ValueType goalRadius);
    virtual DriveController::ValueType genTurnControl(DriveController::StateType x, double t) override;
    virtual DriveController::ValueType genSpeedControl(DriveController::StateType x, double t) override;

protected:
    std::ostream& printSpecificHeaders(std::ostream& s) const override;
    std::ostream& printSpecificData(std::ostream& s) const override;
};

#endif
