#ifndef ROSE499_SFCONTROL_HPP
#define ROSE499_SFCONTROL_HPP

#include <Eigen/Core>
#include "rose499/diffdrive.hpp"
#include "rose499/spline.hpp"

struct SerretFrenetController : public DriveController
{
    SerretFrenetController(DriveSystem&, Eigen::Matrix<ValueType, 2, 1> goal, ValueType goalRadius);
    virtual DriveController::ValueType genTurnControl(DriveController::StateType x, double t) override;

    Eigen::Matrix<DriveController::ValueType, 2, 1> const & linearizedState() const;
    virtual bool hasDiverged() const override;

protected:
    std::ostream& printSpecificHeaders(std::ostream& s) const override;
    std::ostream& printSpecificData(std::ostream& s) const override;

private:
    Eigen::Matrix<DriveController::ValueType, 2, 1> mXi;
};

#endif
