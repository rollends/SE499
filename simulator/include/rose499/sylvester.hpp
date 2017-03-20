#ifndef ROSE499_SYLVESTERCONTROL_HPP
#define ROSE499_SYLVESTERCONTROL_HPP

#include <Eigen/Core>
#include "rose499/diffdrive.hpp"
#include "rose499/spline.hpp"

struct SylvesterController : public DriveController
{
    SylvesterController(DriveSystem&, Eigen::Matrix<ValueType, 2, 1> goal, ValueType goalRadius);
    virtual DriveController::ValueType genTurnControl(DriveController::StateType x, double t) override;

    Eigen::Matrix<DriveController::ValueType, 2, 1> const & linearizedState() const;

protected:
    std::ostream& printSpecificHeaders(std::ostream& s) const override;
    std::ostream& printSpecificData(std::ostream& s) const override;

private:
    Eigen::Matrix<DriveController::ValueType, 2, 1> mXi;
};

Eigen::Matrix<SimulatorTypes::ValueType, 2, 2> quinticHessian(  SimulatorTypes::ValueType A1,
                                                                SimulatorTypes::ValueType A2,
                                                                SimulatorTypes::ValueType A3,
                                                                SimulatorTypes::ValueType A4,
                                                                SimulatorTypes::ValueType A5,
                                                                SimulatorTypes::ValueType A6,
                                                                SimulatorTypes::ValueType B1,
                                                                SimulatorTypes::ValueType B2,
                                                                SimulatorTypes::ValueType B3,
                                                                SimulatorTypes::ValueType B4,
                                                                SimulatorTypes::ValueType B5,
                                                                SimulatorTypes::ValueType B6,
                                                                SimulatorTypes::ValueType x1,
                                                                SimulatorTypes::ValueType x2 );

Eigen::Matrix<SimulatorTypes::ValueType, 1, 2> quinticJacobian( SimulatorTypes::ValueType A1,
                                                                SimulatorTypes::ValueType A2,
                                                                SimulatorTypes::ValueType A3,
                                                                SimulatorTypes::ValueType A4,
                                                                SimulatorTypes::ValueType A5,
                                                                SimulatorTypes::ValueType A6,
                                                                SimulatorTypes::ValueType B1,
                                                                SimulatorTypes::ValueType B2,
                                                                SimulatorTypes::ValueType B3,
                                                                SimulatorTypes::ValueType B4,
                                                                SimulatorTypes::ValueType B5,
                                                                SimulatorTypes::ValueType B6,
                                                                SimulatorTypes::ValueType x1,
                                                                SimulatorTypes::ValueType x2 );

SimulatorTypes::ValueType                      quinticLevelSet( SimulatorTypes::ValueType A1,
                                                                SimulatorTypes::ValueType A2,
                                                                SimulatorTypes::ValueType A3,
                                                                SimulatorTypes::ValueType A4,
                                                                SimulatorTypes::ValueType A5,
                                                                SimulatorTypes::ValueType A6,
                                                                SimulatorTypes::ValueType B1,
                                                                SimulatorTypes::ValueType B2,
                                                                SimulatorTypes::ValueType B3,
                                                                SimulatorTypes::ValueType B4,
                                                                SimulatorTypes::ValueType B5,
                                                                SimulatorTypes::ValueType B6,
                                                                SimulatorTypes::ValueType x1,
                                                                SimulatorTypes::ValueType x2 );
#endif
