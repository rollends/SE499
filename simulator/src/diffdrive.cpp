#include <cmath>

#include "rose499/diffdrive.hpp"

using namespace Eigen;

void DriveSystem::step( StateType x,
                        StateType& dxdt,
                        double t,
                        ValueType turn,
                        ValueType speed )
{
    mState = Map<Matrix<ValueType, 3, 1>>(x.data());
    mFlow << speed * std::cos(mState(2)),
             speed * std::sin(mState(2)),
             turn;
    Map<Matrix<ValueType, 3, 1>>(dxdt.data()) = mFlow;
}

DriveController::DriveController(DriveSystem & driver)
  : mSystem(driver) { }

void DriveController::operator() (StateType x, StateType& dxdt, double t)
{
    mSystem.step(x, dxdt, t, genTurnControl(x, t), genSpeedControl(x, t));
}

DriveController::ValueType DriveController::genTurnControl(StateType, double)
{
    return 0;
}

DriveController::ValueType DriveController::genSpeedControl(StateType, double)
{
    return 1;
}
