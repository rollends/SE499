#include <cmath>

#include "rose499/diffdrive.hpp"

DriveSystem::DriveSystem(DriveController& c)
  : control(c) { }

DriveSystem::DriveSystem(DriveSystem const & sys)
  : control(sys.control) { }

void DriveSystem::operator() (StateType x,
                              StateType& dxdt,
                              double t)
{
    auto speed = control.genSpeedControl(x, t);
    auto turn = control.genTurnControl(x, t);

    dxdt[0] = speed * std::cos(x[3]);
    dxdt[1] = speed * std::sin(x[3]);
    dxdt[2] = turn;
}

DriveController::ValueType DriveController::genTurnControl(StateType, double)
{
    return 0;
}

DriveController::ValueType DriveController::genSpeedControl(StateType, double)
{
    return 1;
}
