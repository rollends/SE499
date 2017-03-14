#include <cmath>

#include "rose499/diffdrive.hpp"

using namespace Eigen;

void DriveSystem::step( StateType x,
                        StateType& dxdt,
                        double t,
                        ValueType turn,
                        ValueType speed )
{
    mTime = t;
    mState = Map<Matrix<ValueType, 3, 1>>(x.data());
    mFlow << speed * std::cos(mState(2)),
             speed * std::sin(mState(2)),
             turn;
    Map<Matrix<ValueType, 3, 1>>(dxdt.data()) = mFlow;
}

DriveController::DriveController(DriveSystem & driver)
  : mSystem(driver) { }


Matrix<DriveSystem::ValueType, 3, 1> const & DriveSystem::state() const { return mState; }
Matrix<DriveSystem::ValueType, 3, 1> const & DriveSystem::flow() const { return mFlow; }
double DriveSystem::time() const { return mTime; }

void DriveController::operator() (StateType x, StateType& dxdt, double t)
{
    mSystem.step(x, dxdt, t, genTurnControl(x, t), genSpeedControl(x, t));
}

std::ostream& DriveController::printHeaders(std::ostream& stream) const
{
    stream << "time, x_1, x_2, x_3, dx_1, dx_2, dx_3";
    return printSpecificHeaders(stream);
}

std::ostream& DriveController::printData(std::ostream& stream) const
{
    auto state = mSystem.state();
    auto flow = mSystem.flow();
    stream << mSystem.time() << "," << state[0] << "," << state[1] << "," << state[2]
                             << "," << flow[0] << "," << flow[1] << "," << flow[2];
    return printSpecificData(stream);
}

std::ostream& DriveController::printSpecificHeaders(std::ostream& s) const { return s; }
std::ostream& DriveController::printSpecificData(std::ostream& s) const { return s; }

DriveController::ValueType DriveController::genTurnControl(StateType, double)
{
    return 0;
}

DriveController::ValueType DriveController::genSpeedControl(StateType, double)
{
    return 1;
}

void DriveController::updateKnownWorld(std::vector<World::Box> const & obstacles)
{
    mKnownWorld.clear();
    mKnownWorld.insert(obstacles.begin(), obstacles.end());
}

void DriveController::path( Spline s ) { mPath = s; }
Spline const & DriveController::path() const { return mPath; }

std::ostream& std::operator << (std::ostream& s, DriveController const & c)
{
    return c.printData(s);
}
