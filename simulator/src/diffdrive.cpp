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

DriveController::DriveController(   DriveSystem & driver,
                                    Matrix<ValueType, 2, 1> goal,
                                    ValueType goalTol )
  : mSystem(driver),
    mGoal(goal),
    mGoalRadius(goalTol),
    mPlanIndex(0),
    mOperatingLambda(0)
{
    replan();
}


Matrix<DriveSystem::ValueType, 3, 1> const & DriveSystem::state() const { return mState; }
Matrix<DriveSystem::ValueType, 3, 1> const & DriveSystem::flow() const { return mFlow; }
double DriveSystem::time() const { return mTime; }

void DriveController::operator() (StateType x, StateType& dxdt, double t)
{
    mSystem.step(x, dxdt, t, genTurnControl(x, t), genSpeedControl(x, t));
}

std::ostream& DriveController::printHeaders(std::ostream& stream) const
{
    stream << "time, x_1, x_2, x_3, dx_1, dx_2, dx_3, plan_ind";
    return printSpecificHeaders(stream);
}

std::ostream& DriveController::printData(std::ostream& stream) const
{
    auto state = mSystem.state();
    auto flow = mSystem.flow();
    stream << mSystem.time() << "," << state[0] << "," << state[1] << "," << state[2]
                             << "," << flow[0] << "," << flow[1] << "," << flow[2]
                             << "," << mPlanIndex;
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

void DriveController::updateKnownWorld(std::list< std::pair<World::Box, int> > const & obstacles)
{
    auto newList = obstacles;
    auto newEnd = std::remove_if(   newList.begin(),
                                    newList.end(),
                                    [this](std::pair<World::Box, int> a) -> bool {
                                        return mKnownBoxes.find(a.second) != mKnownBoxes.end();
                                    } );

    newList.erase(newEnd, newList.end());

    mKnownWorld.insert(newList.begin(), newList.end());

    for(auto&& pair : newList)
        mKnownBoxes.insert(pair.second);

    if( !newList.empty() )
        replan();
}

void DriveController::replan()
{
    auto state = mSystem.state();

    Spline spline;
    std::vector< std::pair< World::Box, int > > setCollisions;
    do
    {
        setCollisions.clear();

        auto path = planRRT( mKnownWorld,
                             World::XMax,
                             World::YMax,
                             World::Point(state[0], state[1]),
                             state[2],
                             World::Point(mGoal[0], mGoal[1]),
                             mGoalRadius );

        /*
        DEBUG: Just to output the path..incase we think we have a crap path generated.
        for(auto&& p : path)
        {
            std::cout << ";" << "(" << p.get<0>() << "," << p.get<1>() << ")";
        }
        std::cout << std::endl;
        */

        Matrix<ValueType, 2, Dynamic> waypoints(2, path.size());
        for(int i = 0; i < waypoints.cols(); ++i)
        {
            waypoints.col(i)[0] = path.front().get<0>();
            waypoints.col(i)[1] = path.front().get<1>();
            path.pop_front();
        }

        spline = Spline(waypoints, state[2]);
        auto const & linApprox = spline.approximation();

        for(auto&& segment : linApprox)
        {
            mKnownWorld.query( geom::index::intersects(segment), std::back_inserter(setCollisions) );

            // Let's not waste time checking ...if there is a collision!
            if(!setCollisions.empty())
                break;
        }

    } while( !setCollisions.empty() );

    path(spline);
    operatingPoint(0.0);
    ++mPlanIndex;
}

void DriveController::path( Spline const & s ) { mPath = s; }
Spline const & DriveController::path() const { return mPath; }

void DriveController::operatingPoint( Spline::ValueType v ) { mOperatingLambda = v; }
Spline::ValueType DriveController::operatingPoint() const { return mOperatingLambda; }

std::ostream& std::operator << (std::ostream& s, DriveController const & c)
{
    return c.printData(s);
}
