#include <cmath>
#include <Eigen/Geometry>
#include "rose499/diffdrive.hpp"

using namespace Eigen;

void DriveSystem::step( StateType x,
                        StateType& dxdt,
                        double t,
                        ValueType turn,
                        ValueType speed )
{
    Map<Matrix<ValueType, 3, 1>>(dxdt.data()) << speed * std::cos(x[2]),
                                                 speed * std::sin(x[2]),
                                                 turn;
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
}


Matrix<DriveSystem::ValueType, 3, 1> const & DriveSystem::state() const { return mState; }
Matrix<DriveSystem::ValueType, 3, 1> const & DriveSystem::flow() const { return mFlow; }
double DriveSystem::time() const { return mTime; }
void DriveSystem::state(Eigen::Matrix<ValueType, 3, 1> x) { mState = x; }
void DriveSystem::flow(Eigen::Matrix<ValueType, 3, 1> f) { mFlow = f; }
void DriveSystem::time(double t) { mTime = t; }

World::Polygon DriveSystem::viewCone() const
{
    World::Polygon viewFrustum;
    geom::append(viewFrustum.outer(), World::Point(mState[0], mState[1]));
    {
        Vector2d dir(std::cos(mState[2]), std::sin(mState[2]));
        Rotation2D<double> right(std::atan2(-1, 1));
        Rotation2D<double> left(std::atan2(1, 1));
        Vector2d viewRight = right * dir;
        Vector2d viewLeft = left * dir;

        geom::append(viewFrustum.outer(), World::Point(mState[0] + 50 * viewLeft[0], mState[1] + 50 * viewLeft[1]));
        geom::append(viewFrustum.outer(), World::Point(mState[0] + 50 * viewRight[0], mState[1] + 50 * viewRight[1]));
    }
    geom::append(viewFrustum.outer(), World::Point(mState[0], mState[1]));
    return viewFrustum;
}

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

bool DriveController::updateKnownWorld(std::list< std::pair<World::Box, int> > const & obstacles, bool forceReplan)
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

    if( !newList.empty() || forceReplan ) {
        replan();
        return true;
    }
    return false;
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

DriveSystem const & DriveController::system() const { return mSystem; }

bool DriveController::hasDiverged() const { return false; }

Matrix<DriveController::ValueType, 2, 1> const & DriveController::goal() const { return mGoal; }

std::ostream& std::operator << (std::ostream& s, DriveController const & c)
{
    return c.printData(s);
}
