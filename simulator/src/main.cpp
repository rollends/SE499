#include <array>
#include <boost/numeric/odeint.hpp>
#include <Eigen/Geometry>
#include <iostream>

#include "rose499/sfcontrol.hpp"

using namespace Eigen;

int main(int argc, char* argv[])
{
    namespace odeint = boost::numeric::odeint;

    constexpr double Tfinal = 15.0;
    double T = 0.0;
    double dt = 0.001;

    // Initial Condition
    DriveSystem::StateType x{ 0, 0, std::atan2(1.0, 1.0) };
    auto ode45 = odeint::make_controlled( 1.0e-12 , 1.0e-12 , odeint::runge_kutta_fehlberg78<DriveSystem::StateType>() );

    World world;
    DriveSystem robot;

    const Eigen::Vector2d goal(95, 95);
    const double goalRadius = 3;

    SerretFrenetController sfcontrol(robot,goal,goalRadius);
    sfcontrol.printHeaders(std::cout) << std::endl;

    // Create World
    int i = 0;
    for(i = 0; i < 6; i++)
        world.addBoxObstacle(i, 10 + 10 * i, 5, 5, 20);
    world.addBoxObstacle(++i, 0, 50, 80, 10);
    world.addBoxObstacle(++i, 90, 50, 10, 10);

    // Plan every 200ms
    constexpr double PlanTime = 0.050;
    double planningClock = 0.0;

    // Prime Planning Algo
    auto hasDiverged = [&sfcontrol](){ return std::abs(sfcontrol.linearizedState()[1]) >= 0.7; };
    sfcontrol.updateKnownWorld(world.observeWorld(robot.viewCone()));
    do
    {
        // If we fail, retry with smaller step size (hopefully it converges, eh?)
        if(odeint::fail == ode45.try_step(std::reference_wrapper<DriveController>(sfcontrol), x, T, dt))
            continue;

        // Update current state of the robot.
        robot.state(Map<Matrix<DriveSystem::ValueType, 3, 1>>(x.data()));
        robot.time(T);

        dt = std::min( dt, 0.001 );

        // If we pass store the system state to file.
        std::cout << sfcontrol << std::endl;

        if(
            (planningClock += dt) >= PlanTime
         || (sfcontrol.operatingPoint() >= 1.0)
         || hasDiverged()
          )
        {
            planningClock = 0.0;

            // Observe the world and update our 'known' (Estimated) world state
            auto observedWorld = world.observeWorld(robot.viewCone());
            sfcontrol.updateKnownWorld(observedWorld, hasDiverged());
        }
    } while( (sfcontrol.goal() - Eigen::Vector2d(x[0], x[1])).norm() > goalRadius );

    return 0;
}
