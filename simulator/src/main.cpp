#include <array>
#include <boost/numeric/odeint.hpp>
#include <Eigen/Geometry>
#include <iostream>

#include "rose499/sfcontrol.hpp"
#include "rose499/sylvester.hpp"

using namespace Eigen;
using namespace SimulatorTypes;

int main(int argc, char* argv[])
{
    namespace odeint = boost::numeric::odeint;

    constexpr ValueType Tfinal = 15.0;
    ValueType T = 0.0;
    ValueType dt = 0.0001;

    // Initial Condition
    DriveSystem::StateType x{ 0, 0, std::atan2(1.0, 1.0) };
    auto ode45 = odeint::make_controlled( 1.0e-12 , 1.0e-3 , odeint::runge_kutta_fehlberg78<DriveSystem::StateType, ValueType>() );
    //auto ode45 = odeint::make_controlled( 1.0e-6 , 1.0e-6 , odeint::runge_kutta_cash_karp54<DriveSystem::StateType>() );
    //auto ode45 = odeint::make_controlled( 1.0e-12 , 1.0e-3 , odeint::runge_kutta_fehlberg78<DriveSystem::StateType>() );
    //odeint::bulirsch_stoer<DriveSystem::StateType> ode45(1e-12, 1e-12);

    World world;
    DriveSystem robot;

    const Vector2T goal(95, 95);
    const ValueType goalRadius = 3;

    //SylvesterController sfcontrol(robot, goal, goalRadius);
    SerretFrenetController sfcontrol(robot, goal, goalRadius);
    sfcontrol.printHeaders(std::cout) << std::endl;

    // Create World
    int i = 0;
    for(i = 0; i < 6; i++)
        world.addBoxObstacle(i, 10 + 10 * i, 5, 5, 20);
    world.addBoxObstacle(++i, 0, 50, 80, 10);
    world.addBoxObstacle(++i, 90, 50, 10, 10);

    // Plan every 200ms
    constexpr double PlanTime = 0.050;
    constexpr ValueType MaxStep = 0.001;
    constexpr double PrintTime = 0.05;
    double planningClock = 0.0;
    double printingClock = 0.0;

    // Prime Planning Algo
    sfcontrol.updateKnownWorld(world.observeWorld(robot.viewCone()));
    do
    {
        // If we fail, retry with smaller step size (hopefully it converges, eh?)
        if(odeint::fail == ode45.try_step(std::reference_wrapper<DriveController>(sfcontrol), x, T, dt))
            continue;

        // Update current state of the robot, and enforce constraints.
        robot.state(Map<Matrix<DriveSystem::ValueType, 3, 1>>(x.data()));
        robot.time(T);

        dt = std::min( dt, MaxStep );

        // If we pass store the system state to file.
        if((printingClock += dt) >= PrintTime)
        {
            printingClock = 0;
            std::cout << sfcontrol << std::endl;
        }

        if(
            (planningClock += dt) >= PlanTime
         || (sfcontrol.operatingPoint() >= 1.0)
         || sfcontrol.hasDiverged()
          )
        {
            planningClock = 0.0;

            // Observe the world and update our 'known' (Estimated) world state
            auto observedWorld = world.observeWorld(robot.viewCone());
            sfcontrol.updateKnownWorld(observedWorld, sfcontrol.hasDiverged());
        }
    } while( (sfcontrol.goal() - Vector2T(x[0], x[1])).norm() > goalRadius );

    return 0;
}
