#include <array>
#include <boost/numeric/odeint.hpp>
#include <Eigen/Geometry>
#include <iostream>

#include "rose499/sfcontrol.hpp"
#include "rose499/sylvester.hpp"

using namespace Eigen;
using namespace SimulatorTypes;

void fillCampusWorld(World& world);

int main(int argc, char* argv[])
{
    namespace odeint = boost::numeric::odeint;

    ValueType T = 0.0;
    ValueType dt = 0.0001;

    // Initial Condition
    DriveSystem::StateType x{ 325, 95, std::atan2(1.0, -1.0) };
    auto ode45 = odeint::make_controlled( 1.0e-12 , 1.0e-12 , odeint::runge_kutta_cash_karp54<DriveSystem::StateType, ValueType>() );
    //auto ode45 = odeint::make_controlled( 1.0e-12 , 1.0e-12 , odeint::runge_kutta_cash_karp54<DriveSystem::StateType, ValueType>() );

    World world;
    DriveSystem robot;
    robot.time(T);
    robot.state(Map<Matrix<DriveSystem::ValueType, 3, 1>>(x.data()));

    const Vector2T goal(218, 600);
    const ValueType goalRadius = 5;

    //SylvesterController sfcontrol(robot, goal, goalRadius);
    SerretFrenetController sfcontrol(robot, goal, goalRadius);
    sfcontrol.printHeaders(std::cout) << std::endl;

    // Create World
    fillCampusWorld(world);

    // Plan every 200ms
    constexpr double PlanTime = 0.05;
    constexpr double PrintTime = 0.01;
    constexpr ValueType MaxStep = 0.005;
    double planningClock = 0.0;
    double printingClock = 0.0;

    // Prime Planning Algo
    sfcontrol.updateKnownWorld(world.observeWorld(robot.viewCone()), true);
    do
    {
        //std::cout << sfcontrol.path().poly() << std::endl;
        try
        {
            // If we fail, retry with smaller step size (hopefully it converges, eh?)
            if(odeint::fail == ode45.try_step(std::reference_wrapper<DriveController>(sfcontrol), x, T, dt))
                continue;
        }
        catch( InvalidPathException const & )
        {
            // Force replan.
            planningClock = 0.0;
            auto observedWorld = world.observeWorld(robot.viewCone());
            sfcontrol.updateKnownWorld(observedWorld, true);
        }

        // Update current state of the robot, and enforce constraints.
        robot.time(T);
        robot.state(Map<Matrix<DriveSystem::ValueType, 3, 1>>(x.data()));
        sfcontrol.operatingPoint(sfcontrol.path().nearestPoint(robot.state().head(2), sfcontrol.operatingPoint()));

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

void fillCampusWorld(World& world)
{
    // Scale , 64 units = 100m
    //      => 1 m /s  = 0.64 units / s
    int i = 0;
    world.addBoxObstacle(i++, 351, 429, 77, 200);
    world.addBoxObstacle(i++, 292, 528, 40, 35);
    world.addBoxObstacle(i++, 244, 525, 40, 35);
    world.addBoxObstacle(i++, 258, 488, 73, 35);
    world.addBoxObstacle(i++, 203, 530, 32, 31);
    world.addBoxObstacle(i++, 200, 487, 37, 36);
    world.addBoxObstacle(i++, 98, 542, 102, 52);
    world.addBoxObstacle(i++, 60, 471, 70, 60);
    world.addBoxObstacle(i++, 92, 404, 55, 61);
    world.addBoxObstacle(i++, 176, 404, 40, 61);
    world.addBoxObstacle(i++, 222, 455, 111, 12);
    world.addBoxObstacle(i++, 227, 393, 42, 30);
    world.addBoxObstacle(i++, 270, 384, 61, 67);
    world.addBoxObstacle(i++, 391, 347, 36, 80);
    world.addBoxObstacle(i++, 298, 294, 70, 79);
    world.addBoxObstacle(i++, 250, 320, 37, 53);
    world.addBoxObstacle(i++, 219, 330, 21, 57);
    world.addBoxObstacle(i++, 121, 360, 62, 41);
    world.addBoxObstacle(i++, 140, 294, 68, 64);
    world.addBoxObstacle(i++, 121, 244, 37, 42);
    world.addBoxObstacle(i++, 180, 224, 38, 38);
    world.addBoxObstacle(i++, 227, 257, 63, 50);
    world.addBoxObstacle(i++, 270, 227, 46, 28);
    world.addBoxObstacle(i++, 294, 264, 93, 20);
    world.addBoxObstacle(i++, 326, 217, 16, 51);
    world.addBoxObstacle(i++, 359, 190, 40, 69);
    world.addBoxObstacle(i++, 283, 197, 69, 16);
    world.addBoxObstacle(i++, 333, 164, 20, 47);
    world.addBoxObstacle(i++, 266, 136, 33, 53);
    world.addBoxObstacle(i++, 254, 200, 17, 15);
    world.addBoxObstacle(i++, 217, 139, 17, 51);
    world.addBoxObstacle(i++, 177, 154, 30, 37);
    world.addBoxObstacle(i++, 113, 193, 48, 23);
    world.addBoxObstacle(i++, 81, 140, 84, 47);
    world.addBoxObstacle(i++, 144, 86, 53, 52);
    world.addBoxObstacle(i++, 85, 49, 54, 79);
}
