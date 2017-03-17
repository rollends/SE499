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
    SerretFrenetController sfcontrol(robot, Eigen::Vector2d(95, 95), 3);
    sfcontrol.printHeaders(std::cout) << std::endl;

    // Create World
    int i = 0;
    for(i = 0; i < 6; i++)
        world.addBoxObstacle(i, 7 + 12 * i, 5, 5, 20);
    world.addBoxObstacle(i++, 0, 50, 80, 10);
    world.addBoxObstacle(i++, 90, 50, 10, 10);

    // Plan every 200ms
    constexpr double PlanTime = 0.050;
    double planningClock = 0.0;

    // Prime Planning Algo
    World::Polygon viewFrustum;
    geom::append(viewFrustum.outer(), World::Point(x[0], x[1]));
    {
        Eigen::Vector2d dir(std::cos(x[3]), std::sin(x[3]));
        Rotation2D<double> right(std::atan2(-1, 1));
        Rotation2D<double> left(std::atan2(1, 1));
        auto viewRight = right * dir;
        auto viewLeft = left * dir;

        geom::append(viewFrustum.outer(), World::Point(x[0] + 7 * viewRight[0], x[1] + 7 * viewRight[1]));
        geom::append(viewFrustum.outer(), World::Point(x[0] + 7 * viewLeft[0], x[1] + 7 * viewLeft[1]));
    }
    geom::append(viewFrustum.outer(), World::Point(x[0], x[1]));
    sfcontrol.updateKnownWorld(world.observeWorld(viewFrustum));

    do
    {

        // If we fail, retry with smaller step size (hopefully it converges, eh?)
        if(odeint::fail == ode45.try_step(std::reference_wrapper<DriveController>(sfcontrol), x, T, dt))
            continue;

        dt = std::min( dt, 0.01 );

        // If we pass store the system state to file.
        std::cout << sfcontrol << std::endl;

        if((planningClock += dt) >= PlanTime)
        {
            planningClock = 0.0;

            // Observe the world, and make sure there are no changes before replanning.
            World::Polygon viewFrustum;
            geom::append(viewFrustum.outer(), World::Point(x[0], x[1]));
            {
                Eigen::Vector2d dir(std::cos(x[3]), std::sin(x[3]));
                Rotation2D<double> right(std::atan2(-1, 1));
                Rotation2D<double> left(std::atan2(1, 1));
                auto viewRight = right * dir;
                auto viewLeft = left * dir;

                geom::append(viewFrustum.outer(), World::Point(x[0] + 10 * viewRight[0], x[1] + 10 * viewRight[1]));
                geom::append(viewFrustum.outer(), World::Point(x[0] + 10 * viewLeft[0], x[1] + 10 * viewLeft[1]));
            }
            geom::append(viewFrustum.outer(), World::Point(x[0], x[1]));
            auto observedWorld = world.observeWorld(viewFrustum);
            sfcontrol.updateKnownWorld(observedWorld);
        }
    } while( (T < Tfinal) || (sfcontrol.operatingPoint() < 1) );

    return 0;
}
