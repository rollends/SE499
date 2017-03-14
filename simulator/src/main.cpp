#include <array>
#include <boost/numeric/odeint.hpp>
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
    DriveSystem::StateType x{ 0, 1, 0 };//std::atan(1.0) * 4.0 / 2.0 };
    auto ode45 = odeint::make_controlled( 1.0e-12 , 1.0e-12 , odeint::runge_kutta_fehlberg78<DriveSystem::StateType>() );

    // Generate a spline
    Matrix<Spline::ValueType, 2, 5> waypoints;
    waypoints(0, 0) = 1; waypoints(1, 0) = 0;
    waypoints(0, 1) = 5; waypoints(1, 1) = 5;
    waypoints(0, 2) = 10; waypoints(1, 2) = 7;
    waypoints(0, 3) = 15; waypoints(1, 3) = 4;
    waypoints(0, 4) = 20; waypoints(1, 4) = 0;

    Spline spline(waypoints);

    DriveSystem robot;
    SerretFrenetController sfcontrol(robot);
    sfcontrol.path(spline);
    sfcontrol.printHeaders(std::cout) << std::endl;
    do
    {
        // If we fail
        if(odeint::fail == ode45.try_step(std::reference_wrapper<DriveController>(sfcontrol), x, T, dt))
            continue;

        dt = std::min( dt, 0.05 );

        // If we pass store the system state to file.
        std::cout << sfcontrol << std::endl;
    } while( (T < Tfinal) || (sfcontrol.operatingPoint() < 1) );

    return 0;
}
