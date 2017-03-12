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
    DriveSystem::StateType x{ 0, 0, 0 };
    auto ode45 = odeint::make_controlled( 1.0e-6 , 1.0e-6 , odeint::runge_kutta_fehlberg78<DriveSystem::StateType>() );

    // Generate a spline
    Matrix<Spline::ValueType, 2, 5> waypoints;
    waypoints(0, 0) = 1; waypoints(1, 0) = 0;
    waypoints(0, 1) = 4; waypoints(1, 1) = 5;
    waypoints(0, 2) = 8; waypoints(1, 2) = 7;
    waypoints(0, 3) = 12; waypoints(1, 3) = 4;
    waypoints(0, 4) = 16; waypoints(1, 4) = -1;

    Spline spline(waypoints);

    DriveSystem robot;
    SerretFrenetController sfcontrol(robot);
    sfcontrol.path(spline);
    do
    {
        // If we fail
        if(odeint::fail == ode45.try_step(std::reference_wrapper<DriveController>(sfcontrol), x, T, dt))
            continue;

        dt = std::min( dt, 0.05 );

        // If we pass store the system state to file.
        std::cout << T << "," << x[0] << ", " << x[1] << "," << x[2] << std::endl;
    } while( T < Tfinal );

    return 0;
}
