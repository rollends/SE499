#include <array>
#include <boost/numeric/odeint.hpp>

#include "rose499/diffdrive.hpp"

int main(int argc, char* argv[])
{
    namespace odeint = boost::numeric::odeint;

    constexpr double Tfinal = 15.0;
    double T = 0.0;
    double dt = 0.05;

    // Initial Condition
    DriveSystem::StateType x{ 10, 0 };
    odeint::controlled_runge_kutta< odeint::runge_kutta_dopri5<DriveSystem::StateType> > ode45;

    DriveSystem robot;
    DriveController defaultControl(robot);
    do
    {
        // If we fail
        if(ode45.try_step<DriveController&>(defaultControl, x, T, dt))
            continue;

        // If we pass store the system state to file.
        std::cout << x[0] << std::endl;
    } while( T < Tfinal );

    return 0;
}