#include <array>
#include <boost/numeric/odeint.hpp>
#include "rose499/sfcontrol.hpp"
#include "test/catch.hpp"

using namespace Eigen;

SCENARIO( "serret frenet controller applied to polynomial path", "[serretfrenet]" )
{
    namespace odeint = boost::numeric::odeint;
    constexpr double Tfinal = 40.0;
    double T = 0.0;
    double dt = 0.001;

    // Initial Condition
    auto ode45 = odeint::make_controlled( 1.0e-12 , 1.0e-12 , odeint::runge_kutta_fehlberg78<DriveSystem::StateType>() );

    // Generate a spline
    Matrix<Spline::ValueType, 2, 2> waypoints;
    waypoints(0, 0) = 1; waypoints(1, 0) = 0;
    waypoints(0, 1) = 100; waypoints(1, 1) = 100;

    // Hopefully the splining works haha
    Spline spline(waypoints);

    GIVEN( "constructed serret frenet controller" )
    {
        DriveSystem robot;
        SerretFrenetController sfcontrol(robot);
        sfcontrol.path(spline);

        WHEN( "simulated with initial condition on the path" )
        {
            DriveSystem::StateType x{ 1, 0, 0 };
            do
            {
                // If we fail
                if(odeint::fail == ode45.try_step(std::reference_wrapper<DriveController>(sfcontrol), x, T, dt))
                    continue;

                dt = std::min( dt, 0.05 );

                // Check if we converged to the path tangentially!
                if( sfcontrol.linearizedState().norm() == Approx(0).margin(1e-6) )
                    break;

            } while( T < Tfinal );

            THEN( "we should converge to the path and stay on it for all time" )
            {
                REQUIRE( sfcontrol.linearizedState().norm() == Approx(0).margin(1e-6) );
                while( T < Tfinal )
                {
                    // If we fail
                    if(odeint::fail == ode45.try_step(std::reference_wrapper<DriveController>(sfcontrol), x, T, dt))
                        continue;

                    dt = std::min( dt, 0.05 );

                    REQUIRE( sfcontrol.linearizedState().norm() == Approx(0).margin(1e-6) );
                }
            }
        }
    }
}
