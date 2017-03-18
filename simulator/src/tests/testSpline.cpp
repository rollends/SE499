#include <Eigen/Core>
#include "test/catch.hpp"
#include "rose499/spline.hpp"

using namespace Eigen;

SCENARIO( "waypoints can be splined with a 5th order C^2 polynomial", "[spline]" )
{
    GIVEN( "trivial linear sequence of waypoints" )
    {
        constexpr auto SegmentCount = 1;
        Matrix<Spline::ValueType, 2, 2> waypoints;
        waypoints(0, 0) = 0; waypoints(1, 0) = 0;
        waypoints(0, 1) = 1; waypoints(1, 1) = 1;

        WHEN( "splined via contructor" )
        {
            Spline spline(waypoints);
            THEN( "waypoints should be on the spline" )
            {
                for(auto i = 0; i < waypoints.cols(); ++i)
                {
                    auto splinePoint = spline((i * 1.0) / SegmentCount);
                    auto waypoint = waypoints.col(i);
                    REQUIRE((splinePoint - waypoint).isZero() );
                }
            }
        }

        WHEN( "splined via contructor with a direction constraint" )
        {
            Spline spline(waypoints, std::atan2(0, 1));
            THEN( "waypoints should be on the spline" )
            {
                for(auto i = 0; i < waypoints.cols(); ++i)
                {
                    auto splinePoint = spline((i * 1.0) / SegmentCount);
                    auto waypoint = waypoints.col(i);
                    REQUIRE((splinePoint - waypoint).isZero() );
                }
            }

            THEN( "derivative of spline at initial point is what was requested" )
            {
                auto slope = spline(0, 1);
                REQUIRE( std::atan2(slope[1], slope[0]) == Approx(std::atan2(0, 1)) );
            }
        }
    }

    GIVEN( "non-trivial sequence of waypoints" )
    {
        constexpr auto SegmentCount = 5 - 1;
        Matrix<Spline::ValueType, 2, 5> waypoints;
        waypoints(0, 0) = 0; waypoints(1, 0) = 0;
        waypoints(0, 1) = 2; waypoints(1, 1) = 1;
        waypoints(0, 2) = 4; waypoints(1, 2) = 0.5;
        waypoints(0, 3) = 6; waypoints(1, 3) = -0.5;
        waypoints(0, 4) = 3; waypoints(1, 4) = -1;

        WHEN( "splined via contructor" )
        {
            Spline spline(waypoints);
            THEN( "waypoints should be on the spline" )
            {
                for(auto i = 0; i < waypoints.cols(); ++i)
                {
                    auto splinePoint = spline((i * 1.0) / SegmentCount);
                    auto waypoint = waypoints.col(i);
                    CAPTURE( spline.poly() );
                    REQUIRE((splinePoint - waypoint).isZero() );
                }
            }

            THEN( "spline should observe C1 continuity on all intermediate points" )
            {
                for(auto i = 1; i < waypoints.cols() - 1; ++i)
                {
                    auto derivativeFromLeft = spline((i * 1.0) / SegmentCount - 1e-12, 1);
                    auto derivativeFromRight = spline((i * 1.0) / SegmentCount + 1e-12, 1);
                    CAPTURE( spline.dpoly() );
                    CAPTURE( i );
                    REQUIRE( (derivativeFromLeft - derivativeFromRight).norm() == Approx(0.0).margin(1e-12) );
                }
            }

            THEN( "spline should observe C2 continuity on all intermediate points" )
            {
                for(auto i = 1; i < waypoints.cols() - 1; ++i)
                {
                    auto derivativeFromLeft = spline((i * 1.0) / SegmentCount - 1e-12, 2);
                    auto derivativeFromRight = spline((i * 1.0) / SegmentCount + 1e-12, 2);
                    CAPTURE( spline.ddpoly() );
                    CAPTURE( i );
                    REQUIRE( (derivativeFromLeft - derivativeFromRight).norm() == Approx(0.0).margin(1e-12) );
                }
            }
        }

        GIVEN( "a working spline" )
        {
            Spline spline(waypoints);

            WHEN( "finding the nearest point on the curve to the end point" )
            {
                Spline::ValueType lambda = spline.nearestPoint(waypoints.col(waypoints.cols()-1), 0.9);

                THEN( "it ought to be 1")
                {
                    REQUIRE( lambda == Approx(1.0).margin(1e-12) );
                }
            }

            WHEN( "finding the nearest point on the curve to some non-trivial point" )
            {
                auto point = Matrix<Spline::ValueType, 2, 1>::Constant(1);
                Spline::ValueType lambda = spline.nearestPoint(point, 0.1);

                THEN( "the tangent of the closest point on the curve should be orthogonal to the error" )
                {
                    CAPTURE( lambda );
                    REQUIRE( spline(lambda, 1).dot(spline(lambda) - point) == Approx(0.0).margin(1e-12) );
                }
            }
        }
    }
}
