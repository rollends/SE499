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
                    REQUIRE((splinePoint - waypoint).isZero() );
                }
            }

            THEN( "spline should observe C1 continuity on all intermediate points" )
            {
                for(auto i = 1; i < waypoints.cols() - 1; ++i)
                {
                    auto derivativeFromLeft = spline((i * 1.0) / SegmentCount, 1);
                    auto derivativeFromRight = spline((i * 1.0) / SegmentCount, 1);
                    REQUIRE( (derivativeFromLeft - derivativeFromRight).isZero() );
                }
            }

            THEN( "spline should observe C2 continuity on all intermediate points" )
            {
                for(auto i = 1; i < waypoints.cols() - 1; ++i)
                {
                    auto derivativeFromLeft = spline((i * 1.0) / SegmentCount, 2);
                    auto derivativeFromRight = spline((i * 1.0) / SegmentCount, 2);
                    REQUIRE( (derivativeFromLeft - derivativeFromRight).isZero() );
                }
            }
        }
    }
}
