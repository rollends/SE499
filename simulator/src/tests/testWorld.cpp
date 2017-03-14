#include "test/catch.hpp"
#include "rose499/world.hpp"

SCENARIO( "a world with boxes", "[world]" )
{
    GIVEN( "a world with a single box" )
    {
        World world;
        world.addBoxObstacle(1, 1, 10, 10);

        WHEN( "observing the world while facing the box" )
        {
            World::Polygon viewFrustum;
            geom::append(viewFrustum.outer(), World::Point(0.0, 0.0));
            geom::append(viewFrustum.outer(), World::Point(2.0, 0.0));
            geom::append(viewFrustum.outer(), World::Point(0.0, 2.0));
            geom::append(viewFrustum.outer(), World::Point(0.0, 0.0));
            auto obstacles = world.observeWorld(viewFrustum);

            THEN( "the box is observed" )
            {
                REQUIRE( !obstacles.empty() );
            }
        }

        WHEN( "observing the world while not facing the box" )
        {
            World::Polygon viewFrustum;
            geom::append(viewFrustum.outer(), World::Point(12.0, 12.0));
            geom::append(viewFrustum.outer(), World::Point(12 + 2.0, 12.0));
            geom::append(viewFrustum.outer(), World::Point(12.0, 12 + 2.0));
            geom::append(viewFrustum.outer(), World::Point(12.0, 12.0));
            auto obstacles = world.observeWorld(viewFrustum);

            THEN( "the box is not observed" )
            {
                REQUIRE( obstacles.empty() );
            }
        }
    }
}
