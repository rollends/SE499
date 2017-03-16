#include <Eigen/Core>
#include <iterator>
#include <set>
#include <vector>

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

    GIVEN( "a world with many boxes" )
    {
        World world;
        for(int i = 0; i < 10; i++)
            world.addBoxObstacle(1 + 12 * i, 5, 5, 5);
        world.addBoxObstacle(0, 20, 95, 10);

        WHEN( "planning a path from [0,0], to [95,95]" )
        {
            auto path = planRRT( world.sceneTree(),     // Scene Tree (R-Tree)
                                 100, 100,              // Bounds of the Scene
                                 World::Point(0, 0),    // Initial condition
                                 0,                     // UNUSED: Direction
                                 World::Point(95, 95),  // Goal Point
                                 5                      // Radius tolerance for achieving goal.
                               );

            THEN( "there exists a non-empty path to the goal" )
            {
                REQUIRE( path.size() > 1 );

                THEN( "the path starts at the initial point and ends in the goal region" )
                {
                    auto start = path.front();
                    auto end = path.back();
                    REQUIRE( Eigen::Vector2d( start.get<0>(), start.get<1>() ).norm() < 1e-6 );
                    REQUIRE( Eigen::Vector2d( end.get<0>() - 95, end.get<1>() - 95 ).norm() <= 5 );
                }

                THEN( "the (segmented) path doesn't intersect any of the boxes in our world" )
                {
                    std::vector< World::Box > setCollisions;
                    auto start = std::begin(path);
                    auto end = std::next(start);
                    do {
                        geom::model::segment< World::Point > segment(*start, *end);
                        world.sceneTree().query(
                            geom::index::intersects(segment),
                            std::back_inserter(setCollisions)
                        );
                        REQUIRE(setCollisions.empty());

                        start = end;
                        end = std::next(end);
                    } while( end != std::end(path) );
                }
            }
        }
    }
}
