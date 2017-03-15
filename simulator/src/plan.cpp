#include <boost/geometry/geometries/register/point.hpp>
#include <Eigen/Core>
#include <list>
#include <map>
#include <random>
#include "rose499/world.hpp"

namespace geom = boost::geometry;

using Point = World::Point;

struct TaggedPoint
{
    double x, y;
    int tag;
};

BOOST_GEOMETRY_REGISTER_POINT_2D(TaggedPoint, double, geom::cs::cartesian, x, y)

std::list< Point > planRRT(     World::RTree knownWorld,
                                double xMax,
                                double yMax,
                                Point initialPoint,
                                double direction,
                                Point goal,
                                double radius )
{
    std::mt19937_64 engine;
    std::uniform_real_distribution<> samplerX(0, xMax), samplerY(0, yMax);

    TaggedPoint lastPoint = TaggedPoint{initialPoint.get<0>(), initialPoint.get<1>(), 0};

    geom::index::rtree< TaggedPoint, geom::index::linear<2, 1> > setPoints;
    std::map< uint32_t, uint32_t > backEdges;
    std::vector< TaggedPoint > listPoints;
    std::vector< TaggedPoint > setNearest;
    std::vector< World::Box > setCollisions;

    setPoints.insert( lastPoint );
    listPoints.push_back( lastPoint );
    int indTaggedPoint = 1;
    while( Eigen::Vector2d(lastPoint.x - goal.get<0>(), lastPoint.y - goal.get<1>()).norm() >= radius )
    {
      RetryRandomPoint:
        setNearest.clear();
        setCollisions.clear();

        auto randomPoint = TaggedPoint{ samplerX(engine), samplerY(engine), indTaggedPoint };

        geom::index::query(setPoints, geom::index::nearest(randomPoint, 1), std::back_inserter(setNearest));
        TaggedPoint nearestPoint = setNearest[0];

        // Need to do a collision check!
        auto start = World::Point(randomPoint.x, randomPoint.y);
        auto end = World::Point(nearestPoint.x, nearestPoint.y);
        geom::model::segment< World::Point > lineTo(start,end);
        World::Box bounds(  World::Point(std::min(randomPoint.x, nearestPoint.x), std::min(randomPoint.y, nearestPoint.y)),
                            World::Point(std::max(randomPoint.x, nearestPoint.x), std::max(randomPoint.y, nearestPoint.y)) );

        geom::index::query( knownWorld,
                            geom::index::intersects(bounds)
                                && !geom::index::contains(start)
                                && !geom::index::contains(end),
                            std::back_inserter(setCollisions));
        for(auto&& box : setCollisions)
        {
            // We actually know that the box doesn't WHOLLY contain our line,
            // so we can make some easing assumptions.
            auto dirVector = bounds.max_corner() - bounds.min_corner();
            if(lineTo, box))
                goto RetryRandomPoint;
        }

        setPoints.insert(randomPoint);
        listPoints.push_back(randomPoint);
        backEdges[randomPoint.tag] = nearestPoint.tag;
        indTaggedPoint++;
    }

    // Build path from back edges
    std::list< Point > path;
    while( lastPoint.tag != 0 )
    {
        path.push_front(Point(lastPoint.x, lastPoint.y));
      if( backEdges.find(lastPoint.tag) == backEdges.end() ) { break; }
        lastPoint = listPoints[backEdges[lastPoint.tag]];
    }
    path.push_front(Point(lastPoint.x, lastPoint.y));
    return path;
}
