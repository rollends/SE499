#include <Eigen/Core>
#include <list>
#include <map>
#include <random>
#include "rose499/world.hpp"

namespace geom = boost::geometry;

using Point = World::Point;
using TaggedPoint = std::pair< Point, int >;

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

    TaggedPoint lastPoint(initialPoint, 0);

    geom::index::rtree< TaggedPoint, geom::index::linear<2, 1> > setPoints;
    std::map< uint32_t, uint32_t > backEdges;
    std::vector< TaggedPoint > listPoints;
    std::vector< TaggedPoint > setNearest;
    std::vector< World::Box > setCollisions;

    setPoints.insert( lastPoint );
    listPoints.push_back( lastPoint );
    int indTaggedPoint = 1;
    while( Eigen::Vector2d( lastPoint.first.get<0>() - goal.get<0>(),
                            lastPoint.first.get<1>() - goal.get<1>() ).norm() >= radius )
    {
        setNearest.clear();
        setCollisions.clear();

        auto randomPoint = TaggedPoint(Point(samplerX(engine), samplerY(engine)), indTaggedPoint);
        Eigen::Vector2d eigRandomPoint(randomPoint.first.get<0>(), randomPoint.first.get<1>());

        // Just make sure this isn't a crap point that is inside an obstacle.
        knownWorld.query( geom::index::intersects(randomPoint.first), std::back_inserter(setCollisions));
        if( !setCollisions.empty() )
            continue;

        // find nearest point in our set
        setPoints.query( geom::index::nearest(randomPoint.first, 10), std::back_inserter(setNearest) );
        std::sort(  setNearest.begin(),
                    setNearest.end(),
                    [eigRandomPoint]( TaggedPoint a, TaggedPoint b){
                        return  (Eigen::Vector2d(a.first.get<0>(), a.first.get<1>()) - eigRandomPoint).norm()
                                <
                                (Eigen::Vector2d(b.first.get<0>(), b.first.get<1>()) - eigRandomPoint).norm();
                    } );
        TaggedPoint nearestPoint = setNearest[0];

        // gotta make sure the path is obstacle free
        geom::model::segment< World::Point > lineTo(nearestPoint.first, randomPoint.first);
        knownWorld.query( geom::index::intersects(lineTo), std::back_inserter(setCollisions));
        if( !setCollisions.empty() )
            continue;

        setPoints.insert(randomPoint);
        listPoints.push_back(randomPoint);
        backEdges[randomPoint.second] = nearestPoint.second;
        indTaggedPoint++;
        lastPoint = randomPoint;
    }

    // Build path from back edges
    std::list< Point > path;
    int tag = 0;
    while( (tag = lastPoint.second) != 0 )
    {
        path.push_front(lastPoint.first);
        lastPoint = listPoints[backEdges[tag]];
    }
    path.push_front(lastPoint.first);
    return path;
}
