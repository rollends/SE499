#include <Eigen/Core>
#include <list>
#include <map>
#include <random>
#include "rose499/world.hpp"

namespace geom = boost::geometry;

using Point = World::Point;
using TaggedPoint = std::pair< Point, int >;

namespace
{
    std::mt19937_64 engine(0x1EA5ECED);
}

std::list< Point > planRRT(     World::RTree knownWorld,
                                double xMax,
                                double yMax,
                                Point initialPoint,
                                double direction,
                                Point goal,
                                double radius )
{
    std::uniform_real_distribution<> samplerX(0, xMax), samplerY(0, yMax);

    TaggedPoint lastPoint(initialPoint, 0);
    double closestDist = std::max(xMax, yMax);

    geom::index::rtree< TaggedPoint, geom::index::linear<5, 1> > setPoints;
    std::map< uint32_t, uint32_t > backEdges;
    std::vector< TaggedPoint > listPoints;
    std::vector< TaggedPoint > setNearest;
    std::vector< std::pair<World::Box, int> > setCollisions;

    setPoints.insert( lastPoint );
    listPoints.push_back( lastPoint );
    int indTaggedPoint = 1;
    while( Eigen::Vector2d( lastPoint.first.get<0>() - goal.get<0>(),
                            lastPoint.first.get<1>() - goal.get<1>() ).norm() >= radius )
    {
        setNearest.clear();
        setCollisions.clear();

        std::normal_distribution<> importanceSamplerX(goal.get<0>(), closestDist),
                                   importanceSamplerY(goal.get<1>(), closestDist);
        auto x = importanceSamplerX(engine);
        auto y = importanceSamplerY(engine);

      if( x > xMax || x < 0 || y > yMax || y < 0 ) { continue; }

        auto randomPoint = TaggedPoint(Point(x, y), indTaggedPoint);
        Eigen::Vector2d eigRandomPoint(randomPoint.first.get<0>(), randomPoint.first.get<1>());

        // Just make sure this isn't a crap point that is inside an obstacle.
        knownWorld.query( geom::index::intersects(randomPoint.first), std::back_inserter(setCollisions));

      if( !setCollisions.empty() ) { continue; }

        // find nearest point in our set
        setPoints.query( geom::index::nearest(randomPoint.first, 10), std::back_inserter(setNearest) );

        // Filter out points that are either too close OR do not fan out from the source point.
        auto endSetNearest =
            std::remove_if(
                setNearest.begin(),
                setNearest.end(),
                [direction, eigRandomPoint, &backEdges, &listPoints]( TaggedPoint nearest ) -> bool {
                    Eigen::Vector2d eigNearest(nearest.first.get<0>(), nearest.first.get<1>());
                    Eigen::Vector2d eigFrom(std::cos(direction), std::sin(direction));
                    Eigen::Vector2d vecLeaving = eigRandomPoint - eigNearest;
                    Eigen::Vector2d vecEntering = eigFrom;

                    if( nearest.second != 0 )
                    {
                        auto fromPoint = listPoints[backEdges[nearest.second]];
                        eigFrom = Eigen::Vector2d(fromPoint.first.get<0>(), fromPoint.first.get<1>());
                        vecEntering = eigNearest - eigFrom;
                    }

                    double proj = vecLeaving.dot(vecEntering) / (vecLeaving.norm() * vecEntering.norm());

                    return (vecLeaving.dot(vecEntering) < 0)
                        || (std::acos(proj) > std::atan2(1, 1))
                        ;//|| (vecLeaving.norm() < 1);
                }
            );
        setNearest.erase(endSetNearest, setNearest.end());

      if( setNearest.empty() ) { continue; }

        std::sort(  setNearest.begin(),
                    setNearest.end(),
                    [eigRandomPoint](TaggedPoint a, TaggedPoint b) -> bool {
                        return  (Eigen::Vector2d(a.first.get<0>(), a.first.get<1>()) - eigRandomPoint).squaredNorm()
                                <
                                (Eigen::Vector2d(b.first.get<0>(), b.first.get<1>()) - eigRandomPoint).squaredNorm();
                    } );
        TaggedPoint nearestPoint = setNearest.front();

        // gotta make sure the path is obstacle free
        geom::model::segment< World::Point > lineTo(nearestPoint.first, randomPoint.first);
        knownWorld.query( geom::index::intersects(lineTo), std::back_inserter(setCollisions));

      if( !setCollisions.empty() ) { continue; }

        closestDist = std::min( closestDist,
                        std::max( std::abs(nearestPoint.first.get<0>() - goal.get<0>()),
                                  std::abs(nearestPoint.first.get<1>() - goal.get<1>()) )
                      );

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
