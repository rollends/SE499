#ifndef ROSE499_WORLD_HPP
#define ROSE499_WORLD_HPP

#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <boost/geometry/geometries/geometries.hpp>
#include <vector>

namespace geom = boost::geometry;

struct World
{
    using Point = geom::model::point<double, 2, geom::cs::cartesian>;
    using Polygon = geom::model::polygon<Point>;
    using Box = geom::model::box<Point>;
    using RTree = geom::index::rtree< Box, geom::index::linear<2, 1> >;

    /** Construct empty world **/
    World() = default;

    void addBoxObstacle(double x, double y, double w, double h);
    std::vector<Box> observeWorld(Polygon viewregion) const;

private:
    RTree mMap;
};

#endif
