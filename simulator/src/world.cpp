#include "rose499/world.hpp"

using Point = geom::model::point<double, 2, geom::cs::cartesian>;
using Box = geom::model::box<Point>;
using Polygon = geom::model::polygon<Point>;

constexpr double World::XMax;
constexpr double World::YMax;

void World::addBoxObstacle(int id, double x, double y, double w, double h)
{
    Point p1(x, y);
    Point p2(x + w, y + h);

    Box obstacle(p1, p2);
    mMap.insert(std::make_pair(obstacle, id));
}

std::list<std::pair<Box, int>> World::observeWorld(Polygon viewregion) const
{
    std::list< std::pair<Box, int> > resultSet;
    geom::index::query(mMap, geom::index::intersects(viewregion), std::back_inserter(resultSet));
    return resultSet;
}

World::RTree const & World::sceneTree() const { return mMap; }
