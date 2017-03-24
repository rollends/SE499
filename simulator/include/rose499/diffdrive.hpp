#ifndef ROSE499_DIFFDRIVE_HPP
#define ROSE499_DIFFDRIVE_HPP

#include <array>
#include <Eigen/Core>
#include <set>
#include "rose499/types.hpp"
#include "rose499/spline.hpp"
#include "rose499/world.hpp"

struct InvalidPathException{ };

struct DriveSystem
{
    typedef SimulatorTypes::ValueType ValueType;
    typedef std::array<ValueType, 3> StateType;

    void step(StateType x, StateType& dxdt, double t, ValueType turnControl, ValueType speedControl);

    double time() const;
    Eigen::Matrix<ValueType, 3, 1> const & state() const;
    Eigen::Matrix<ValueType, 3, 1> const & flow() const;
    void state(Eigen::Matrix<ValueType, 3, 1> x);
    void flow(Eigen::Matrix<ValueType, 3, 1> f);
    void time(double t);

    World::Polygon viewCone() const;

private:
    double mTime;
    Eigen::Matrix<ValueType, 3, 1> mState;
    Eigen::Matrix<ValueType, 3, 1> mFlow;
};

struct DriveController
{
    typedef DriveSystem::ValueType ValueType;
    typedef DriveSystem::StateType StateType;

    DriveController(DriveSystem &, Eigen::Matrix<ValueType, 2, 1> goal, ValueType goalRadius);

    void operator() (StateType x, StateType& dxdt, double t);

    void operatingPoint( Spline::ValueType );
    Spline::ValueType operatingPoint() const;
    void path( Spline const & );
    Spline const & path() const;
    Eigen::Matrix<ValueType, 2, 1> const & goal() const;

    std::ostream& printData(std::ostream&) const;
    std::ostream& printHeaders(std::ostream&) const;

    bool updateKnownWorld(std::list<std::pair<World::Box, int>> const & obstacles, bool forceReplan = false);
    void replan();

    virtual bool hasDiverged() const;

protected:
    virtual ValueType genSpeedControl(StateType x, double t);
    virtual ValueType genTurnControl(StateType x, double t);

    virtual std::ostream& printSpecificHeaders(std::ostream&) const;
    virtual std::ostream& printSpecificData(std::ostream&) const;

private:
    int mPlanIndex;
    Spline mPath;
    Spline::ValueType mOperatingLambda;
    DriveSystem& mSystem;
    World::RTree mKnownWorld;
    std::set<int> mKnownBoxes;
    Eigen::Matrix<ValueType, 2, 1> mGoal;
    ValueType mGoalRadius;
};

namespace std
{
    std::ostream& operator << (std::ostream&, DriveController const &);
}

#endif
