#ifndef ROSE499_DIFFDRIVE_HPP
#define ROSE499_DIFFDRIVE_HPP

#include <array>
#include <Eigen/Core>
#include "rose499/types.hpp"
#include "rose499/spline.hpp"
#include "rose499/world.hpp"

struct DriveSystem
{
    typedef SimulatorTypes::ValueType ValueType;
    typedef std::array<ValueType, 3> StateType;

    void step(StateType x, StateType& dxdt, double t, ValueType turnControl, ValueType speedControl);

    double time() const;
    Eigen::Matrix<ValueType, 3, 1> const & state() const;
    Eigen::Matrix<ValueType, 3, 1> const & flow() const;

private:
    double mTime;
    Eigen::Matrix<ValueType, 3, 1> mState;
    Eigen::Matrix<ValueType, 3, 1> mFlow;
};

struct DriveController
{
    typedef DriveSystem::ValueType ValueType;
    typedef DriveSystem::StateType StateType;

    DriveController(DriveSystem &);

    void operator() (StateType x, StateType& dxdt, double t);

    void path( Spline );
    Spline const & path() const;

    std::ostream& printData(std::ostream&) const;
    std::ostream& printHeaders(std::ostream&) const;

    void updateKnownWorld(std::vector<World::Box> const & obstacles);

protected:
    virtual ValueType genSpeedControl(StateType x, double t);
    virtual ValueType genTurnControl(StateType x, double t);

    virtual std::ostream& printSpecificHeaders(std::ostream&) const;
    virtual std::ostream& printSpecificData(std::ostream&) const;

private:
    Spline mPath;
    DriveSystem& mSystem;
    World::RTree mKnownWorld;
};

namespace std
{
    std::ostream& operator << (std::ostream&, DriveController const &);
}

#endif
