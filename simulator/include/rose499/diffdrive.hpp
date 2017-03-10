#ifndef ROSE499_DIFFDRIVE_HPP
#define ROSE499_DIFFDRIVE_HPP

#include <array>
#include <Eigen/Core>
#include "rose499/types.hpp"

struct DriveSystem
{
    typedef SimulatorTypes::ValueType ValueType;
    typedef std::array<ValueType, 3> StateType;

    void step(StateType x, StateType& dxdt, double t, ValueType turnControl, ValueType speedControl);

private:
    Eigen::Matrix<ValueType, 3, 1> mState;
    Eigen::Matrix<ValueType, 3, 1> mFlow;
};

struct DriveController
{
    typedef DriveSystem::ValueType ValueType;
    typedef DriveSystem::StateType StateType;

    DriveController(DriveSystem &);

    void operator() (StateType x, StateType& dxdt, double t);

protected:
    virtual ValueType genSpeedControl(StateType x, double t);
    virtual ValueType genTurnControl(StateType x, double t);

private:
    DriveSystem& mSystem;
};

#endif
