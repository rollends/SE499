#ifndef ROSE499_DIFFDRIVE_HPP
#define ROSE499_DIFFDRIVE_HPP

#include <array>

#include "rose499/types.hpp"

struct DriveController;

struct DriveSystem
{
    typedef SimulatorTypes::ValueType ValueType;
    typedef std::array<ValueType, 3> StateType;

    DriveSystem(DriveController&);
    DriveSystem(DriveSystem const &);
    DriveSystem(DriveSystem&&) = delete;

    void operator() (StateType x, StateType& dxdt, double t);

  private:
    DriveController& control;
};

struct DriveController
{
    typedef DriveSystem::ValueType ValueType;
    typedef DriveSystem::StateType StateType;

    virtual ValueType genSpeedControl(StateType x, double t);
    virtual ValueType genTurnControl(StateType x, double t);
};

#endif
