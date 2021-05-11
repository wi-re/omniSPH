#pragma once
#include <atomic>

inline bool inletSwitch = false;
constexpr bool outletSwitch = false;
constexpr bool gravitySwitch = true;
constexpr bool backgroundPressureSwitch = false;
enum struct particleConfig {
    Domain, DamBreak, None
};
enum struct boundaryConfig {
    Box, ObstacleBT, ObstacleBTLow, Trapezoid, CenterBox
};
constexpr boundaryConfig bc = boundaryConfig::ObstacleBTLow;
constexpr particleConfig pc = particleConfig::DamBreak;

