#pragma once
#include <atomic>

inline bool inletSwitch = false;
inline bool outletSwitch = false;
inline  bool gravitySwitch = false;
inline bool backgroundPressureSwitch = true;
enum struct particleConfig {
    Domain, DamBreak, None
};
enum struct boundaryConfig {
    Box, ObstacleBT, ObstacleBTLow, Trapezoid, CenterBox
};
constexpr boundaryConfig bc = boundaryConfig::ObstacleBTLow;
constexpr particleConfig pc = particleConfig::Domain;

enum struct scenario {
    dambreak, lid, sphere
};
constexpr scenario scenarioConfig = scenario::lid;