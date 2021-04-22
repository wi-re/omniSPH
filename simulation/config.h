#pragma once
#include <atomic>

constexpr bool inletSwitch = true;
constexpr bool outletSwitch = true;
constexpr bool gravitySwitch = false;
constexpr bool backgroundPressureSwitch = false;
enum struct particleConfig {
    Domain, DamBreak, None
};
enum struct boundaryConfig {
    Box, ObstacleBT, ObstacleBTLow, Trapezoid, CenterBox
};
constexpr boundaryConfig bc = boundaryConfig::Trapezoid;
constexpr particleConfig pc = particleConfig::Domain;

