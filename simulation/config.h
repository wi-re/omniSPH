#pragma once
#include <atomic>

constexpr bool inletSwitch = false;
constexpr bool outletSwitch = false;
constexpr bool gravitySwitch = true;
constexpr bool backgroundPressureSwitch = false;
enum struct particleConfig {
    Domain, DamBreak, None
};
enum struct boundaryConfig {
    Box, ObstacleBT, ObstacleBTLow, Trapezoid, CenterBox
};
constexpr boundaryConfig bc = boundaryConfig::ObstacleBT;
constexpr particleConfig pc = particleConfig::DamBreak;

