#pragma once
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#define _USE_MATH_DEFINES
#include <array>
#include <chrono>
#include <eigen3/Eigen/Dense>
#include <tools/ParameterManager.h>
#include <iostream>
#include <math.h>
#include <sstream>
#include <vector>
#include <glm/glm.hpp>

struct v4 { float x, y, z, w; };
GLuint create1DTexture(v4* colorMap, int32_t elements);
GLuint createProgram(std::string vertexSource, std::string fragmentSource);

void renderMarching();
void renderRays();

std::pair<bool, glm::vec2> intersectChebyshev(glm::vec2 origin, glm::vec2 direction, scalar hScaled, scalar gScaled);
std::pair<bool, glm::vec2> intersectImplicit(glm::vec2 origin, glm::vec2 direction, scalar hScaled, scalar gScaled);
std::pair<bool, glm::vec2> intersectExplicit(glm::vec2 origin, glm::vec2 direction, scalar hScaled, scalar gScaled);
std::tuple<bool, float, float> rayIntersectAABB(glm::vec2 orig, glm::vec2 dir, glm::vec2 b_min, glm::vec2 b_max);

inline auto square(scalar x) { return x * x; }
inline auto cubic(scalar x) { return x * x * x; }
inline auto renderKernel(scalar q) {
	return std::max(0.0, cubic(1.0 - q * q));
}
inline auto positionKernel(scalar q) {
	return std::max(0.0, 1.0 - cubic(q));
}
inline std::size_t cellsMX, cellsMY;
#include <simulation/2DMath.h>
inline auto gridPoint(vec gridPos, scalar h, scalar g) {
	vec xAvg(0.0, 0.0);
	scalar wSum = 0.0;
	scalar value = 0.0;


	auto [xi, yi] = getCellIdx(gridPos.x(), gridPos.y());
	auto range = (int32_t) ::ceil(h / scale);
	for (int32_t ix = -range; ix <= range; ++ix) {
		for (int32_t iy = -range; iy <= range; ++iy) {
			if (xi + ix < 0 || xi + ix >= cellsX || yi + iy < 0 || yi + iy >= cellsY)
				continue;
			const auto& cellData = getCell(xi + ix, yi + iy);
			for (const auto& j : cellData) {
				value += area * W(gridPos, particles[j].pos);
				float wi = renderKernel((gridPos - particles[j].pos).norm() / (h));
				xAvg += wi * particles[j].pos;
				wSum += wi;
			}
		}
	}
	//return -value;
	if (wSum > 0.0) {
		//return 0.0;
		xAvg /= wSum;
		return -(gridPos - xAvg).norm();
	}
	else
		return -1.0;
	//return cubic(h);


}
