#include <glad/glad.h>
#include <GLFW/glfw3.h>
// include order for gl has to be left like this
#include "SPH.h"
#include "SPHRender.h"
#include "Math.h"
#include "tinycolormap.h"
#include <chrono>
#include <iostream>
#include <sstream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtx/matrix_transform_2d.hpp>



int8_t sgn(float x) { return x > 0.f ? 1 : (x < 0.f ? -1 : 0); }

// Returns 1 if the lines intersect, otherwise 0. In addition, if the lines 
// intersect the intersection point may be stored in the floats i_x and i_y.
auto get_line_intersection(glm::vec2 p0, glm::vec2 p1, glm::vec2 p2, glm::vec2 p3)
{
	float s1_x, s1_y, s2_x, s2_y;
	s1_x = p1.x - p0.x;     s1_y = p1.y - p0.y;
	s2_x = p3.x - p2.x;     s2_y = p3.y - p2.y;

	float s, t;
	s = (-s1_y * (p0.x - p2.x) + s1_x * (p0.y - p2.y)) / (-s2_x * s1_y + s1_x * s2_y);
	t = (s2_x * (p0.y - p2.y) - s2_y * (p0.x - p2.x)) / (-s2_x * s1_y + s1_x * s2_y);

	if (s >= -1e-3 && s <= 1+1e-3 && t >= -1e-3 && t <= 1 + 1e-3)
	{
		return std::make_pair(true, glm::vec2(p0.x + t * s1_x, p0.y + t * s1_y));
	}

	return std::make_pair(false, glm::vec2(0,0)); // No collision
}


auto getCellIdxRender(float x, float y, float gScale) {
	//static auto& gScale = ParameterManager::instance().get<scalar>("ray.gScaleExplicit");
	static auto& hScale = ParameterManager::instance().get<scalar>("ray.hScale");

	if (x != x)
		x = (scalar)0.0;
	if (y != y)
		y = (scalar)0.0;
	std::size_t xi = static_cast<std::size_t>(::floor(std::clamp((double) x, (scalar)0.0, domainWidth) / gScale));
	std::size_t yi = static_cast<std::size_t>(::floor(std::clamp((double)y, (scalar)0.0, domainHeight) / gScale));
	xi = std::clamp(xi, (std::size_t)0, cellsMX);
	yi = std::clamp(yi, (std::size_t)0, cellsMY);

	//std::cout << x << " -> " << std::clamp((double)x, (scalar)0.0, domainWidth) << " -> " << std::clamp((double)x, (scalar)0.0, domainWidth) / gScale << " -> " << ::floor(std::clamp((double)x, (scalar)0.0, domainWidth) / gScale) << " --> " << static_cast<std::size_t>(::floor(std::clamp((double)x, (scalar)0.0, domainWidth) / gScale)) << " = " << xi << std::endl;;
	//std::cout << x << " -> " << std::clamp((double)y, (scalar)0.0, domainHeight) << " -> " << std::clamp((double)y, (scalar)0.0, domainHeight) / gScale << " -> " << ::floor(std::clamp((double)y, (scalar)0.0, domainHeight) / gScale) << " --> " << static_cast<std::size_t>(::floor(std::clamp((double)y, (scalar)0.0, domainHeight) / gScale)) << " = " << yi << std::endl;

	return std::make_pair((int32_t)xi, (int32_t)yi);
}

std::pair<bool, glm::vec2> intersectExplicit(glm::vec2 origin, glm::vec2 direction, scalar hScaled, scalar gScaled) {
	//std::cout << "Checking ray [ " << origin.x << " " << origin.y << " ] + t x [ " << direction.x << " " << direction.y << " ]\n";
	static auto& isoLevel = ParameterManager::instance().get<scalar>("ray.iso");
	static auto& hScale = ParameterManager::instance().get<scalar>("ray.hScale");

	// Explicit MC surfacing
	scalar h = support;

	cellsMX = (std::size_t) ::ceil((domainWidth) / gScaled);
	cellsMY = (std::size_t) ::ceil((domainHeight) / gScaled);

	auto [xi, yi] = getCellIdxRender(origin.x, origin.y, gScaled);
	glm::ivec2 voxel{ xi,yi };
	int32_t i = 0;
	auto getVoxelPosition = [&](auto voxel) {
		scalar xPos = 0.0 + (scalar)voxel.x * gScaled;
		scalar yPos = 0.0 + (scalar)voxel.y * gScaled;
		return glm::vec2(xPos, yPos);
	};
	auto getGrid = [&](auto x, auto y) {
		scalar xPos = 0.0 + (scalar)x * gScaled;
		scalar yPos = 0.0 + (scalar)y * gScaled;
		auto pos = glm::vec2(xPos, yPos);
		return std::make_pair(pos, gridPoint(vec(xPos, yPos), hScaled, gScaled));
	};
	auto vertexInterp = [](auto p1, auto p2, auto valp1, auto valp2) {
		float mu;
		glm::vec2 p;
		if (fabsf(isoLevel - valp1) < 0.00001f)
			return glm::vec2{ (float)p1.x, (float)p1.y };
		if (fabsf(isoLevel - valp2) < 0.00001f)
			return glm::vec2(p2.x, p2.y);
		if (fabsf(valp1 - valp2) < 0.00001f)
			return glm::vec2(p1.x, p1.y);
		mu = (isoLevel - valp1) / (valp2 - valp1);
		p.x = p1.x + mu * (p2.x - p1.x);
		p.y = p1.y + mu * (p2.y - p1.y);
		return p;
	};
	auto checkVoxel = [&](auto voxel) {
		return voxel.x > cellsMX || voxel.x < 0 || voxel.y > cellsMY || voxel.y < 0;
	};
	glLineWidth(0.5f);

	auto nxt = voxel + glm::ivec2{ direction.x > 0 ? 1 : 0, direction.y > 0 ? 1 : 0 };
	auto nxtB = getVoxelPosition(nxt);
	auto tMax = glm::vec2{
		fabsf(direction.x) <= 1e-5f ? FLT_MAX : (nxtB.x - origin.x) / direction.x,
		fabsf(direction.y) <= 1e-5f ? FLT_MAX : (nxtB.y - origin.y) / direction.y
	};
	do {
		if (checkVoxel(voxel))
			break;

		auto [posDL, sDL] = getGrid(voxel.x, voxel.y);
		auto [posDR, sDR] = getGrid(voxel.x + 1, voxel.y);
		auto [posTL, sTL] = getGrid(voxel.x, voxel.y + 1);
		auto [posTR, sTR] = getGrid(voxel.x + 1, voxel.y + 1);
		auto lerpL = vertexInterp(posDL, posTL, sDL, sTL);
		auto lerpR = vertexInterp(posDR, posTR, sDR, sTR);
		auto lerpD = vertexInterp(posDL, posDR, sDL, sDR);
		auto lerpT = vertexInterp(posTL, posTR, sTL, sTR);
		auto march = (sDL > isoLevel ? 0x1 : 0x0) | (sDR > isoLevel ? 0x2 : 0x0) | (sTR > isoLevel ? 0x4 : 0x0) | (sTL > isoLevel ? 0x8 : 0x0);
		const auto& DR = posDR;
		const auto& DL = posDL;
		const auto& TR = posTR;
		const auto& TL = posTL;
		auto intersection = [&](auto v0, auto v1) {
			return get_line_intersection(origin, origin + direction * 200.f, v0, v1);
		};
#define intersectLine(v0, v1) \
if(auto [lhit, pos] = intersection(v0,v1);lhit) { \
return std::make_pair(true, pos);\
}
		bool hit = false;
		switch (march) {
		case 0: break;
		case 1: {
			intersectLine(lerpD, lerpL);
		}break;
		case 2: {
			intersectLine(lerpR, lerpD);
		} break;
		case 3: {
			intersectLine(lerpR, lerpL);
		} break;
		case 4: {
			intersectLine(lerpT, lerpR);
		} break;
		case 5: {
			if (sDL + sDR + sTL + sTR > 4.0 * isoLevel) {
				intersectLine(lerpD, lerpR);
				intersectLine(lerpL, lerpT);
			}
			else {
				intersectLine(lerpL, lerpD);
				intersectLine(lerpT, lerpR);
			}
		} break;
		case 6: {
			intersectLine(lerpT, lerpD);
		} break;
		case 7: {
			intersectLine(lerpT, lerpL);
		} break;
		case 8: {
			intersectLine(lerpL, lerpT);
		} break;
		case 9: {
			intersectLine(lerpD, lerpT);
		} break;
		case 10: {
			if (sDL + sDR + sTL + sTR > 4.0 * isoLevel) {
				intersectLine(lerpD, lerpL);
				intersectLine(lerpR, lerpT);
			}
			else {
				intersectLine(lerpL, lerpT);
				intersectLine(lerpD, lerpR);
			}
		} break;
		case 11: {
			intersectLine(lerpR, lerpT);
		} break;
		case 12: {
			intersectLine(lerpL, lerpR);
		} break;
		case 13: {
			intersectLine(lerpD, lerpR);
		} break;
		case 14: {
			intersectLine(lerpL, lerpD);
		} break;
		case 15: {
		} break;
		}
		if (hit) break;

		if (tMax.x < tMax.y) {
			if (checkVoxel(voxel)) break;
			voxel.x += sgn(direction.x);
			tMax.x += sgn(direction.x) / direction.x * gScaled;
		}
		else {
			if (checkVoxel(voxel)) break;
			voxel.y += sgn(direction.y);
			tMax.y += sgn(direction.y) / direction.y * gScaled;
		}





	} while (true);
	//std::cout << "No intersection!" << std::endl;
	return std::make_pair(false, glm::vec2(0, 0));

}

void intersectExplicitVis(){
	static auto& orig = ParameterManager::instance().get<vec>("ray.origin");
	static auto& target = ParameterManager::instance().get<vec>("ray.target");
	static auto& fov = ParameterManager::instance().get<scalar>("ray.fov");
	static auto& res = ParameterManager::instance().get<int32_t>("ray.resolution");
	static auto& isoLevel = ParameterManager::instance().get<scalar>("ray.iso");
	static auto& gScale = ParameterManager::instance().get<scalar>("ray.gScaleExplicit");
	static auto& hScale = ParameterManager::instance().get<scalar>("ray.hScale");
	vec dir = target - orig;
	if (dir.norm() > 1e-5)
		dir.normalize();

	glm::vec2 origin(orig.x(), orig.y());
	glm::vec2 direction(dir.x(), dir.y());
	// Explicit MC surfacing
	scalar h = support;
	scalar hScaled = h * hScale;
	scalar gScaled = h * gScale;

	cellsMX = (std::size_t) ::ceil((domainWidth) / gScaled);
	cellsMY = (std::size_t) ::ceil((domainHeight) / gScaled);

	auto [xi, yi] = getCellIdxRender(origin.x, origin.y, gScaled);
	glm::ivec2 voxel{ xi,yi };
	int32_t i = 0;
	auto getVoxelPosition = [&](auto voxel) {
		scalar xPos = 0.0 + (scalar)voxel.x * gScaled;
		scalar yPos = 0.0 + (scalar)voxel.y * gScaled;
		return glm::vec2(xPos, yPos);
	};
	auto getGrid = [&](auto x, auto y) {
		scalar xPos = 0.0 + (scalar)x * gScaled;
		scalar yPos = 0.0 + (scalar)y * gScaled;
		auto pos = glm::vec2(xPos, yPos);
		return std::make_pair(pos, gridPoint(vec(xPos, yPos), hScaled, gScaled));
	};
	auto vertexInterp = [](auto p1, auto p2, auto valp1, auto valp2) {
		float mu;
		glm::vec2 p;
		if (fabsf(isoLevel - valp1) < 0.00001f)
			return glm::vec2{ (float)p1.x, (float)p1.y };
		if (fabsf(isoLevel - valp2) < 0.00001f)
			return glm::vec2(p2.x, p2.y);
		if (fabsf(valp1 - valp2) < 0.00001f)
			return glm::vec2(p1.x, p1.y);
		mu = (isoLevel - valp1) / (valp2 - valp1);
		p.x = p1.x + mu * (p2.x - p1.x);
		p.y = p1.y + mu * (p2.y - p1.y);
		return p;
	};
	auto checkVoxel = [&](auto voxel) {
		return voxel.x > cellsMX || voxel.x < 0 || voxel.y > cellsMY || voxel.y < 0;
	};
	glLineWidth(0.5f);

	auto nxt = voxel + glm::ivec2{ direction.x > 0 ? 1 : 0, direction.y > 0 ? 1 : 0 };
	auto nxtB = getVoxelPosition(nxt);
	auto tMax = glm::vec2{
		fabsf(direction.x) <= 1e-5f ? FLT_MAX : (nxtB.x - origin.x) / direction.x,
		fabsf(direction.y) <= 1e-5f ? FLT_MAX : (nxtB.y - origin.y) / direction.y
	};
	do {
		if (checkVoxel(voxel))
			break;

		auto [posDL, sDL] = getGrid(voxel.x, voxel.y);
		auto [posDR, sDR] = getGrid(voxel.x + 1, voxel.y);
		auto [posTL, sTL] = getGrid(voxel.x, voxel.y + 1);
		auto [posTR, sTR] = getGrid(voxel.x + 1, voxel.y + 1);
		auto lerpL = vertexInterp(posDL, posTL, sDL, sTL);
		auto lerpR = vertexInterp(posDR, posTR, sDR, sTR);
		auto lerpD = vertexInterp(posDL, posDR, sDL, sDR);
		auto lerpT = vertexInterp(posTL, posTR, sTL, sTR);
		auto march = (sDL > isoLevel ? 0x1 : 0x0) | (sDR > isoLevel ? 0x2 : 0x0) | (sTR > isoLevel ? 0x4 : 0x0) | (sTL > isoLevel ? 0x8 : 0x0);
		const auto& DR = posDR;
		const auto& DL = posDL;
		const auto& TR = posTR;
		const auto& TL = posTL;
		glBegin(GL_LINES);
		if (i++ == 0)
			glColor3f(1.f, 0.f, 0.f);
		else
			glColor3f(1.0f, 0.647f, 0.f);

		glVertex2f(posDL.x, posDL.y);
		glVertex2f(posDR.x, posDR.y);
		glVertex2f(posDR.x, posDR.y);
		glVertex2f(posTR.x, posTR.y);
		glVertex2f(posTR.x, posTR.y);
		glVertex2f(posTL.x, posTL.y);
		glVertex2f(posTL.x, posTL.y);
		glVertex2f(posDL.x, posDL.y);
		glEnd();

		auto intersection = [&](auto v0, auto v1) {
			return get_line_intersection(origin, origin + direction * 200.f, v0, v1);
		};
#define intersectLine(v0, v1) \
if(auto [lhit, pos] = intersection(v0,v1);lhit) { \
	glBegin(GL_POINTS);\
	glColor3f(0.4588f, 0.4902f, 0.4588f);\
	glVertex2f(pos.x, pos.y);\
	glEnd();\
	hit = true;\
	break;\
}
		bool hit = false;
		switch (march) {
		case 0: break;
		case 1: {
			intersectLine(lerpD, lerpL);
		}break;
		case 2: {
			intersectLine(lerpR, lerpD);
		} break;
		case 3: {
			intersectLine(lerpR, lerpL);
		} break;
		case 4: {
			intersectLine(lerpT, lerpR);
		} break;
		case 5: {
			if (sDL + sDR + sTL + sTR > 4.0 * isoLevel) {
				intersectLine(lerpD, lerpR);
				intersectLine(lerpL, lerpT);
			}
			else {
				intersectLine(lerpL, lerpD);
				intersectLine(lerpT, lerpR);
			}
		} break;
		case 6: {
			intersectLine(lerpT, lerpD);
		} break;
		case 7: {
			intersectLine(lerpT, lerpL);
		} break;
		case 8: {
			intersectLine(lerpL, lerpT);
		} break;
		case 9: {
			intersectLine(lerpD, lerpT);
		} break;
		case 10: {
			if (sDL + sDR + sTL + sTR > 4.0 * isoLevel) {
				intersectLine(lerpD, lerpL);
				intersectLine(lerpR, lerpT);
			}
			else {
				intersectLine(lerpL, lerpT);
				intersectLine(lerpD, lerpR);
			}
		} break;
		case 11: {
			intersectLine(lerpR, lerpT);
		} break;
		case 12: {
			intersectLine(lerpL, lerpR);
		} break;
		case 13: {
			intersectLine(lerpD, lerpR);
		} break;
		case 14: {
			intersectLine(lerpL, lerpD);
		} break;
		case 15: {
		} break;
		}
		if (hit) break;

		if (tMax.x < tMax.y) {
			if (checkVoxel(voxel)) break;
			voxel.x += sgn(direction.x);
			tMax.x += sgn(direction.x) / direction.x * gScale;
		}
		else {
			if (checkVoxel(voxel)) break;
			voxel.y += sgn(direction.y);
			tMax.y += sgn(direction.y) / direction.y * gScale;
		}





	} while (true);
}


std::tuple<bool, float, float> rayIntersectAABB(glm::vec2 orig, glm::vec2 dir, glm::vec2 b_min, glm::vec2 b_max) {
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))
	float t1 = (b_min.x - orig.x) / dir.x;
	float t2 = (b_max.x - orig.x) / dir.x;

	float tmin = MIN(t1, t2);
	float tmax = MAX(t1, t2);

	t1 = (b_min.y - orig.y) / dir.y;
	t2 = (b_max.y - orig.y) / dir.y;

	tmin = MAX(tmin, MIN(t1, t2));
	tmax = MIN(tmax, MAX(t1, t2));

	return std::make_tuple(tmax > MAX(tmin, 0.f), tmin, tmax);
}

void intersectImplicitVis() {
	static auto& subSteps = ParameterManager::instance().get<int32_t>("ray.subSteps");
	static auto& orig = ParameterManager::instance().get<vec>("ray.origin");
	static auto& target = ParameterManager::instance().get<vec>("ray.target");
	static auto& fov = ParameterManager::instance().get<scalar>("ray.fov");
	static auto& res = ParameterManager::instance().get<int32_t>("ray.resolution");
	static auto& isoLevel = ParameterManager::instance().get<scalar>("ray.iso");
	static auto& gScale = ParameterManager::instance().get<scalar>("ray.gScaleExplicit");
	static auto& hScale = ParameterManager::instance().get<scalar>("ray.hScale");
	vec dir = target - orig;
	if (dir.norm() > 1e-5)
		dir.normalize();

	glm::vec2 origin(orig.x(), orig.y());
	glm::vec2 direction(dir.x(), dir.y());
	// Explicit MC surfacing
	scalar h = support;
	scalar hScaled = h * hScale;
	scalar gScaled = h * gScale;

	cellsMX = (std::size_t) ::ceil((domainWidth) / gScaled);
	cellsMY = (std::size_t) ::ceil((domainHeight) / gScaled);

	auto [xi, yi] = getCellIdxRender(origin.x, origin.y, gScaled);
	glm::ivec2 voxel{ xi,yi };
	int32_t i = 0;
	auto getVoxelPosition = [&](auto voxel) {
		scalar xPos = 0.0 + (scalar)voxel.x * gScaled;
		scalar yPos = 0.0 + (scalar)voxel.y * gScaled;
		return glm::vec2(xPos, yPos);
	};
	auto getGrid = [&](auto x, auto y) {
		scalar xPos = 0.0 + (scalar)x * gScaled;
		scalar yPos = 0.0 + (scalar)y * gScaled;
		auto pos = glm::vec2(xPos, yPos);
		return std::make_pair(pos, gridPoint(vec(xPos, yPos), hScaled, gScaled));
	};
	auto vertexInterp = [](auto p1, auto p2, auto valp1, auto valp2) {
		float mu;
		glm::vec2 p;
		if (fabsf(isoLevel - valp1) < 0.00001f)
			return glm::vec2{ (float)p1.x, (float)p1.y };
		if (fabsf(isoLevel - valp2) < 0.00001f)
			return glm::vec2(p2.x, p2.y);
		if (fabsf(valp1 - valp2) < 0.00001f)
			return glm::vec2(p1.x, p1.y);
		mu = (isoLevel - valp1) / (valp2 - valp1);
		p.x = p1.x + mu * (p2.x - p1.x);
		p.y = p1.y + mu * (p2.y - p1.y);
		return p;
	};
	auto checkVoxel = [&](auto voxel) {
		return voxel.x > cellsMX || voxel.x < 0 || voxel.y > cellsMY || voxel.y < 0;
	};
	glLineWidth(0.5f);

	auto nxt = voxel + glm::ivec2{ direction.x > 0 ? 1 : 0, direction.y > 0 ? 1 : 0 };
	auto nxtB = getVoxelPosition(nxt);
	auto tMax = glm::vec2{
		fabsf(direction.x) <= 1e-5f ? FLT_MAX : (nxtB.x - origin.x) / direction.x,
		fabsf(direction.y) <= 1e-5f ? FLT_MAX : (nxtB.y - origin.y) / direction.y
	};
	auto inside = gridPoint(vec(origin.x, origin.y), hScaled, gScaled) < isoLevel ? 0 : 1;
	//std::cout << "Ray started " << (inside ? "inside" : "outside") << std::endl;
	//std::cout << "Starting traversal at voxel [ " << voxel.x << " " << voxel.y << " ]" << std::endl;
	do {
		//std::cout << "Traversal at voxel [ " << voxel.x << " " << voxel.y << " ]" << std::endl;
		if (checkVoxel(voxel))
			break;
		auto vmin = getVoxelPosition(voxel);
		auto vmax = getVoxelPosition(voxel + glm::ivec2(1, 1));

		//glBegin(GL_LINES);
		//if (i++ == 0)
		//	glColor3f(1.f, 0.f, 0.f);
		//else
		//	glColor3f(1.0f, 0.647f, 0.f);
		//glVertex2f(vmin.x, vmin.y);
		//glVertex2f(vmax.x, vmin.y);
		//glVertex2f(vmax.x, vmin.y);
		//glVertex2f(vmax.x, vmax.y);
		//glVertex2f(vmax.x, vmax.y);
		//glVertex2f(vmin.x, vmax.y);
		//glVertex2f(vmin.x, vmax.y);
		//glVertex2f(vmin.x, vmin.y);
		//glEnd();

		auto [hit, tmin, tmax] = rayIntersectAABB(origin, direction, vmin, vmax);
		tmin = std::max(0.f, tmin);
		std::vector<double> valArray;
		valArray.resize(subSteps);
		//std::cout << "Values: " << std::endl;
		for (int32_t i = 0; i < subSteps; ++i) {
			float t = tmin + (tmax - tmin) * (((float)(i)) / ((double)subSteps - 1));
			auto pos = origin + t * direction;
			valArray[i] = gridPoint(vec(pos.x, pos.y), hScaled, gScaled);
			//std::cout << "\t" << valArray[i] << "\n";
		}
		if (!inside) {
			for (int32_t i = 1; i < subSteps; ++i) {
				if (valArray[i] > isoLevel) {
					//std::cout << "Found ISO-Surface" << std::endl;
					auto y0 = valArray[i - 1];
					auto y1 = valArray[i];
					//std::cout << "Change from " << y0 << " to " << y1 << std::endl;
					auto dy = y1 - y0;
					auto alpha = (isoLevel- y0) / dy;
					alpha = std::clamp(alpha, 0.0, 1.0);
					//std::cout << "Alpha value: " << alpha << std::endl;
					auto t0 = tmin + (tmax - tmin) * (((float)(i - 1)) / ((double)subSteps - 1));
					auto t1 = tmin + (tmax - tmin) * (((float)(i)) / ((double)subSteps - 1));;
					//std::cout << "t0: " << t0 << "\tt1:" << t1 << std::endl;
					float t = t0 * (1.f - alpha) + t1 * alpha;
					//std::cout << "Lerp t: " << t << std::endl;
					auto pos = origin + t * direction;
					//std::cout << "Final Position: " << pos.x << " x " << pos.y << std::endl;
					glBegin(GL_POINTS); 
					glColor3f(0.2667f, 0.4235f, 0.8118f);  // han blue
					glVertex2f(pos.x, pos.y); 
					glEnd(); 
					return;
				}
			}
		}else{
			for (int32_t i = 1; i < subSteps; ++i) {
				if (valArray[i] < isoLevel) {
					//std::cout << "Found ISO-Surface" << std::endl;
					auto y0 = valArray[i - 1];
					auto y1 = valArray[i];
					//std::cout << "Change from " << y0 << " to " << y1 << std::endl;
					auto dy = y1 - y0;
					auto alpha = (isoLevel - y0) / dy;
					alpha = std::clamp(alpha, 0.0, 1.0);
					//std::cout << "Alpha value: " << alpha << std::endl;
					//alpha = 0.f;
					auto t0 = tmin + (tmax - tmin) * (((float)(i - 1)) / ((double)subSteps - 1));
					auto t1 = tmin + (tmax - tmin) * (((float)(i)) / ((double)subSteps - 1));;
					//std::cout << "t0: " << t0 << "\tt1:" << t1 << std::endl;
					float t = t0 * (1.f - alpha) + t1 * alpha;
					//std::cout << "Lerp t: " << t << std::endl;
					auto pos = origin + t * direction;
					//std::cout << "Final Position: " << pos.x << " x " << pos.y << std::endl;
					glBegin(GL_POINTS);
					glColor3f(0.2667f, 0.4235f, 0.8118f);
					glVertex2f(pos.x, pos.y);
					glEnd();
					return;
				}
			}
		}
	

		//if (hit) break;
		//std::cout << "Moving to next voxel" << std::endl;
		if (tMax.x < tMax.y) {
			if (checkVoxel(voxel)) break;
			//std::cout << "Increment in x direction" << std::endl;
			voxel.x += sgn(direction.x);
			tMax.x += sgn(direction.x) / direction.x * gScale;
		}
		else {
			if (checkVoxel(voxel)) break;
			//std::cout << "Increment in x direction" << std::endl;
			voxel.y += sgn(direction.y);
			tMax.y += sgn(direction.y) / direction.y * gScale;
		}
		//std::cout << "Traversal step is done" << std::endl;
	} while (true);
}

std::pair<bool, glm::vec2> intersectImplicit(glm::vec2 origin, glm::vec2 direction, scalar hScaled, scalar gScaled) {
	static auto& subSteps = ParameterManager::instance().get<int32_t>("ray.subSteps");
	static auto& orig = ParameterManager::instance().get<vec>("ray.origin");
	static auto& target = ParameterManager::instance().get<vec>("ray.target");
	static auto& fov = ParameterManager::instance().get<scalar>("ray.fov");
	static auto& res = ParameterManager::instance().get<int32_t>("ray.resolution");
	static auto& isoLevel = ParameterManager::instance().get<scalar>("ray.iso");
	static auto& gScale = ParameterManager::instance().get<scalar>("ray.gScaleExplicit");
	static auto& hScale = ParameterManager::instance().get<scalar>("ray.hScale");

	// Explicit MC surfacing
	scalar h = support;

	cellsMX = (std::size_t) ::ceil((domainWidth) / gScaled);
	cellsMY = (std::size_t) ::ceil((domainHeight) / gScaled);

	auto [xi, yi] = getCellIdxRender(origin.x, origin.y, gScaled);
	glm::ivec2 voxel{ xi,yi };
	int32_t i = 0;
	auto getVoxelPosition = [&](auto voxel) {
		scalar xPos = 0.0 + (scalar)voxel.x * gScaled;
		scalar yPos = 0.0 + (scalar)voxel.y * gScaled;
		return glm::vec2(xPos, yPos);
	};
	auto getGrid = [&](auto x, auto y) {
		scalar xPos = 0.0 + (scalar)x * gScaled;
		scalar yPos = 0.0 + (scalar)y * gScaled;
		auto pos = glm::vec2(xPos, yPos);
		return std::make_pair(pos, gridPoint(vec(xPos, yPos), hScaled, gScaled));
	};
	auto vertexInterp = [](auto p1, auto p2, auto valp1, auto valp2) {
		float mu;
		glm::vec2 p;
		if (fabsf(isoLevel - valp1) < 0.00001f)
			return glm::vec2{ (float)p1.x, (float)p1.y };
		if (fabsf(isoLevel - valp2) < 0.00001f)
			return glm::vec2(p2.x, p2.y);
		if (fabsf(valp1 - valp2) < 0.00001f)
			return glm::vec2(p1.x, p1.y);
		mu = (isoLevel - valp1) / (valp2 - valp1);
		p.x = p1.x + mu * (p2.x - p1.x);
		p.y = p1.y + mu * (p2.y - p1.y);
		return p;
	};
	auto checkVoxel = [&](auto voxel) {
		return voxel.x > cellsMX || voxel.x < 0 || voxel.y > cellsMY || voxel.y < 0;
	};
	glLineWidth(0.5f);

	auto nxt = voxel + glm::ivec2{ direction.x > 0 ? 1 : 0, direction.y > 0 ? 1 : 0 };
	auto nxtB = getVoxelPosition(nxt);
	auto tMax = glm::vec2{
		fabsf(direction.x) <= 1e-5f ? FLT_MAX : (nxtB.x - origin.x) / direction.x,
		fabsf(direction.y) <= 1e-5f ? FLT_MAX : (nxtB.y - origin.y) / direction.y
	};
	auto inside = gridPoint(vec(origin.x, origin.y), hScaled, gScaled) < isoLevel ? 0 : 1;
	//std::cout << "Ray started " << (inside ? "inside" : "outside") << std::endl;
	//std::cout << "Starting traversal at voxel [ " << voxel.x << " " << voxel.y << " ]" << std::endl;
	do {
		//std::cout << "Traversal at voxel [ " << voxel.x << " " << voxel.y << " ]" << std::endl;
		if (checkVoxel(voxel))
			break;
		auto vmin = getVoxelPosition(voxel);
		auto vmax = getVoxelPosition(voxel + glm::ivec2(1, 1));

		auto [hit, tmin, tmax] = rayIntersectAABB(origin, direction, vmin, vmax);
		tmin = std::max(0.f, tmin);
		std::vector<double> valArray;
		valArray.resize(subSteps);
		//std::cout << "Values: " << std::endl;
		for (int32_t i = 0; i < subSteps; ++i) {
			float t = tmin + (tmax - tmin) * (((float)(i)) / ((double)subSteps - 1));
			auto pos = origin + t * direction;
			valArray[i] = gridPoint(vec(pos.x, pos.y), hScaled, gScaled);
			//std::cout << "\t" << valArray[i] << "\n";
		}
		if (!inside) {
			for (int32_t i = 1; i < subSteps; ++i) {
				if (valArray[i] > isoLevel) {
					//std::cout << "Found ISO-Surface" << std::endl;
					auto y0 = valArray[i - 1];
					auto y1 = valArray[i];
					//std::cout << "Change from " << y0 << " to " << y1 << std::endl;
					auto dy = y1 - y0;
					auto alpha = (isoLevel - y0) / dy;
					alpha = std::clamp(alpha, 0.0, 1.0);
					//std::cout << "Alpha value: " << alpha << std::endl;
					auto t0 = tmin + (tmax - tmin) * (((float)(i - 1)) / ((double)subSteps - 1));
					auto t1 = tmin + (tmax - tmin) * (((float)(i)) / ((double)subSteps - 1));;
					//std::cout << "t0: " << t0 << "\tt1:" << t1 << std::endl;
					float t = t0 * (1.f - alpha) + t1 * alpha;
					//std::cout << "Lerp t: " << t << std::endl;
					auto pos = origin + t * direction;
					//std::cout << "Final Position: " << pos.x << " x " << pos.y << std::endl;
					return std::make_pair(true, pos);
				}
			}
		}
		else {
			for (int32_t i = 1; i < subSteps; ++i) {
				if (valArray[i] < isoLevel) {
					//std::cout << "Found ISO-Surface" << std::endl;
					auto y0 = valArray[i - 1];
					auto y1 = valArray[i];
					//std::cout << "Change from " << y0 << " to " << y1 << std::endl;
					auto dy = y1 - y0;
					auto alpha = (isoLevel - y0) / dy;
					alpha = std::clamp(alpha, 0.0, 1.0);
					//std::cout << "Alpha value: " << alpha << std::endl;
					//alpha = 0.f;
					auto t0 = tmin + (tmax - tmin) * (((float)(i - 1)) / ((double)subSteps - 1));
					auto t1 = tmin + (tmax - tmin) * (((float)(i)) / ((double)subSteps - 1));;
					//std::cout << "t0: " << t0 << "\tt1:" << t1 << std::endl;
					float t = t0 * (1.f - alpha) + t1 * alpha;
					//std::cout << "Lerp t: " << t << std::endl;
					auto pos = origin + t * direction;
					//std::cout << "Final Position: " << pos.x << " x " << pos.y << std::endl;
					return std::make_pair(true, pos);
				}
			}
		}


		//if (hit) break;
		//std::cout << "Moving to next voxel" << std::endl;
		if (tMax.x < tMax.y) {
			if (checkVoxel(voxel)) break;
			//std::cout << "Increment in x direction" << std::endl;
			voxel.x += sgn(direction.x);
			tMax.x += sgn(direction.x) / direction.x * gScale;
		}
		else {
			if (checkVoxel(voxel)) break;
			//std::cout << "Increment in x direction" << std::endl;
			voxel.y += sgn(direction.y);
			tMax.y += sgn(direction.y) / direction.y * gScale;
		}
		//std::cout << "Traversal step is done" << std::endl;
	} while (true);
	return std::make_pair(false, glm::vec2(0,0));
}
#include <cheb/cheb.h>
#include <iomanip>
std::pair<bool, glm::vec2> intersectChebyshev(glm::vec2 origin, glm::vec2 direction, scalar hScaled, scalar gScaled) {
	static auto& subSteps = ParameterManager::instance().get<int32_t>("ray.subSteps");
	static auto& orig = ParameterManager::instance().get<vec>("ray.origin");
	static auto& target = ParameterManager::instance().get<vec>("ray.target");
	static auto& fov = ParameterManager::instance().get<scalar>("ray.fov");
	static auto& res = ParameterManager::instance().get<int32_t>("ray.resolution");
	static auto& isoLevel = ParameterManager::instance().get<scalar>("ray.iso");
	static auto& gScale = ParameterManager::instance().get<scalar>("ray.gScaleExplicit");
	static auto& hScale = ParameterManager::instance().get<scalar>("ray.hScale");

	// Explicit MC surfacing
	scalar h = support;

	cellsMX = (std::size_t) ::ceil((domainWidth) / gScaled);
	cellsMY = (std::size_t) ::ceil((domainHeight) / gScaled);

	auto [xi, yi] = getCellIdxRender(origin.x, origin.y, gScaled);
	glm::ivec2 voxel{ xi,yi };
	int32_t i = 0;
	auto getVoxelPosition = [&](auto voxel) {
		scalar xPos = 0.0 + (scalar)voxel.x * gScaled;
		scalar yPos = 0.0 + (scalar)voxel.y * gScaled;
		return glm::vec2(xPos, yPos);
	};
	auto getGrid = [&](auto x, auto y) {
		scalar xPos = 0.0 + (scalar)x * gScaled;
		scalar yPos = 0.0 + (scalar)y * gScaled;
		auto pos = glm::vec2(xPos, yPos);
		return std::make_pair(pos, gridPoint(vec(xPos, yPos), hScaled, gScaled));
	};
	auto vertexInterp = [](auto p1, auto p2, auto valp1, auto valp2) {
		float mu;
		glm::vec2 p;
		if (fabsf(isoLevel - valp1) < 0.00001f)
			return glm::vec2{ (float)p1.x, (float)p1.y };
		if (fabsf(isoLevel - valp2) < 0.00001f)
			return glm::vec2(p2.x, p2.y);
		if (fabsf(valp1 - valp2) < 0.00001f)
			return glm::vec2(p1.x, p1.y);
		mu = (isoLevel - valp1) / (valp2 - valp1);
		p.x = p1.x + mu * (p2.x - p1.x);
		p.y = p1.y + mu * (p2.y - p1.y);
		return p;
	};
	auto checkVoxel = [&](auto voxel) {
		return voxel.x > cellsMX || voxel.x < 0 || voxel.y > cellsMY || voxel.y < 0;
	};
	glLineWidth(0.5f);

	auto nxt = voxel + glm::ivec2{ direction.x > 0 ? 1 : 0, direction.y > 0 ? 1 : 0 };
	auto nxtB = getVoxelPosition(nxt);
	auto tMax = glm::vec2{
		fabsf(direction.x) <= 1e-5f ? FLT_MAX : (nxtB.x - origin.x) / direction.x,
		fabsf(direction.y) <= 1e-5f ? FLT_MAX : (nxtB.y - origin.y) / direction.y
	};
	auto inside = gridPoint(vec(origin.x, origin.y), hScaled, gScaled) < isoLevel ? 0 : 1;
	//std::cout << "Ray started " << (inside ? "inside" : "outside") << std::endl;
	//std::cout << "Starting traversal at voxel [ " << voxel.x << " " << voxel.y << " ]" << std::endl;
	do {
		//std::cout << "Traversal at voxel [ " << voxel.x << " " << voxel.y << " ]" << std::endl;
		if (checkVoxel(voxel))
			break;
		auto vmin = getVoxelPosition(voxel);
		auto vmax = getVoxelPosition(voxel + glm::ivec2(1, 1));

		auto [hit, tmin, tmax] = rayIntersectAABB(origin, direction, vmin, vmax);
		tmin = std::max(0.f, tmin)+1e-6;
		if (tmax - 1e-6 > tmin) {
			std::vector<double> valArray;
			//std::cout << "Interval from " << tmin << " to " << tmax << std::endl;
			auto gridC = origin + direction * (tmax + tmin) * 0.5f;
			auto [xi, yi] = getCellIdx(gridC.x, gridC.y);
			auto range = (int32_t)1;

			auto gridPos = vec((scalar)origin.x + (tmax + tmin) * 0.5f * (scalar)direction.x, (scalar)origin.y + (tmax + tmin) * 0.5f * (scalar)direction.y);
			auto closest = hScaled * hScaled * 2;
			auto closestIdx = INT_MAX;
			for (int32_t ix = -range; ix <= range; ++ix) {
				for (int32_t iy = -range; iy <= range; ++iy) {
					if (xi + ix < 0 || xi + ix >= cellsX || yi + iy < 0 || yi + iy >= cellsY)
						continue;
					const auto& cellData = getCell(xi + ix, yi + iy);

					for (const auto& j : cellData) {
						if ((gridPos - particles[j].pos).norm() < closest) {
							closest = (gridPos - particles[j].pos).norm();
							closestIdx = j;
						}
					}
				}
			}
			auto sigma = 0.2;
			auto gauss = [=](scalar q) {
				return 1. / (sigma * std::sqrt(cheb::pi)) * std::exp(-.5 * square(q / sigma));
			};

			cheb::Function fn([=](cheb::scalar t) {
				auto gridPos = vec((scalar)origin.x + t * (scalar)direction.x, (scalar) origin.y + t * (scalar)direction.y);

				vec xAvg(0.0, 0.0);
				scalar wSum = 0.0;
				scalar value = 0.0;

				for (int32_t ix = -range; ix <= range; ++ix) {
					for (int32_t iy = -range; iy <= range; ++iy) {
						if (xi + ix < 0 || xi + ix >= cellsX || yi + iy < 0 || yi + iy >= cellsY)
							continue;
						const auto& cellData = getCell(xi + ix, yi + iy);
						
						for (const auto& j : cellData) {
							//if ((gridPos - particles[j].pos).norm() < closest)
							//	closest = (gridPos - particles[j].pos).norm();
							value += area * W(gridPos, particles[j].pos);
							float wi = gauss((gridPos - particles[j].pos).norm() / (hScaled));
							xAvg += wi * particles[j].pos;
							wSum += wi;
						}
					}
				}
				//return value + isoLevel;
				if (wSum > 0.0) {
					//return 0.0;
					xAvg /= wSum;
					return -(gridPos - xAvg).norm() - isoLevel;
				}
				else
					return closestIdx != INT_MAX ? -(gridPos - particles[closestIdx].pos).norm() - isoLevel : -1.0;


				}, { tmin, tmax - 1e-6 });
			if (fn.funs[0].size() > 1)
				printf("Chebyshev polynomial size: %d from %f to %f\n", fn.funs[0].size(), tmin, tmax);

			if (fn.funs[0].size() > 4000) {
				auto coeffs = fn.funs[0].coeffs();
				auto values = fn.funs[0].values();
				auto points = cheb::chebpts2(4097);
				for (int32_t i = 0; i < 4097; ++i)
					std::cout << std::setprecision(17) << points[i] << " " << values[i] << " " << coeffs[i] << std::endl;
			}
				/*std::cout << "Chebyshev size : " << fn.funs[0].size() << std::endl;*/
			auto roots = fn.roots();
			std::cout << "Number of roots on interval: " << roots.size() << std::endl;
			if (roots.size() > 0) {
				std::cout << "First root at: " << roots[0] << std::endl;
				auto pos = origin + glm::vec2(roots[0] * direction.x, roots[0] * direction.y);
				return std::make_pair(true, pos);
			}
		}

		if (tMax.x < tMax.y) {
			if (checkVoxel(voxel)) break;
			//std::cout << "Increment in x direction" << std::endl;
			voxel.x += sgn(direction.x);
			tMax.x += sgn(direction.x) / direction.x * gScale;
		}
		else {
			if (checkVoxel(voxel)) break;
			//std::cout << "Increment in x direction" << std::endl;
			voxel.y += sgn(direction.y);
			tMax.y += sgn(direction.y) / direction.y * gScale;
		}
		//std::cout << "Traversal step is done" << std::endl;
	} while (true);
	return std::make_pair(false, glm::vec2(0, 0));
}
#include <numbers>
void renderRays() {
	static auto& renderI = ParameterManager::instance().get<bool>("ray.renderImplicit");
	static auto& renderE = ParameterManager::instance().get<bool>("ray.renderExplicit");
	static auto& render = ParameterManager::instance().get<bool>("ray.render");
	if (!render) return;

	static auto& orig = ParameterManager::instance().get<vec>("ray.origin");
	static auto& target = ParameterManager::instance().get<vec>("ray.target");
	static auto& fov = ParameterManager::instance().get<scalar>("ray.fov");
	static auto& res = ParameterManager::instance().get<int32_t>("ray.resolution");
	static auto& isoLevel = ParameterManager::instance().get<scalar>("ray.iso");
	static auto& gScale = ParameterManager::instance().get<scalar>("ray.gScaleExplicit");
	static auto& hScale = ParameterManager::instance().get<scalar>("ray.hScale");
	static auto& subSteps = ParameterManager::instance().get<int32_t>("ray.subSteps");
	scalar gScaled = support * gScale;
	scalar hScaled = support * hScale;
	vec dir = target - orig;
	if (dir.norm() > 1e-5)
		dir.normalize();

	auto degToRad = [](auto angle) {
		return angle / 360.f * std::numbers::pi * 2.f;
	};
	glm::vec2 origin(orig.x(), orig.y());
	glm::vec2 direction(dir.x(), dir.y());

	auto fovRad = degToRad(fov);
	glm::vec2 rotated_p{ cos(fovRad / 2.0) * direction.x - sin(fovRad / 2.0) * direction.y,sin(fovRad / 2.0) * direction.x + cos(fovRad / 2.0) * direction.y };
	glm::vec2 rotated_n{ cos(-fovRad / 2.0) * direction.x - sin(-fovRad / 2.0) * direction.y,sin(-fovRad / 2.0) * direction.x + cos(-fovRad / 2.0) * direction.y };

	//std::cout << direction.x << " " << direction.y << " -> " << rotated_p.x << " " << rotated_p.y << " | " << rotated_n.x << " " << rotated_n.y << std::endl;

	glUseProgram(0);
	glLoadIdentity();
	glOrtho(0, domainWidth, 0, domainHeight, 0, 1);
	glBegin(GL_POINTS);
	glColor3f(1.f, 1.f, 1.f);
	glVertex2f(orig.x(), orig.y());
	glVertex2f(target.x(), target.y());
	glEnd();
	glLineWidth(1.5f);
	glBegin(GL_LINES);
	glColor3f(0.5f, 0.f, 0.f);
	glVertex2f(orig.x(), orig.y());
	glVertex2f(orig.x() + dir.x() * 200.0, orig.y() + dir.y() * 200.0);
	glEnd();
	glLineWidth(1.5f);

	glBegin(GL_LINES);
	glColor3f(0.f, 0.5f, 0.f);
	glVertex2f(orig.x(), orig.y());
	glVertex2f(orig.x() + rotated_p.x * 200.0, orig.y() + rotated_p.y * 200.0);
	glColor3f(0.f, 0.f, 0.5f);
	glVertex2f(orig.x(), orig.y());
	glVertex2f(orig.x() + rotated_n.x * 200.0, orig.y() + rotated_n.y * 200.0);
	glEnd();
	//if(renderI)
	//	intersectImplicitVis();
	if (renderE)
		intersectExplicitVis();

	static std::vector<glm::vec2> intersectionsExplicit;
	static std::vector<glm::vec2> intersectionsImplicit;
	static std::vector<glm::vec2> intersectionsCheb;
	static vec sOrig, sTarget;
	static 	scalar sFov, sIso, sG, sH;
	static int32_t sRes, sFrame, sSteps;
	static auto& frame = ParameterManager::instance().get<int32_t>("sim.frame");
	if (sSteps != subSteps || sOrig.x() != orig.x() || sTarget.x() != target.x() || sOrig.y() != orig.y() || sTarget.y() != target.y() || sFov != fov || sRes != res || sIso != isoLevel || sG != gScale || sH != hScale || sFrame != frame ) {
		sOrig = ParameterManager::instance().get<vec>("ray.origin");
		sTarget = ParameterManager::instance().get<vec>("ray.target");
		sFov = ParameterManager::instance().get<scalar>("ray.fov");
		sRes = ParameterManager::instance().get<int32_t>("ray.resolution");
		sIso = ParameterManager::instance().get<scalar>("ray.iso");
		sG = ParameterManager::instance().get<scalar>("ray.gScaleExplicit");
		sH = ParameterManager::instance().get<scalar>("ray.hScale");
		sFrame = frame;
		sSteps = subSteps;

		intersectionsCheb.clear();
		intersectionsImplicit.clear();
		intersectionsExplicit.clear();
		glm::vec2 direction = glm::vec2(dir.x(), dir.y());
		glm::vec2 origin = glm::vec2(orig.x(), orig.y());
		for (int32_t r = 0; r <= res; ++r) {
			float angle = -degToRad(fov / 2.0) + degToRad(fov) * (float)r / (float)res;
			glm::vec2 rDir{ cos(angle) * direction.x - sin(angle) * direction.y,sin(angle) * direction.x + cos(angle) * direction.y };
			auto [hit, pos ] = intersectExplicit(origin, rDir, hScaled, gScaled / 32);
			//std::cout << "Explicit: " << std::boolalpha << hit << " -> " << pos.x << " " << pos.y << std::endl;
			if (hit) intersectionsExplicit.push_back(pos);
		}
		for (int32_t r = 0; r <= res; ++r) {
			float angle = -degToRad(fov / 2.0) + degToRad(fov) * (float)r / (float)res;
			glm::vec2 rDir{ cos(angle) * direction.x - sin(angle) * direction.y,sin(angle) * direction.x + cos(angle) * direction.y };
			auto [hit, pos] = intersectImplicit(origin, rDir, hScaled, gScaled);
			//std::cout << "Implicit: " << std::boolalpha << hit << " -> " << pos.x << " " << pos.y << std::endl;
			if (hit) intersectionsImplicit.push_back(pos);
		}
////#pragma omp parallel for
//		for (int32_t r = 0; r <= res; ++r) {
//			float angle = -degToRad(fov / 2.0) + degToRad(fov) * (float)r / (float)res;
//			glm::vec2 rDir{ cos(angle) * direction.x - sin(angle) * direction.y,sin(angle) * direction.x + cos(angle) * direction.y };
//			auto [hit, pos] = intersectChebyshev(origin, rDir, hScaled, gScaled);
//			//std::cout << "Chebyshev: " << std::boolalpha << hit << " -> " << pos.x << " " << pos.y << std::endl;
////#pragma omp critical
//			if (hit) intersectionsCheb.push_back(pos);
//		}
	}
	auto& grid = ParameterManager::instance().get<bool>("render.showGrid");
	if (grid) {
		glLineWidth(0.5f);
		glColor4f(0.3f, 0.3f, 0.3f, 1);
		glBegin(GL_LINES);
		for (int32_t x = 1; x < (int32_t)cellsMX; x += 1) {
			for (int32_t y = 1; y < (int32_t)cellsMY; y += 1) {
				if ((y % 5 == 0)) {
					glLineWidth(1.f);
					glColor4f(0.7098f, 0.4941f, 0.8627f, 1);
				}
				glVertex2f(0.0, (float)y * gScaled);
				glVertex2f(domainWidth, (float)y * gScaled);
				if ((y % 5 == 0)) {
					glLineWidth(0.5f);
					glColor4f(0.5882f, 0.2824f, 0.8039f, 1);
				}
			}
			if ((x % 5 == 0)) {
				glLineWidth(1.f);
				glColor4f(0.7098f, 0.4941f, 0.8627f, 1);
			}
			glVertex2f((float)x * gScaled, 0.0);
			glVertex2f((float)x * gScaled, domainHeight);

			if ((x % 5 == 0)) {
				glLineWidth(0.5f);
				glColor4f(0.5882f, 0.2824f, 0.8039f, 1);
			}
		}
		glEnd();
	}
	glPointSize(8.f);
	glBegin(GL_POINTS);
	glColor3f(0.f, 1.f, 0.498f); // guppie green
	if(renderI)
	for (const auto& p : intersectionsImplicit) {
		//std::cout << p.x << " " << p.y << " \n";
		glVertex2f(p.x, p.y);
	}
	glColor3f(1.f, 0.f, 0.4235f); // vivid raspberry
	if(renderE)
	for (const auto& p : intersectionsExplicit) {
		//std::cout << p.x << " " << p.y << " \n";
		glVertex2f(p.x, p.y);
	}
	glColor3f(0.4235f, 0.f, 1.f); // vivid raspberry
	if (true)
		for (const auto& p : intersectionsCheb) {
			//std::cout << p.x << " " << p.y << " \n";
			glVertex2f(p.x, p.y);
		}
	glEnd();
	glPointSize(16.f * scale * pointScale/** 5.f*/);

}