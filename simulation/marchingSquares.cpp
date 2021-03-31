#include <glad/glad.h>
#include <GLFW/glfw3.h>
// include order for gl has to be left like this
#include "SPH.h"
#include "Math.h"
#include "tinycolormap.h"
#include "SPHRender.h"
#include <chrono>
#include <iostream>
#include <sstream>
#include <vector>



#define PAR_MSQUARES_IMPLEMENTATION
#include <tools/par_msquares.h>
#include <glm/glm.hpp>
struct Mesh {
	std::vector<glm::vec2> positions;
	std::vector<glm::ivec3> triangles;
};

auto updateMarchingGrid() {
	static auto& isoLevel = ParameterManager::instance().get<scalar>("marching.iso");
	static auto& gScale = ParameterManager::instance().get<scalar>("marching.gScale");
	static auto& hScale = ParameterManager::instance().get<scalar>("marching.hScale");
	static auto& method = ParameterManager::instance().get<int32_t>("marching.method");

	scalar h = support;
	scalar hScaled = h * hScale;
	scalar gScaled = h * gScale;

	cellsMX = (std::size_t) ::ceil((domainWidth) / gScaled);
	cellsMY = (std::size_t) ::ceil((domainHeight) / gScaled);
	float* data = (float*)malloc(sizeof(float) * cellsMX * cellsMY);

	std::cout << "Starting grid evaluation of dimension " << cellsMX << " x " << cellsMY << std::endl;
	std::cout << "For scale parameters: " << hScaled << " x " << gScaled << std::endl;
#pragma omp parallel for
	for (int32_t x = 0; x < cellsMX; ++x) {
		for (int32_t y = 0; y < cellsMY; ++y) {
			scalar xPos = 0.0 + (scalar)x / (scalar)cellsMX * (scalar)domainWidth;
			scalar yPos = 0.0 + (scalar)y / (scalar)cellsMY * (scalar)domainHeight;
			data[x + y * cellsMX] = (float)gridPoint(vec(xPos, yPos), hScaled, gScaled);
		}
	}
	std::cout << "Running marching squares" << std::endl;

	// https://en.wikipedia.org/wiki/Marching_squares
	auto getGridPos = [&](auto x, auto y) {
		scalar xPos = 0.0 + (scalar)x / (scalar)cellsMX * (scalar)domainWidth;
		scalar yPos = 0.0 + (scalar)y / (scalar)cellsMY * (scalar)domainHeight;
		return std::make_tuple(glm::vec2(xPos, yPos), x + y * cellsMX, data[x + y * cellsMX]);
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

	std::vector<glm::vec2> positions;
	std::vector<glm::ivec3> triangles;
#define pushTrianglePositions(p1, p2 ,p3)\
	positions.push_back(p1);\
	positions.push_back(p2);\
	positions.push_back(p3);\
	//auto pushTrianglePositions = [&](const auto& p1, const auto& p2, const auto& p3) {
	//	positions.push_back(p1);
	//	positions.push_back(p2);
	//	positions.push_back(p3);
	//};
	for (int32_t x = 0; x < cellsMX - 1; ++x) {
		for (int32_t y = 0; y < cellsMY - 1; ++y) {
			
			auto [posDL, idxDL, sDL] = getGridPos(x, y);
			auto [posDR, idxDR, sDR] = getGridPos(x+1, y);
			auto [posTR, idxTR, sTR] = getGridPos(x+1, y+1);
			auto [posTL, idxTL, sTL] = getGridPos(x, y+1);
			auto lerpL = vertexInterp(posDL, posTL, sDL, sTL);
			auto lerpR = vertexInterp(posDR, posTR, sDR, sTR);
			auto lerpD = vertexInterp(posDL, posDR, sDL, sDR);
			auto lerpT = vertexInterp(posTL, posTR, sTL, sTR);

			auto march = (sDL > isoLevel ? 0x1 : 0x0) | (sDR > isoLevel ? 0x2 : 0x0) | (sTR > isoLevel ? 0x4 : 0x0) | (sTL > isoLevel ? 0x8 : 0x0);

			//std::cout << "Grid cell " << x << " x " << y << "\n" <<
			//	sTL << " --- " << sTR << "\n|     |\n|     |\n" <<
			//	sDL << " --- " << sDR << "\n" <<
			//	march << " [ " << (sDL > isoLevel ? 0x1 : 0x0) << (sDR > isoLevel ? 0x2 : 0x0) << (sTR > isoLevel ? 0x4 : 0x0) << (sTL > isoLevel ? 0x8 : 0x0) << "]" << std::endl;
			const auto& DR = posDR;
			const auto& DL = posDL;
			const auto& TR = posTR;
			const auto& TL = posTL;

			//std::cout << x << " [" << cellsMX << "] x " << y << " [" << cellsMY << "] -> " << positions.size() << " -> " << march << std::endl;

			switch (march) {
			case 0: break;
			case 1: {
				// Marching Squares case:
				// 0 --- 0
				// |     |
				// |     |
				// 1 --- 0
				// Output: 1 Triangle (LR + lerpD + lerpL)
				// 1st triangle
				//pushTrianglePositions(posDL, lerpD, lerpL);
				pushTrianglePositions(posDL, lerpD, lerpL);
			}break;
			case 2: {
				// Marching Squares case:
				// 0 --- 0
				// |     |
				// |     |
				// 0 --- 1
				// Output: 1 Triangle (DR + lerpR + lerpD)
				pushTrianglePositions(DR, lerpR, lerpD);
			} break;
			case 3: { 
				// Marching Squares case:
				// 0 --- 0
				// |     |
				// |     |
				// 1 --- 1
				// Output: 2 Triangles
				// 1:		(DL + DR    + lerpL)
				// 2:		(DR + lerpR + lerpL)
				pushTrianglePositions(DL, DR, lerpL);
				pushTrianglePositions(DR, lerpR, lerpL);
			} break;
			case 4: {
				// Marching Squares case:
				// 0 --- 1
				// |     |
				// |     |
				// 0 --- 0
				// Output: 1 Triangle (TR + lerpT + lerpR)
				pushTrianglePositions(TR, lerpT, lerpR);
			} break;
			case 5: { 
				// Marching Squares case:
				// 0 --- 1
				// |     |
				// |     |
				// 1 --- 0
				// Output: 4 Triangles
				// 1:		(DL + lerpD + lerpL)
				// 2:		(lerpD + lerpR + lerpL)
				// 3:		(lerpL + lerpR + lerpT)
				// 4:		(TR + lerpT + lerpR)
				pushTrianglePositions(DL, lerpD, lerpL);
				if (sDL + sDR + sTL + sTR > 4.0 * isoLevel) {
					pushTrianglePositions(lerpD, lerpR, lerpL);
					pushTrianglePositions(lerpL, lerpR, lerpT);
				}
				pushTrianglePositions(TR, lerpT, lerpR);
			} break;
			case 6: {
				// Marching Squares case:
				// 0 --- 1
				// |     |
				// |     |
				// 0 --- 1
				// Output: 2 Triangles
				// 1:		(DR + lerpT + lerpD)
				// 2:		(TR + lerpT + DR   )
				pushTrianglePositions(DR, lerpT, lerpD);
				pushTrianglePositions(TR, lerpT, DR);
			} break;
			case 7: { 
				// Marching Squares case:
				// 0 --- 1
				// |     |
				// |     |
				// 1 --- 1
				// Output: 3 Triangles
				// 1:		(DL + DR    + TR   )
				// 2:		(DL + lerpT + lerpL)
				// 3:		(TR + lerpT + DL   )
				pushTrianglePositions(DL, DR, TR);
				pushTrianglePositions(DL, lerpT, lerpL);
				pushTrianglePositions(TR, lerpT, DL);
			} break;
			case 8: {
				// Marching Squares case:
				// 1 --- 0
				// |     |
				// |     |
				// 0 --- 0
				// Output: 1 Triangle (TR + lerpL + lerpT)
				pushTrianglePositions(TL, lerpL, lerpT);
			} break;
			case 9: { 
				// Marching Squares case:
				// 1 --- 0
				// |     |
				// |     |
				// 1 --- 0
				// Output: 2 Triangles
				// 1:		(DL + lerpD + lerpT)
				// 2:		(TL + DL    + lerpT)
				pushTrianglePositions(DL, lerpD, lerpT);
				pushTrianglePositions(TL, DL, lerpT);
			} break;
			case 10: {
				// Marching Squares case:
				// 1 --- 0
				// |     |
				// |     |
				// 0 --- 1
				// Output: 4 Triangles
				// 1:		(DR + lerpR + lerpD)
				// 2:		(lerpR + lerpT + lerpD)
				// 3:		(lerpT + lerpL + lerpD)
				// 4:		(TL + lerpL + lerpT)
				pushTrianglePositions(DR, lerpR, lerpD);
				if (sDL + sDR + sTL + sTR > 4.0 * isoLevel) {
					pushTrianglePositions(lerpR, lerpT, lerpD);
					pushTrianglePositions(lerpT, lerpL, lerpD);
				}
				pushTrianglePositions(TL, lerpL, lerpT);
			} break;
			case 11: { 
				// Marching Squares case:
				// 1 --- 0
				// |     |
				// |     |
				// 1 --- 1
				// Output: 3 Triangles
				// 1:		(DL + DR    + TL   )
				// 2:		(TL + lerpR + lerpT)
				// 3:		(DR + lerpR + TL   )
				pushTrianglePositions(DL, DR, TL);
				pushTrianglePositions(TL, lerpR, lerpT);
				pushTrianglePositions(DR, lerpR, TL);
			} break;
			case 12: { 
				// Marching Squares case:
				// 1 --- 1
				// |     |
				// |     |
				// 0 --- 0
				// Output: 2 Triangles
				// 1:		(TR + TL    + lerpL)
				// 2:		(TR + lerpL + lerpR)
				pushTrianglePositions(TR, TL, lerpL);
				pushTrianglePositions(TR, lerpL, lerpR);
			} break;
			case 13: { 
				// Marching Squares case:
				// 1 --- 1
				// |     |
				// |     |
				// 1 --- 0
				// Output: 3 Triangles
				// 1:		(DL + TR    + TL   )
				// 2:		(DL + lerpD + lerpR)
				// 3:		(DL + lerpR + TL   )
				pushTrianglePositions(DL, TR, TL);
				pushTrianglePositions(DL, lerpD, lerpR);
				pushTrianglePositions(DL, lerpR, TL);
			} break;
			case 14: {
				// Marching Squares case:
				// 1 --- 1
				// |     |
				// |     |
				// 0 --- 1
				// Output: 3 Triangles
				// 1:		(TR + TL    + DR   )
				// 2:		(TL + lerpL + lerpD)
				// 3:		(DR + TL    + lerpD)
				pushTrianglePositions(TR, TL, DR);
				pushTrianglePositions(TL, lerpL, lerpD);
				pushTrianglePositions(DR, TL, lerpD);
			} break;
			case 15: {
				// Marching Squares case:
				// 1 --- 1
				// |     |
				// |     |
				// 1 --- 1
				// Output: 2 Triangles
				// 1:		(DL + DR    + TL   )
				// 2:		(DR + TR    + TL   )
				pushTrianglePositions(posDL, posDR, posTL);
				pushTrianglePositions(posDR, posTR, posTL);
			} break;
			}
		}
	}
	for (int32_t i = 0; i < positions.size() / 3; ++i) {
		triangles.push_back(glm::ivec3{ i * 3, i * 3 + 1, i * 3 + 2 });
	}




	std::cout << "Finished marching squares" << std::endl;
	Mesh m{positions,triangles};
	//std::vector<const par_msquares_mesh*> meshes;
	//for (auto i = 0; i < mesh_ptr->nmeshes; ++i) {
	//    std::cout << "Extracting mesh " << i << std::endl;
	//    meshes.push_back(par_msquares_get_mesh(mesh_ptr, i));
	//    for (int32_t j = 0; j < meshes[i]->npoints; ++j) {
	//       // meshes[i]->points[j] *=  ((scalar)std::max(cellsMX, cellsMY));
	//    }
	//}

	return std::make_pair(data, m);
}
auto updateMarchingOutline() {
	static auto& isoLevel = ParameterManager::instance().get<scalar>("marching.iso");
	static auto& gScale = ParameterManager::instance().get<scalar>("marching.gScale");
	static auto& hScale = ParameterManager::instance().get<scalar>("marching.hScale");
	static auto& method = ParameterManager::instance().get<int32_t>("marching.method");

	scalar h = support;
	scalar hScaled = h * hScale;
	scalar gScaled = h * gScale;

	cellsMX = (std::size_t) ::ceil((domainWidth) / gScaled);
	cellsMY = (std::size_t) ::ceil((domainHeight) / gScaled);
	float* data = (float*)malloc(sizeof(float) * cellsMX * cellsMY);

	std::cout << "Starting grid evaluation of dimension " << cellsMX << " x " << cellsMY << std::endl;
	std::cout << "For scale parameters: " << hScaled << " x " << gScaled << std::endl;
#pragma omp parallel for
	for (int32_t x = 0; x < cellsMX; ++x) {
		for (int32_t y = 0; y < cellsMY; ++y) {
			scalar xPos = 0.0 + (scalar)x / (scalar)cellsMX * (scalar)domainWidth;
			scalar yPos = 0.0 + (scalar)y / (scalar)cellsMY * (scalar)domainHeight;
			data[x + y * cellsMX] = (float)gridPoint(vec(xPos, yPos), hScaled, gScaled);
		}
	}
	std::cout << "Running marching squares" << std::endl;

	// https://en.wikipedia.org/wiki/Marching_squares
	auto getGridPos = [&](auto x, auto y) {
		scalar xPos = 0.0 + (scalar)x / (scalar)cellsMX * (scalar)domainWidth;
		scalar yPos = 0.0 + (scalar)y / (scalar)cellsMY * (scalar)domainHeight;
		return std::make_tuple(glm::vec2(xPos, yPos), x + y * cellsMX, data[x + y * cellsMX]);
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

	std::vector<glm::vec2> positions;
	std::vector<glm::ivec3> triangles;
#define pushLinePositions(p1, p2)\
	positions.push_back(p1);\
	positions.push_back(p2);
	//auto pushTrianglePositions = [&](const auto& p1, const auto& p2, const auto& p3) {
	//	positions.push_back(p1);
	//	positions.push_back(p2);
	//	positions.push_back(p3);
	//};
	for (int32_t x = 0; x < cellsMX - 1; ++x) {
		for (int32_t y = 0; y < cellsMY - 1; ++y) {

			auto [posDL, idxDL, sDL] = getGridPos(x, y);
			auto [posDR, idxDR, sDR] = getGridPos(x + 1, y);
			auto [posTR, idxTR, sTR] = getGridPos(x + 1, y + 1);
			auto [posTL, idxTL, sTL] = getGridPos(x, y + 1);
			auto lerpL = vertexInterp(posDL, posTL, sDL, sTL);
			auto lerpR = vertexInterp(posDR, posTR, sDR, sTR);
			auto lerpD = vertexInterp(posDL, posDR, sDL, sDR);
			auto lerpT = vertexInterp(posTL, posTR, sTL, sTR);

			auto march = (sDL > isoLevel ? 0x1 : 0x0) | (sDR > isoLevel ? 0x2 : 0x0) | (sTR > isoLevel ? 0x4 : 0x0) | (sTL > isoLevel ? 0x8 : 0x0);

			//std::cout << "Grid cell " << x << " x " << y << "\n" <<
			//	sTL << " --- " << sTR << "\n|     |\n|     |\n" <<
			//	sDL << " --- " << sDR << "\n" <<
			//	march << " [ " << (sDL > isoLevel ? 0x1 : 0x0) << (sDR > isoLevel ? 0x2 : 0x0) << (sTR > isoLevel ? 0x4 : 0x0) << (sTL > isoLevel ? 0x8 : 0x0) << "]" << std::endl;
			const auto& DR = posDR;
			const auto& DL = posDL;
			const auto& TR = posTR;
			const auto& TL = posTL;

			//std::cout << x << " [" << cellsMX << "] x " << y << " [" << cellsMY << "] -> " << positions.size() << " -> " << march << std::endl;

			switch (march) {
			case 0: break;
			case 1: {
				// Marching Squares case:
				// 0 --- 0
				// |     |
				// |     |
				// 1 --- 0
				// Output: 1 Triangle (LR + lerpD + lerpL)
				// 1st triangle
				//pushTrianglePositions(posDL, lerpD, lerpL);
				pushLinePositions(lerpD, lerpL);
			}break;
			case 2: {
				// Marching Squares case:
				// 0 --- 0
				// |     |
				// |     |
				// 0 --- 1
				// Output: 1 Triangle (DR + lerpR + lerpD)
				pushLinePositions(lerpR, lerpD);
			} break;
			case 3: {
				// Marching Squares case:
				// 0 --- 0
				// |     |
				// |     |
				// 1 --- 1
				// Output: 2 Triangles
				// 1:		(DL + DR    + lerpL)
				// 2:		(DR + lerpR + lerpL)
				pushLinePositions(lerpR, lerpL);
			} break;
			case 4: {
				// Marching Squares case:
				// 0 --- 1
				// |     |
				// |     |
				// 0 --- 0
				// Output: 1 Triangle (TR + lerpT + lerpR)
				pushLinePositions(lerpT, lerpR);
			} break;
			case 5: {
				// Marching Squares case:
				// 0 --- 1
				// |     |
				// |     |
				// 1 --- 0
				// Output: 4 Triangles
				// 1:		(DL + lerpD + lerpL)
				// 2:		(lerpD + lerpR + lerpL)
				// 3:		(lerpL + lerpR + lerpT)
				// 4:		(TR + lerpT + lerpR)
				if (sDL + sDR + sTL + sTR > 4.0 * isoLevel) {
					pushLinePositions(lerpD, lerpR);
					pushLinePositions(lerpL, lerpT);
				}
				else {
					pushLinePositions(lerpL, lerpD);
					pushLinePositions(lerpT, lerpR);
				}
			} break;
			case 6: {
				// Marching Squares case:
				// 0 --- 1
				// |     |
				// |     |
				// 0 --- 1
				// Output: 2 Triangles
				// 1:		(DR + lerpT + lerpD)
				// 2:		(TR + lerpT + DR   )
				pushLinePositions(lerpT, lerpD);
			} break;
			case 7: {
				// Marching Squares case:
				// 0 --- 1
				// |     |
				// |     |
				// 1 --- 1
				// Output: 3 Triangles
				// 1:		(DL + DR    + TR   )
				// 2:		(DL + lerpT + lerpL)
				// 3:		(TR + lerpT + DL   )
				pushLinePositions(lerpT, lerpL);
			} break;
			case 8: {
				// Marching Squares case:
				// 1 --- 0
				// |     |
				// |     |
				// 0 --- 0
				// Output: 1 Triangle (TR + lerpL + lerpT)
				pushLinePositions(lerpL, lerpT);
			} break;
			case 9: {
				// Marching Squares case:
				// 1 --- 0
				// |     |
				// |     |
				// 1 --- 0
				// Output: 2 Triangles
				// 1:		(DL + lerpD + lerpT)
				// 2:		(TL + DL    + lerpT)
				pushLinePositions(lerpD, lerpT);
			} break;
			case 10: {
				// Marching Squares case:
				// 1 --- 0
				// |     |
				// |     |
				// 0 --- 1
				// Output: 4 Triangles
				// 1:		(DR + lerpR + lerpD)
				// 2:		(lerpR + lerpT + lerpD)
				// 3:		(lerpT + lerpL + lerpD)
				// 4:		(TL + lerpL + lerpT)
				if (sDL + sDR + sTL + sTR > 4.0 * isoLevel) {
					pushLinePositions(lerpD, lerpL);
					pushLinePositions(lerpR, lerpT);
				}
				else {
					pushLinePositions(lerpL, lerpT);
					pushLinePositions(lerpD, lerpR);
				}
			} break;
			case 11: {
				// Marching Squares case:
				// 1 --- 0
				// |     |
				// |     |
				// 1 --- 1
				// Output: 3 Triangles
				// 1:		(DL + DR    + TL   )
				// 2:		(TL + lerpR + lerpT)
				// 3:		(DR + lerpR + TL   )
				pushLinePositions(lerpR, lerpT);
			} break;
			case 12: {
				// Marching Squares case:
				// 1 --- 1
				// |     |
				// |     |
				// 0 --- 0
				// Output: 2 Triangles
				// 1:		(TR + TL    + lerpL)
				// 2:		(TR + lerpL + lerpR)
				pushLinePositions(lerpL, lerpR);
			} break;
			case 13: {
				// Marching Squares case:
				// 1 --- 1
				// |     |
				// |     |
				// 1 --- 0
				// Output: 3 Triangles
				// 1:		(DL + TR    + TL   )
				// 2:		(DL + lerpD + lerpR)
				// 3:		(DL + lerpR + TL   )
				pushLinePositions(lerpD, lerpR);
			} break;
			case 14: {
				// Marching Squares case:
				// 1 --- 1
				// |     |
				// |     |
				// 0 --- 1
				// Output: 3 Triangles
				// 1:		(TR + TL    + DR   )
				// 2:		(TL + lerpL + lerpD)
				// 3:		(DR + TL    + lerpD)
				pushLinePositions(lerpL, lerpD);
			} break;
			case 15: {
				// Marching Squares case:
				// 1 --- 1
				// |     |
				// |     |
				// 1 --- 1
				// Output: 2 Triangles
				// 1:		(DL + DR    + TL   )
				// 2:		(DR + TR    + TL   )
			} break;
			}
		}
	}
	free(data);
	std::cout << "Finished marching squares" << std::endl;
	return positions;




}
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

void renderMarchingOutline() {
	static scalar g = 0.0, h = 0.0, i = 0.0;
	static int32_t m = -1, f = -1;
	static auto& iso = ParameterManager::instance().get<scalar>("marching.iso");
	static auto& gScale = ParameterManager::instance().get<scalar>("marching.gScale");
	static auto& hScale = ParameterManager::instance().get<scalar>("marching.hScale");
	static auto& method = ParameterManager::instance().get<int32_t>("marching.method");
	static auto& frame = ParameterManager::instance().get<int32_t>("sim.frame");
	static std::vector<glm::vec2> sMesh;
	static GLuint VAO = 0xDEADBEEF, VBO;
	static GLuint program = 0xDEADBEEF;
	if (frame != f || method != m || hScale != h || gScale != g || iso != i) {
		i = iso;
		f = frame;
		g = gScale;
		h = hScale;
		m = method;
		std::cout << "Updating marching squares" << std::endl;
		auto mesh = updateMarchingOutline();
		if (VAO != 0xDEADBEEF) {
			glDeleteVertexArrays(1, &VAO);
			glDeleteBuffers(1, &VBO);
		}
		sMesh = mesh;

		glGenVertexArrays(1, &VAO);
		glGenBuffers(1, &VBO);

		std::cout << "Mesh " << 0 << " has " << sMesh.size() /2 << " lines and " << sMesh.size() << " vertices" << std::endl;

		//for (int32_t t = 0; t < sMesh.triangles.size(); ++t) {
		//	std::cout << "Triangle " << t << " -> [ " << sMesh.triangles[t].y << " " << sMesh.triangles[t].x << " " << sMesh.triangles[t].z << " ] \n"
		//		<< "vtx0: " << sMesh.positions[sMesh.triangles[t].x].x << " x " << sMesh.positions[sMesh.triangles[t].x].y << "\n"
		//		<< "vtx0: " << sMesh.positions[sMesh.triangles[t].y].x << " x " << sMesh.positions[sMesh.triangles[t].y].y << "\n"
		//		<< "vtx0: " << sMesh.positions[sMesh.triangles[t].z].x << " x " << sMesh.positions[sMesh.triangles[t].z].y << "\n"
		//		;
		//}
		glBindVertexArray(VAO);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, sMesh.size() * sizeof(float) * 2, sMesh.data(), GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glVertexAttribPointer(
			0,                  // attribute
			2,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);
		// Index buffer

	}
	if (program == 0xDEADBEEF) {

		static const GLchar* vertex_shader_source =
			R"(#version 150 
in vec2 position;
uniform mat4 transform;
void main() {
    gl_Position = transform * vec4(position, 0.f, 1.0f);
}
)";
		static const GLchar* fragment_shader_source =
			R"(#version 150
out vec4 color;
void main() {
    color = vec4(1.f,0.5f,0.2f, 1.0f);
};
)";
		program = createProgram(vertex_shader_source, fragment_shader_source);
		glUseProgram(program);
		auto loc = glGetUniformLocation(program, "transform");
		glUniformMatrix4fv(loc, 1, GL_FALSE, &glm::ortho(0.f, (float)domainWidth, 0.f, (float)domainHeight, 0.f,1.f)[0][0]);
		glUseProgram(0);
	}
	glUseProgram(program);


	glBindVertexArray(VAO);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//std::cout << "Rendering VAO " << m << " with " << sMeshes[m]->ntriangles * 3 << " indices " << std::endl;
	glDrawArrays(
		GL_LINES,      // mode
		0,
		sMesh.size() * 2
	);	
	//glDrawElements(
	//	GL_TRIANGLES,      // mode
	//	sMesh.triangles.size() * 3,    // count
	//	GL_UNSIGNED_SHORT,   // type
	//	(void*)0           // element array buffer offset
	//);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBindVertexArray(0);

	glUseProgram(0);
}


void renderMarchingSolid() {
	static scalar g = 0.0, h = 0.0, i = 0.0;
	static int32_t m = -1, f = -1;
	static auto& iso = ParameterManager::instance().get<scalar>("marching.iso");
	static auto& gScale = ParameterManager::instance().get<scalar>("marching.gScale");
	static auto& hScale = ParameterManager::instance().get<scalar>("marching.hScale");
	static auto& method = ParameterManager::instance().get<int32_t>("marching.method");
	static auto& frame = ParameterManager::instance().get<int32_t>("sim.frame");
	static Mesh sMesh;
	static GLuint VAO = 0xDEADBEEF, VBO, IBO;
	static GLuint program = 0xDEADBEEF;
	if (frame != f || method != m || hScale != h || gScale != g || iso != i) {
		i = iso;
		f = frame;
		g = gScale;
		h = hScale;
		m = method;
		std::cout << "Updating marching squares" << std::endl;
		auto [data, mesh] = updateMarchingGrid();
		if (VAO != 0xDEADBEEF) {
			glDeleteVertexArrays(1, &VAO);
			glDeleteBuffers(1, &VBO);
			glDeleteBuffers(1, &IBO);
		}
		sMesh = mesh;

		glGenVertexArrays(1, &VAO);
		glGenBuffers(1, &VBO);
		glGenBuffers(1, &IBO);

		std::cout << "Mesh " << 0 << " has " << sMesh.triangles.size() << " triangles and " << sMesh.positions.size() << " vertices" << std::endl;

		//for (int32_t t = 0; t < sMesh.triangles.size(); ++t) {
		//	std::cout << "Triangle " << t << " -> [ " << sMesh.triangles[t].y << " " << sMesh.triangles[t].x << " " << sMesh.triangles[t].z << " ] \n"
		//		<< "vtx0: " << sMesh.positions[sMesh.triangles[t].x].x << " x " << sMesh.positions[sMesh.triangles[t].x].y << "\n"
		//		<< "vtx0: " << sMesh.positions[sMesh.triangles[t].y].x << " x " << sMesh.positions[sMesh.triangles[t].y].y << "\n"
		//		<< "vtx0: " << sMesh.positions[sMesh.triangles[t].z].x << " x " << sMesh.positions[sMesh.triangles[t].z].y << "\n"
		//		;
		//}
		glBindVertexArray(VAO);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, sMesh.triangles.size() * 3 * sizeof(int32_t), sMesh.triangles.data(), GL_STATIC_DRAW);

		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glBufferData(GL_ARRAY_BUFFER, sMesh.positions.size() * sizeof(float) * 2, sMesh.positions.data(), GL_STATIC_DRAW);

		glEnableVertexAttribArray(0);
		glBindBuffer(GL_ARRAY_BUFFER, VBO);
		glVertexAttribPointer(
			0,                  // attribute
			2,                  // size
			GL_FLOAT,           // type
			GL_FALSE,           // normalized?
			0,                  // stride
			(void*)0            // array buffer offset
		);
		// Index buffer
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IBO);
		glBindVertexArray(0);

		free(data);
	}
	if (program == 0xDEADBEEF) {

		static const GLchar* vertex_shader_source =
			R"(#version 150 
in vec2 position;
uniform mat4 transform;
void main() {
    gl_Position = transform * vec4(position, 0.f, 1.0f);
}
)";
		static const GLchar* fragment_shader_source =
			R"(#version 150
out vec4 color;
void main() {
    color = vec4(1.f,0.5f,0.2f, 1.0f);
};
)";
		program = createProgram(vertex_shader_source, fragment_shader_source);
		glUseProgram(program);
		auto loc = glGetUniformLocation(program, "transform");
		glUniformMatrix4fv(loc, 1, GL_FALSE, &glm::ortho(0.f, (float)domainWidth, 0.f, (float)domainHeight, 0.f, 1.f)[0][0]);
		glUseProgram(0);
	}
	glUseProgram(program);


	glBindVertexArray(VAO);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//std::cout << "Rendering VAO " << m << " with " << sMeshes[m]->ntriangles * 3 << " indices " << std::endl;
	glDrawElements(
		GL_TRIANGLES,      // mode
		sMesh.triangles.size() * 3,    // count
		GL_UNSIGNED_INT,   // type
		(void*)0           // element array buffer offset
	);
	//glDrawElements(
	//	GL_TRIANGLES,      // mode
	//	sMesh.triangles.size() * 3,    // count
	//	GL_UNSIGNED_SHORT,   // type
	//	(void*)0           // element array buffer offset
	//);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBindVertexArray(0);

	glUseProgram(0);
}


void renderMarching() {
	static auto& solid = ParameterManager::instance().get<bool>("marching.solid");
	if (solid)
		renderMarchingSolid();
	else 
		renderMarchingOutline();
}