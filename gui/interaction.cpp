#include "glui.h"
#include <imgui/imgui.h>
#include <iostream>

void GUI::keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
	static bool emit = false;
	//if (emit) emitParticles();


	if (action != GLFW_PRESS) { 
		if (key == GLFW_KEY_E && action == GLFW_RELEASE)
			emit = false;
		
		return; }
	switch (key) {
	case GLFW_KEY_P:
		if (mods & GLFW_MOD_SHIFT)
			simulationState.timestep();
		else
			simulationRunning = !simulationRunning;
		break;
	case GLFW_KEY_H:
		m_showText = !m_showText;
		break;
	case GLFW_KEY_G:
	{
		auto& grid = simulationState.pm.get<bool>("render.showGrid");
		grid = !grid;
		break;
	}
	case GLFW_KEY_O:
	{
		auto& grid = simulationState.pm.get<bool>("render.showPtcls");
		grid = !grid;
		break;
	}
	case GLFW_KEY_1: simulationState.pm.get<int32_t>("colorMap.map") = 0; break;
	case GLFW_KEY_2: simulationState.pm.get<int32_t>("colorMap.map") = 1; break;
	case GLFW_KEY_3: simulationState.pm.get<int32_t>("colorMap.map") = 2; break;
	case GLFW_KEY_T: simulationState.pm.get<bool>("colorMap.auto") = !simulationState.pm.get<bool>("colorMap.auto"); break;
	case GLFW_KEY_F: simulationState.pm.get<bool>("field.render") = !simulationState.pm.get<bool>("field.render"); break;
	case GLFW_KEY_M: simulationState.pm.get<bool>("marching.render") = !simulationState.pm.get<bool>("marching.render"); break;
	case GLFW_KEY_N: simulationState.pm.get<bool>("marching.solid") = !simulationState.pm.get<bool>("marching.solid"); break;
	
	case GLFW_KEY_E: emit=true; break;
	}
}
bool trackingLeft = false;
bool trackingRight = false;

void GUI::mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
	if(action == GLFW_RELEASE) {
		if (button == GLFW_MOUSE_BUTTON_LEFT)
			trackingLeft = false;
		if (button == GLFW_MOUSE_BUTTON_RIGHT)
			trackingRight = false;
	}
}

//#include <glrender/glparticleIndexRender/particleIndexRender.h>
#include <imgui/imgui_internal.h>
void GUI::cursorPositionCallback(GLFWwindow* window, double xpos, double ypos) {
  m_cursorPosition = vec(xpos, ypos);



}
void GUI::scrollCallback(GLFWwindow* window, double xoffset, double yoffset) { }
