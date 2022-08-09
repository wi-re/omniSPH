#include "glui.h"
void GUI::OSD(){
    if (m_showText) {
        const float DISTANCE = 10.0f;
        static int corner = 0;
        ImGuiIO& io = ImGui::GetIO();
        if (corner != -1)
        {
            ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImVec2 work_area_pos = viewport->GetWorkPos();   // Instead of using viewport->Pos we use GetWorkPos() to avoid menu bars, if any!
            ImVec2 work_area_size = viewport->GetWorkSize();
            ImVec2 window_pos = ImVec2((corner & 1) ? (work_area_pos.x + work_area_size.x - DISTANCE) : (work_area_pos.x + DISTANCE), (corner & 2) ? (work_area_pos.y + work_area_size.y - DISTANCE) : (work_area_pos.y + DISTANCE));
            ImVec2 window_pos_pivot = ImVec2((corner & 1) ? 1.0f : 0.0f, (corner & 2) ? 1.0f : 0.0f);
            ImGui::SetNextWindowPos(window_pos, ImGuiCond_Always, window_pos_pivot);
            ImGui::SetNextWindowViewport(viewport->ID);
        }
        ImGui::SetNextWindowBgAlpha(0.35f); // Transparent background
        if (ImGui::Begin("Stats overlay", &m_showText, (corner != -1 ? ImGuiWindowFlags_NoMove : 0) | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav))
        {
            auto addParameter = [&](auto paramName) {
                auto& param = simulationState.pm.getParameter(paramName);
                ImGui::PushID(param.identifier.c_str());
                simulationState.pm.uiFunctions[param.type](param);
                ImGui::PopID();
            };

            addParameter("time");
            addParameter("dt");
            //addParameter("simTime");
            //addParameter("renderTime");
            addParameter("numPtcls");
            addParameter("densityIterations");
            addParameter("densityError");
            addParameter("divergenceIterations");
            addParameter("divergenceError");
            //addParameter("color_map.buffer");
            addParameter("colorMap.min");
            addParameter("colorMap.max");
            //if (dualView) {
            //    addParameter("colorMap.minDual");
            //    addParameter("colorMap.maxDual");

            //}

        }
        ImGui::End();
    }
    if(m_showText){
        static auto& numPtcls = simulationState.pm.get<int32_t>("props.numPtcls");
        if (numPtcls == 0 ||
            simulationState.pm.get<int32_t>("sim.frame") == 0)return;


        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);

        auto x = std::clamp(m_cursorPosition.x(), 0.0, (scalar)display_w) / (scalar)display_w;
        auto y =1. * std::clamp(m_cursorPosition.y(), 0.0, (scalar)display_h) / (scalar)display_h;
        if (y > 1.0)
            y -= 1.0;
        y = 1.0 - y;
        auto baseSupport = simulationState.pm.get<scalar>("props.baseSupport");
        auto domainMin = simulationState.pm.get<vec>("domain.min");// - vec(baseSupport,baseSupport);
        auto domainMax = simulationState.pm.get<vec>("domain.max");// + vec(baseSupport,baseSupport);
        auto domainWidth = domainMax.x() - domainMin.x();
        auto domainHeight = domainMax.y() - domainMin.y();

        x *= domainWidth;
        y *= domainHeight;
        x += domainMin.x();
        y += domainMin.y();
        vec position(x, y);
        auto [ix, iy] = simulationState.getCellIdx(x, y);
        int32_t closestIdx = -1;
        scalar closestDist = DBL_MAX;
        bool found = false;
        auto cellsX = simulationState.pm.get<int32_t>("props.cellsX");
        auto cellsY = simulationState.pm.get<int32_t>("props.cellsY");
        for (int32_t xi = -1; xi <= 1; ++xi) {
            for (int32_t yi = -1; yi <= 1; ++yi) {
                if (ix + xi < 0 || ix + xi > cellsX - 1 || iy + yi < 0 || iy + yi > cellsY - 1)
                    continue;
                const auto& cell = simulationState.getCell(ix + xi, iy + yi);
                for (auto j : cell) {
                    scalar dist = (simulationState.fluidPosition[j] - position).norm();
                    if (dist < closestDist) {
                        closestDist = dist;
                        closestIdx = j;
                    }
                }
            }
        }
        if (closestIdx != -1)found = true;
        const float DISTANCE = 10.0f;
        static int corner = 3;
        ImGuiIO& io = ImGui::GetIO();
        if (corner != -1)
        {
            ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImVec2 work_area_pos = viewport->GetWorkPos();   // Instead of using viewport->Pos we use GetWorkPos() to avoid menu bars, if any!
            ImVec2 work_area_size = viewport->GetWorkSize();
            ImVec2 window_pos = ImVec2((corner & 1) ? (work_area_pos.x + work_area_size.x - DISTANCE) : (work_area_pos.x + DISTANCE), (corner & 2) ? (work_area_pos.y + work_area_size.y - DISTANCE) : (work_area_pos.y + DISTANCE));
            ImVec2 window_pos_pivot = ImVec2((corner & 1) ? 1.0f : 0.0f, (corner & 2) ? 1.0f : 0.0f);
            ImGui::SetNextWindowPos(window_pos, ImGuiCond_Always, window_pos_pivot);
            ImGui::SetNextWindowViewport(viewport->ID);
        }
        ImGui::SetNextWindowBgAlpha(0.35f); // Transparent background
        if (ImGui::Begin("Particle Info", &m_showInfo, (corner != -1 ? ImGuiWindowFlags_NoMove : 0) | ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_AlwaysAutoResize | ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing | ImGuiWindowFlags_NoNav))
        {
            auto addScalar = [&](auto name, scalar value, auto description) {
                float vcp = value;
                auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
                ImGui::DragFloat(name, &vcp, 0, vcp, vcp,"%.6f");
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
                if(ImGui::IsItemHovered())
                    ImGui::SetTooltip("%s",description);
                return;
            };
            struct float2{float x,y;};
            auto addScalar2 = [&](auto name, vec value, auto description) {
                float2 vcp{(float)value.x(), (float)value.y()};
                auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
                ImGui::DragFloat2(name, &vcp.x, 0, 0, 0, "%.6f");
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("%s",description);
                return;
            };
            auto addInteger = [&](auto name, int32_t value, auto description) {
                int32_t vcp = value;
                auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
                ImGui::DragInt(name, &vcp, 0, vcp, vcp);
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip("%s",description);
                return;
            };




            addScalar2("Cursor        ", position, "");
            addInteger("Index         ", closestIdx, "");
            int32_t i = closestIdx;
            if(found){
            addInteger("UID           ", simulationState.fluidUID[i], "");
            addScalar2("Position      ", simulationState.fluidPosition[i], "");
            addScalar2("Velocity      ", simulationState.fluidVelocity[i], "");
            addScalar2("Acceleration  ", simulationState.fluidAccel[i], "");
            addScalar ("Density       ", simulationState.fluidDensity[i], "");
            addScalar ("Neighbors     ", simulationState.fluidNeighborList[i].size(), "");
            //addScalar ("Vorticity     ", p.vorticity, "");
            //addScalar ("Angular Vel   ", p.angularVelocity, "");
            //addScalar ("dfsph: alpha  ", dp.alpha, "");
            //addScalar ("dfsph: area   ", dp.area, "");
            //addScalar ("dfsph: p1     ", dp.pressure1, "");
            //addScalar ("dfsph: p2     ", dp.pressure2, "");
            //addScalar ("dfsph: pb     ", dp.pressureBoundary, "");
            //addScalar ("dfsph: source ", dp.source, "");
            //addScalar ("dfsph: dpdt   ", dp.dpdt, "");
            //addScalar ("dfsph: rhostar", dp.rhoStar, "");
            //addScalar2("dfsph: vel    ", dp.vel, "");
            //addScalar2("dfsph: accel  ", dp.accel, "");
            }




        }
        ImGui::End();
    }
}
