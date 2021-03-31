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
                auto& param = ParameterManager::instance().getParameter(paramName);
                ImGui::PushID(param.identifier.c_str());
                ParameterManager::instance().uiFunctions[param.type](param);
                ImGui::PopID();
            };

            addParameter("time");
            addParameter("dt");
            //addParameter("simTime");
            //addParameter("renderTime");
            addParameter("numptcls");
            addParameter("densityIterations");
            addParameter("densityError");
            addParameter("divergenceIterations");
            addParameter("divergenceError");
            //addParameter("color_map.buffer");
            addParameter("colorMap.min");
            addParameter("colorMap.max");

        }
        ImGui::End();
    }
    if(m_showInfo){
        if (particles.size() == 0 ||
            ParameterManager::instance().get<int32_t>("sim.frame") == 0)return;


        auto x = std::clamp(m_cursorPosition.x(), 0.0, (scalar)screenWidth) / (scalar)screenWidth;
        auto y = 2.0 * std::clamp(m_cursorPosition.y(), 0.0, (scalar)screenHeight) / (scalar)screenHeight;
        if (y > 1.0)
            y -= 1.0;
        y = 1.0 - y;
        x *= domainWidth;
        y *= domainHeight;
        vec position(x, y);
        auto [ix, iy] = getCellIdx(x, y);
        int32_t closestIdx = -1;
        scalar closestDist = DBL_MAX;
        for (int32_t xi = -1; xi <= 1; ++xi) {
            for (int32_t yi = -1; yi <= 1; ++yi) {
                if (ix + xi < 0 || ix + xi > cellsX - 1 || iy + yi < 0 || iy + yi > cellsY - 1)
                    continue;
                const auto& cell = getCell(ix + xi, iy + yi);
                for (auto j : cell) {
                    scalar dist = (particles[j].pos - position).norm();
                    if (dist < closestDist) {
                        closestDist = dist;
                        closestIdx = j;
                    }
                }
            }
        }
        if (closestIdx == -1)return;
        const float DISTANCE = 10.0f;
        static int corner = 2;
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
                    ImGui::SetTooltip(description);
                return;
            };
            auto addScalar2 = [&](auto name, vec value, auto description) {
                float2 vcp(value.x(), value.y());
                auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
                ImGui::DragFloat2(name, &vcp.x(), 0, 0, 0, "%.6f");
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip(description);
                return;
            };
            auto addInteger = [&](auto name, int32_t value, auto description) {
                int32_t vcp = value;
                auto col = ImGui::GetStyle().Colors[ImGuiCol_FrameBg];
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = ImVec4(0.8f, 0.8f, 0.8f, 1.0f);
                ImGui::DragInt(name, &vcp, 0, vcp, vcp);
                ImGui::GetStyle().Colors[ImGuiCol_FrameBg] = col;
                if (ImGui::IsItemHovered())
                    ImGui::SetTooltip(description);
                return;
            };



            auto p = particles[closestIdx];
            auto dp = particlesDFSPH[closestIdx];

            addScalar2("Cursor        ", position, "");
            addInteger("Index         ", closestIdx, "");
            addInteger("UID           ", p.uid, "");
            addScalar2("Position      ", p.pos, "");
            addScalar2("Velocity      ", p.vel, "");
            addScalar2("Acceleration  ", p.accel, "");
            addScalar ("Density       ", p.rho, "");
            addScalar ("Neighbors     ", p.neighbors.size(), "");
            addScalar ("Vorticity     ", p.vorticity, "");
            addScalar ("Angular Vel   ", p.angularVelocity, "");
            addScalar ("dfsph: alpha  ", dp.alpha, "");
            addScalar ("dfsph: area   ", dp.area, "");
            addScalar ("dfsph: p1     ", dp.pressure1, "");
            addScalar ("dfsph: p2     ", dp.pressure2, "");
            addScalar ("dfsph: pb     ", dp.pressureBoundary, "");
            addScalar ("dfsph: source ", dp.source, "");
            addScalar ("dfsph: dpdt   ", dp.dpdt, "");
            addScalar ("dfsph: rhostar", dp.rhoStar, "");
            addScalar2("dfsph: vel    ", dp.vel, "");
            addScalar2("dfsph: accel  ", dp.accel, "");





        }
        ImGui::End();
    }
}
