#include <glad/glad.h> 
#include "glui.h"
#include <tools/Timer.h>
#include <imgui/implot.h>
#include <simulation/SPHRender.h>

void GUI::RayWindow(bool* p_open) {
    if (!ParameterManager::instance().get<bool>("ray.render")) return;
    ImGui::SetNextWindowSizeConstraints(ImVec2(1920, 1080), ImVec2(FLT_MAX, FLT_MAX));
    if (!ImGui::Begin("Ray Window", p_open))
    {
        ImGui::End();
        return;
    }
    int32_t i = 0;

    static auto& orig = ParameterManager::instance().get<vec>("ray.origin");
    static auto& target = ParameterManager::instance().get<vec>("ray.target");
    static auto& fov = ParameterManager::instance().get<scalar>("ray.fov");
    static auto& res = ParameterManager::instance().get<int32_t>("ray.resolution");
    static auto& isoLevel = ParameterManager::instance().get<scalar>("ray.iso");
    static auto& gScale = ParameterManager::instance().get<scalar>("ray.gScaleExplicit");
    static auto& gScaleFine = ParameterManager::instance().get<scalar>("ray.gScaleFine");
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



    static std::vector<glm::vec2> intersectionsExplicit;
    static std::vector<glm::vec2> intersectionsExplicitFine;
    static std::vector<glm::vec2> intersectionsImplicit;
    static std::vector<float> distancesExplicit;
    static std::vector<float> distancesExplicitFine;
    static std::vector<float> distancesImplicit;
    static std::vector<float> angles;
    static vec sOrig, sTarget;
    static 	scalar sFov, sIso, sG, sH, sGF;
    static int32_t sRes, sFrame, sSteps;
    static auto& frame = ParameterManager::instance().get<int32_t>("sim.frame");
    if (sGF != gScaleFine || sSteps != subSteps || sOrig.x() != orig.x() || sTarget.x() != target.x() || sOrig.y() != orig.y() || sTarget.y() != target.y() || sFov != fov || sRes != res || sIso != isoLevel || sG != gScale || sH != hScale || sFrame != frame) {
        sOrig = ParameterManager::instance().get<vec>("ray.origin");
        sTarget = ParameterManager::instance().get<vec>("ray.target");
        sFov = ParameterManager::instance().get<scalar>("ray.fov");
        sRes = ParameterManager::instance().get<int32_t>("ray.resolution");
        sIso = ParameterManager::instance().get<scalar>("ray.iso");
        sG = ParameterManager::instance().get<scalar>("ray.gScaleExplicit");
        sGF = gScaleFine;
        sH = ParameterManager::instance().get<scalar>("ray.hScale");
        sFrame = frame;
        sSteps = subSteps;

        intersectionsImplicit.clear();
        intersectionsExplicit.clear();
        intersectionsExplicitFine.clear();
        distancesImplicit.clear();
        distancesExplicit.clear();
        distancesExplicitFine.clear();
        angles.clear();
        glm::vec2 direction = glm::vec2(dir.x(), dir.y());
        glm::vec2 origin = glm::vec2(orig.x(), orig.y());
        for (int32_t r = 0; r <= res; ++r) {
            float angle = -degToRad(fov / 2.0) + degToRad(fov) * (float)r / (float)res;
            glm::vec2 rDir{ cos(angle) * direction.x - sin(angle) * direction.y,sin(angle) * direction.x + cos(angle) * direction.y };
            angles.push_back(angle);
            auto [hit, pos] = intersectExplicit(origin, rDir, hScaled, gScaled);
            //std::cout << std::boolalpha << hit << " -> " << pos.x << " " << pos.y << std::endl;
            if (hit) intersectionsExplicit.push_back(pos);
            else {
                auto [hit, tmin, tmax] = rayIntersectAABB(origin, direction, glm::vec2(0, 0), glm::vec2(domainWidth, domainHeight));
                intersectionsExplicit.push_back(origin + tmin * direction);
            }

            auto [hitFine, posFine] = intersectExplicit(origin, rDir, hScaled, support / gScaleFine);
            //std::cout << std::boolalpha << hit << " -> " << pos.x << " " << pos.y << std::endl;
            if (hitFine) intersectionsExplicitFine.push_back(posFine);
            else {
                auto [hit, tmin, tmax] = rayIntersectAABB(origin, direction, glm::vec2(0, 0), glm::vec2(domainWidth, domainHeight));
                intersectionsExplicitFine.push_back(origin + tmin * direction);
            }
            
        }
        for (int32_t r = 0; r <= res; ++r) {
            float angle = -degToRad(fov / 2.0) + degToRad(fov) * (float)r / (float)res;
            glm::vec2 rDir{ cos(angle) * direction.x - sin(angle) * direction.y,sin(angle) * direction.x + cos(angle) * direction.y };
            auto [hit, pos] = intersectImplicit(origin, rDir, hScaled, gScaled);
            //std::cout << std::boolalpha << hit << " -> " << pos.x << " " << pos.y << std::endl;
            if (hit) intersectionsImplicit.push_back(pos);
            else {
                auto [hit, tmin, tmax] = rayIntersectAABB(origin, direction, glm::vec2(0, 0), glm::vec2(domainWidth, domainHeight));
                intersectionsImplicit.push_back(origin + tmin * direction);
            }
        }
//#pragma omp parallel for
//        for (int32_t r = 0; r <= res; ++r) {
//            float angle = -degToRad(fov / 2.0) + degToRad(fov) * (float)r / (float)res;
//            glm::vec2 rDir{ cos(angle) * direction.x - sin(angle) * direction.y,sin(angle) * direction.x + cos(angle) * direction.y };
//            auto [hit, pos] = intersectChebyshev(origin, rDir, hScaled, gScaled);
//            //std::cout << "Chebyshev: " << std::boolalpha << hit << " -> " << pos.x << " " << pos.y << std::endl;
//#pragma omp critical
//            if (hit) intersectionsCheb.push_back(pos);
//        }

        float minDist = FLT_MAX;
        float maxDist = -FLT_MAX;
        for (const auto& inter : intersectionsExplicit) {
            auto d = glm::length(inter - origin);
            minDist = std::min(minDist, d);
            maxDist = std::max(maxDist, d);
            distancesExplicit.push_back(d);
        }
        for (const auto& inter : intersectionsImplicit) {
            auto d = glm::length(inter - origin);
            minDist = std::min(minDist, d);
            maxDist = std::max(maxDist, d);
            distancesImplicit.push_back(d);
        }
        for (const auto& inter : intersectionsExplicitFine) {
            auto d = glm::length(inter - origin);
            minDist = std::min(minDist, d);
            maxDist = std::max(maxDist, d);
            distancesExplicitFine.push_back(d);
        }
        ImPlot::SetNextPlotLimits(-degToRad(fov / 2.0), degToRad(fov / 2.0), minDist - 1e-5f, maxDist +1e-5f, ImGuiCond_Always);
    }

    if (ImPlot::BeginPlot("Absolute plot", "x", "f(x)", ImVec2(1900,500))) {
        ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 7.5f);
        ImPlot::PlotLine("Implicit", angles.data(), distancesImplicit.data(), angles.size());
        ImPlot::PopStyleVar();
        ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 2.5f);
        ImPlot::PlotLine("Explicit Base", angles.data(), distancesExplicit.data(), angles.size());
        ImPlot::PlotLine("Explicit Fine", angles.data(), distancesExplicitFine.data(), angles.size());
        ImPlot::PopStyleVar();
        ImPlot::EndPlot();
    }
    std::vector<float> diffFine, diffBase;
    diffFine.resize(distancesExplicit.size());
    diffBase.resize(distancesExplicit.size());
    float minDist = FLT_MAX;
    float maxDist = -FLT_MAX;
    for (auto i = 0; i < distancesExplicit.size(); ++i) {
        diffFine[i] = distancesImplicit[i] - distancesExplicitFine[i];
        diffBase[i] = distancesImplicit[i] - distancesExplicit[i];
        minDist = std::min(minDist, diffFine[i]);
        maxDist = std::max(maxDist, diffFine[i]);
        minDist = std::min(minDist, diffBase[i]);
        maxDist = std::max(maxDist, diffBase[i]);
    }
    ImPlot::SetNextPlotLimits(-degToRad(fov / 2.0), degToRad(fov / 2.0), minDist - 1e-5f, maxDist + 1e-5f, ImGuiCond_Always);
    if (ImPlot::BeginPlot("Error Plot", "x", "f(x)", ImVec2(1900, 500))) {
        ImPlot::PushStyleVar(ImPlotStyleVar_LineWeight, 2.5f);
        ImPlot::PlotLine("Implicit - Fine", angles.data(), diffFine.data(), angles.size());
        ImPlot::PlotLine("Implicit - Base", angles.data(), diffBase.data(), angles.size());
        ImPlot::PopStyleVar();
        ImPlot::EndPlot();
    }
    //for (auto tptr : TimerManager::getTimers()) {
    //    auto& t = *tptr;
    //    ImGui::PushStyleColor(ImGuiCol_Header, hexToCol(t.getColor()));
    //    auto id = t.getDecriptor() + (t.getSamples().size() > 0 ? 
    //        ("\t" + std::to_string(std::floorf(*t.getSamples().rbegin() * 100.f) * 0.01f) + "ms") : std::string(""));
    //    if (ImGui::CollapsingHeader((id + "###" + (t.getDecriptor() + std::to_string(i++))).c_str())) {
    //        ImGui::PushID(i);
    //        static char buf1[512] = "";
    //        strcpy_s(buf1, 512, t.getDecriptor().c_str());
    //        ImGui::PushItemWidth(-100);
    //        ImGui::InputText("Identifier", buf1, 64);
    //        if (t.getSourceLocation() != "") {
    //            strcpy_s(buf1, 512, (t.getSourceLocation() + ":" + std::to_string(t.getLine())).c_str());
    //            ImGui::InputText("Location", buf1, 64);
    //        }
    //        ImGui::PopItemWidth();
    //        if (auto statsopt = t.getStats()) {
    //            auto stats = statsopt.value();
    //            static char buf[512] = "";
    //            sprintf_s(buf, 512, "%.7g", std::floorf(stats.min * 100.f) * 0.01f);
    //            ImGui::Text("min"); ImGui::SameLine();
    //            ImGui::Button(buf, ImVec2(66, 0)); ImGui::SameLine();
    //            ImGui::Text("max"); ImGui::SameLine();
    //            sprintf_s(buf, 512, "%.7g", std::floorf(stats.max * 100.f) * 0.01f);
    //            ImGui::Button(buf, ImVec2(66, 0)); ImGui::SameLine();
    //            ImGui::Text("avg"); ImGui::SameLine();
    //            sprintf_s(buf, 512, "%.7g", std::floorf(stats.avg * 100.f) * 0.01f);
    //            ImGui::Button(buf, ImVec2(66, 0)); ImGui::SameLine();
    //            ImGui::Text("med"); ImGui::SameLine();
    //            sprintf_s(buf, 512, "%.7g", std::floorf(stats.median * 100.f) * 0.01f);
    //            ImGui::Button(buf, ImVec2(66, 0)); ImGui::SameLine();
    //            ImGui::Text("dev"); ImGui::SameLine();
    //            sprintf_s(buf, 512, "%.7g", std::floorf(stats.stddev * 100.f) * 0.01f);
    //            ImGui::Button(buf, ImVec2(66, 0));

    //            ImGui::PushItemWidth(-100);
    //            ImGui::PlotLines("Timings", t.getSamples().data(), t.getSamples().size(), 0.f, NULL, stats.min, stats.max, ImVec2(0, 80));
    //            auto bmax = 0.f;
    //            for (auto& b : stats.histogram)
    //                bmax = std::max(b, bmax);
    //            ImGui::PlotHistogram("Histogram", stats.histogram.data(), stats.histogram.size(), 0, NULL, 0, bmax, ImVec2(0, 80));
    //            ImGui::PopItemWidth();
    //        }

    //        ImGui::PopID();
    //    }
    //    ImGui::PopStyleColor();
    //}
    ImGui::End();
}
