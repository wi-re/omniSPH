#include <tools/exceptionHandling.h>
#include "gui/glui.h"
#include <sstream>
#include <thread>
#include <windows.h>
#include <iostream> 

BOOL CtrlHandler(DWORD fdwCtrlType) {
    std::clog << "Caught signal " << fdwCtrlType << std::endl;
    switch (fdwCtrlType)
    {
    case CTRL_CLOSE_EVENT:
        GUI::instance().quit();
        GUI::instance().render_lock.unlock();
        std::this_thread::sleep_for(std::chrono::milliseconds(2000));
        return(TRUE);

    default:
        return FALSE;
    }
}
#include <test.h>
#include <simulation/2DMath.h>
#include <cheb/cheb.h>
#include <iomanip>


scalar sdpolygon2(std::vector<vec> v, vec p) {
    scalar d = (p - v[0]).dot(p - v[0]);
    scalar s = 1.0;
    for (int i = 0, j = v.size() - 1; i < v.size(); j = i, i++) {
        vec e = v[j] - v[i];
        vec w = p - v[i];
        vec b = w - e * std::clamp(w.dot(e) / e.dot(e), 0.0, 1.0);
        d = std::min(d, b.dot(b));
        bool b1 = p.y() >= v[i].y();
        bool b2 = p.y() < v[j].y();
        bool b3 = e.x() * w.y() > e.y() * w.x();
        if ((b1 && b2 && b3) || (!b1 && !b2 && !b3))
            s *= -1.0;
    }
    return s * sqrt(d);
}

int main(int argc, char* argv[]) 
try 
{
    auto getDelta = [](scalar packing) {
        auto cubeMax = [](auto x) {
            return std::max(0.0, x * x * x);
        };
        cheb::scalar delta = 0.0;
        for (int32_t i = -5; i <= 5; ++i) {
            for (int32_t j = -5; j <= 5; ++j) {
                cheb::scalar x = packing * (2.0 * (cheb::scalar)i + (cheb::scalar)(j % 2));
                cheb::scalar y = packing * std::sqrt(3.0) * cheb::scalar(j);
                cheb::scalar dist = std::sqrt(x * x + y * y);
                cheb::scalar q = dist / std::sqrt(20.0 / M_PI);
                //std::cout << i << " : " << j << " -> " << x << " : " << y << " => " << dist << " : " << q << " ==> " << 4. / 7. * (cubeMax(1 - q) - 4. * cubeMax(0.5 - q)) << std::endl;
                delta += 4. / 7. * (cubeMax(1 - q) - 4. * cubeMax(0.5 - q));
            }
        }
        return delta;
    };
    auto getDelta2 = [](scalar packing) {
        auto cubeMax = [](auto x) {
            return std::max(0.0, x * x * x);
        };
        cheb::scalar delta = 0.0;
        for (int32_t i = -5; i <= 5; ++i) {
            for (int32_t j = 0; j <= 5; ++j) {
                cheb::scalar x = packing * (2.0 * (cheb::scalar)i + (cheb::scalar)(j % 2));
                cheb::scalar y = packing * std::sqrt(3.0) * cheb::scalar(j);
                cheb::scalar dist = std::sqrt(x * x + y * y);
                cheb::scalar q = dist / std::sqrt(20.0 / M_PI);
                //std::cout << i << " : " << j << " -> " << x << " : " << y << " => " << dist << " : " << q << " ==> " << 4. / 7. * (cubeMax(1 - q) - 4. * cubeMax(0.5 - q)) << std::endl;
                delta += 4. / 7. * (cubeMax(1 - q) - 4. * cubeMax(0.5 - q));
            }
        }
        return delta;
    };


    auto desiredCompression = 1.00;
    auto desiredParticles = 100.0;

    cheb::Function packFun([&](cheb::scalar s) {
        return getDelta(s * std::sqrt(20) / std::sqrt(M_PI)) - desiredCompression;
        }, cheb::Domain{ 0.01,0.5 });
    auto chebPacking = packFun.roots()[0];

    cheb::Function distFun([&](cheb::scalar s) {
        return getDelta2(chebPacking * std::sqrt(20) / std::sqrt(M_PI)) + boundaryKernel(-s) - desiredCompression;
        }, cheb::Domain{ 0.01,0.5 });
    auto chebDist = distFun.roots()[0];
    std::cout << "Packing: " << std::setprecision(17) << chebPacking << std::endl;
    std::cout << "Spacing: " << std::setprecision(17) << chebDist << std::endl;

    cheb::Function spaceFun([&](cheb::scalar s) {
        return ((domainHeight - 2.0 * domainEpsilon) / s - 2.0 * chebDist) / ( std::sqrt(3) * chebPacking) - desiredParticles;
        }, cheb::Domain{ 0.01,1.0 });
    std::cout << "Scaling: " << std::setprecision(17) << spaceFun.roots()[0] << std::endl;

    inletPacking = chebPacking;
    inletSpacing = chebDist;
    std::cout << "Error: " << getDelta(packing_2D * std::sqrt(20) / std::sqrt(M_PI)) << std::endl;

    cheb::scalar delta = 0.0;
    for (int32_t i = -5; i <= 5; ++i) {
        for (int32_t j = -5; j <= 5; ++j) {
            cheb::scalar x = packing_2D * (2.0 * (cheb::scalar)i + (cheb::scalar)(j % 2));
            cheb::scalar y = packing_2D * std::sqrt(3.0) * cheb::scalar(j);
            cheb::scalar dist = std::sqrt(x * x + y * y);
            cheb::scalar q = dist / support;
            //std::cout << i << " : " << j << " -> " << x << " : " << y << " => " << dist << " : " << q << " ==> " << 4. / 7. * (cubeMax(1 - q) - 4. * cubeMax(0.5 - q)) << std::endl;
            delta += area * W(vec(0, 0), vec(x, y));
        }
    }
    std::cout << "Delta: " << delta << std::endl;
    std::cout << "\n\n\n\n";
    std::cout << "constexpr inline scalar scale = " << std::setprecision(17) << spaceFun.roots()[0] << "; // desired particles: " << desiredParticles << std::endl;
    std::cout << "constexpr inline scalar packing_2D = " << std::setprecision(17) << chebPacking << " * scale; // desired compression: " << desiredCompression << std::endl;
    std::cout << "constexpr inline scalar spacing_2D = " << std::setprecision(17) << chebDist << " * scale; // actual delta: " << delta << std::endl;
    std::cout << "\n\n\n\n";




    vec P1(domainEpsilon, domainEpsilon);
    vec P2(domainWidth - domainEpsilon, domainEpsilon);
    vec P3(domainWidth - domainEpsilon, domainHeight - domainEpsilon);
    vec P4(domainEpsilon, domainHeight - domainEpsilon);
    std::vector<vec> polygon;
    polygon.push_back(P2);
    polygon.push_back(P3);
    polygon.push_back(P4);
    polygon.push_back(P1);
    auto v = vec(2 * domainEpsilon, domainEpsilon + spacing_2D);
    scalar d =  sdpolygon2(polygon, v);
    std::cout << d << " -> " << d / scale << std::endl;
    auto [hit, pb, dis, k, gk] = interactLines(v, polygon, true);
    std::cout << hit << " -> " << dis << ", " << k << " @ " << dis / scale << " -> " << boundaryKernel(dis/scale) << std::endl;


    std::cout << "r = " << 1 / std::sqrt(M_PI) << std::endl;
    std::cout << "h = " << std::sqrt(20) / std::sqrt(M_PI) << std::endl;
    std::cout << "packing = " << 0.21314955649168332 * std::sqrt(20) / std::sqrt(M_PI) << std::endl;

    std::cout << getDelta(0.21314955649168332 * std::sqrt(20) / std::sqrt(M_PI)) << std::endl;

    // found by maxima
    auto h = std::sqrt(20) / std::sqrt(M_PI);
    auto dist = 0.1977501895407778 * h; 
    auto packing = 0.2131495564916834 * h;

    std::cout << getDelta2(packing) << std::endl;
    std::cout << boundaryKernel(-dist / h) << " : " << boundaryKernel(dist / h) << std::endl;
    std::cout << getDelta2(packing) + boundaryKernel(-dist / h) << std::endl;

    //return 0;





    //cheb::Function testFn([](cheb::scalar s) {
    //    if (s < 0) return -1.;
    //    return 1.;
    //    });
    //testInterval();
    //testDomain();
    //testChebTech();
    //testChebTech_functions();
    //testArithmetic();
    //testRoots();
    //testUtilities();

    //testFunction_Construction();
    //testFunction_Properties();
    //testFunction_ClassUsage();
    //testFunction_Evaluation();
    //testFunction_Calculus();
    //testFunction_Roots();
    //testFunction_Arithmetic();
    //testFunction();
    if (!SetConsoleCtrlHandler((PHANDLER_ROUTINE)CtrlHandler, TRUE)) {
        std::cerr << "Could not set CtrlHandler. Exiting." << std::endl;
        return 0;
    }
    //_set_se_translator(translator);

    auto& gui = GUI::instance();
    gui.render_lock.lock();
    gui.initParameters(argc, argv);
    gui.initGVDB();
    gui.initSimulation();
    gui.initGL();
    gui.renderLoop();
    return 0;
}
CATCH_DEFAULT