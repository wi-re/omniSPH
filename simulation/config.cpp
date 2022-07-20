#include <simulation/SPH.h>
#include <config/config.h>
#include <yaml-cpp/yaml.h>
#include "2DMath.h"
#include <imgui/imgui.h>


namespace YAML {
template<>
struct convert<vec> {
  static Node encode(const vec& rhs) {
    Node node;
    node.push_back(rhs.x());
    node.push_back(rhs.y());
    return node;
  }

  static bool decode(const Node& node, vec& rhs) {
    if(!node.IsSequence() || node.size() != 2) {
      return false;
    }

    rhs.x() = node[0].as<double>();
    rhs.y() = node[1].as<double>();
    return true;
  }
};
}

void SPHSimulation::initializeParameters(){
  auto working_directory = std::filesystem::current_path(); 

  pm.newParameter("internal.working_directory", std::string(working_directory.string()), {.constant = true});
  pm.newParameter("internal.source_directory", std::string(sourceDirectory), {.constant = true});
  pm.newParameter("internal.build_directory", std::string(binaryDirectory), {.constant = true});

  std::string fileName = "cfg/viridis.png";
  std::string stylePath = resolveFile(fileName).string();
  pm.newParameter("stylePath", stylePath, {.hidden = true});

    
    static auto vectorModes = std::vector<detail::iAny>{
            std::string("magnitude"),
            std::string("x"),
            std::string("y") };
    static auto buffers = std::vector<detail::iAny>{        
        std::string("position"),
        std::string("velocity"),
        std::string("accel"),
        std::string("predVelocity"),
        std::string("predAccel"),
        std::string("predPosition"),

        std::string("density"),
        std::string("vorticity"),
        std::string("angularVelocity"),
        std::string("area"),
        std::string("restDensity"),
        std::string("support"),

        std::string("alpha"),
        std::string("actualArea"),
        std::string("pressure1"),
        std::string("pressure2"),
        std::string("boundaryPressure"),
        std::string("sourceTerm"),
        std::string("dpdt"),
        std::string("rhoStar"),
        std::string("priorPressure"),

        std::string("UID"),
        std::string("neighbors")
    };

    static auto cMapPresets = std::vector<std::string>{[this]() {
    std::vector<std::string> colorMaps;
    auto f = std::filesystem::path(pm.get<std::string>("stylePath"));
    auto p = f.parent_path().string();
    if (*(p.end() - 1) == '/' || *(p.end() - 1) == '\\')
      p = p.substr(0, p.length() - 1);
    std::replace(p.begin(), p.end(), '\\', '/');
    for (auto &p : std::filesystem::directory_iterator(p))
      if (p.path().extension().string().find(".png") != std::string::npos)
        colorMaps.push_back(p.path().filename().replace_extension("").string());
    return colorMaps;
  }()};
  static auto colorMaps = std::vector<detail::iAny>{};
  for (auto c : cMapPresets)
    colorMaps.push_back(std::string(c));


    pm.newParameter("colorMap.vectorMode", std::string("magnitude"), { .constant = false
        });
    pm.newParameter("colorMap.buffer", std::string("neighbors"), { .constant = false, .presets = buffers
        });


    auto genRange = [](auto min, auto max) { return Range{ min, max }; };

    pm.newParameter("sim.frame", 0, { .constant = false });

    pm.newParameter("colorMap.min", scalar(0.99), { .constant = false , .range = genRange(-10.0,10.0)});
    pm.newParameter("colorMap.cmap", std::string("viridis"), { .constant = false , .presets = colorMaps});
    pm.newParameter("colorMap.map", 0, { .constant = false , .range = genRange(0,3)});
    pm.newParameter("colorMap.limit", false, { .constant = false });
    pm.newParameter("colorMap.flip", false, { .constant = false });
    pm.newParameter("colorMap.auto", true, { .constant = false });

    pm.newParameter("colorMap.max", scalar(1.03), { .constant = false , .range = genRange(-10.0,10.0)});

    pm.newParameter("vorticity.nu_t", 0.025, { .constant = false, .range = genRange(0.0, 1.0)});
    pm.newParameter("vorticity.inverseInertia", 0.5, { .constant = false, .range = genRange(0.0, 1.0)});
    pm.newParameter("vorticity.angularViscosity", 0.005, { .constant = false, .range = genRange(0.0, 1.0)});

    pm.newParameter("sim.time", 0., { .constant = true });
    pm.newParameter("props.numPtcls", 0, { .constant = true });
    pm.newParameter("props.maxnumptcls", 1024*128, { .constant = true });
    //pm.newParameter("props.scale", simulationState.baseSupport, { .constant = true });
    pm.newParameter("ptcl.render", true, { .constant = false });
    pm.newParameter("props.baseRadius", 0., { .constant = true });
    pm.newParameter("props.baseSupport", 0., { .constant = true });
    pm.newParameter("props.cellsX", 0, { .constant = true });
    pm.newParameter("props.cellsY", 0, { .constant = true });

    pm.newParameter("export.active", false, { .constant = true });
    pm.newParameter("export.interval", 1, { .constant = true });
    pm.newParameter("export.limit", -1., { .constant = true });
    pm.newParameter("video.active", false, { .constant = true });
    pm.newParameter("video.fps", -1., { .constant = true });

    //pm.newParameter("ptcl.radius", radius, { .constant = true });
    pm.newParameter("sim.minDt", 0.001, { .constant = false, .range = genRange(0.0001, 0.016)});
    pm.newParameter("sim.maxDt", 0.002, { .constant = false, .range = genRange(0.0001, 0.016)});
    pm.newParameter("sim.barycentricPressure", false, { .constant = false});
    pm.newParameter("sim.dt", 0.001, { .constant = true });
    pm.newParameter("props.backgroundPressure", false, { .constant = true });
    pm.newParameter("domain.min", vec(-1.,-1.), { .constant = true });
    pm.newParameter("domain.max", vec(1.,1.), { .constant = true });
    pm.newParameter("domain.epsilon", 0.02, {.constant=true});
    //pm.newParameter("props.cellsX", cellsX, { .constant = true });
    //pm.newParameter("props.cellsY", cellsY, { .constant = true });
    //pm.newParameter("props.packing_2D", packing_2D, { .constant = true });
    //pm.newParameter("props.support", support, { .constant = true });
    pm.newParameter("props.damping", 0.0, { .constant = false, .range = genRange(0.00,1.00)});
    //pm.newParameter("ptcl.area", area, { .constant = true });
    //pm.newParameter("ptcl.mass", mass, { .constant = true });
    pm.newParameter("ptcl.maxVelocity", scalar(0), { .constant = true });
    pm.newParameter("ptcl.viscosityConstant", 0.0001, { .constant = false, .range = genRange(0.01, 0.05)});
    pm.newParameter("ptcl.boundaryViscosity", 0.50, { .constant = false, .range = genRange(0.01, 1.0)});

    pm.newParameter("state.kappa", 2.16, { .constant = false, .range = genRange(0.1, 10.0)});
    pm.newParameter("state.gamma", 7.15, { .constant = false, .range = genRange(0.1, 100.0)});
    pm.newParameter("state.restPressure", 101325., { .constant = false, .range = genRange(0., 1e6)});
    pm.newParameter("state.useEOS", false, { .constant = false});

    pm.newParameter("dfsph.densityEta", scalar(0.001), { .constant = true });
    pm.newParameter("dfsph.divergenceEta", scalar(0.001), { .constant = true });
    pm.newParameter("dfsph.densityIterations", 0, { .constant = true });
    pm.newParameter("dfsph.divergenceIterations", 0, { .constant = true });
    pm.newParameter("dfsph.densityError", scalar(0), { .constant = true });
    pm.newParameter("dfsph.divergenceError", scalar(0), { .constant = true });
    pm.newParameter("dfsph.divergenceSolve", true, { .constant = false });

    pm.get<scalar>("ptcl.boundaryViscosity") = 0.;
    pm.get<scalar>("ptcl.viscosityConstant") = 0.001;
    pm.get<scalar>("sim.maxDt") = 0.001;
    pm.get<scalar>("sim.minDt") = 0.000125;
    pm.get<scalar>("vorticity.nu_t") = 0.0125;
    pm.get<scalar>("vorticity.angularViscosity") = 0.1;
    pm.get<bool>("dfsph.divergenceSolve") = true;




#define callui(member) callUI(std::decay_t<decltype(val.member)>, val.member, #member)
#define encode(ty, x) encoders[typeid(ty)](x);
    #define vectorEncoder(ty, encoder)\
	pm.addEncoder(typeid(std::vector<ty>), [this](const detail::iAny& any) {\
		auto var = boost::any_cast<std::vector<ty>>(any);\
		auto node = YAML::Node();\
		for (const auto& elem : var){\
            node.push_back(encoder(elem));\
		}\
		return node;\
		});\
	pm.addDecoder(typeid(std::vector<ty>), [this](const YAML::Node& node) {\
		std::vector<ty> vec;\
		for (const auto& child : node)\
			vec.push_back((ty)pm.decoders[typeid(ty)](child));\
		return detail::iAny(vec);\
		}); 


    static auto encoders = pm.encoders;

    static auto fluidSourceEncoder = [](const fluidSource& var){
	    auto node  = YAML::Node();
        node["min"] = encode(vec, var.emitterMin);
        node["max"] = encode(vec, var.emitterMax);
        node["density"] = encode(scalar, var.emitterDensity);
        node["velocity"] = encode(vec, var.emitterVelocity);
        node["timeLimit"] = encode(scalar, var.timeLimit);
        node["ramp"] = encode(scalar, var.emitterRampTime);
        node["radius"] = encode(scalar, var.emitterRadius);
        node["shape"] = var.shape == shape_t::rectangular ? "rectangular" : "spherical";
        switch(var.emitter){
            case emitter_t::inlet: node["type"] = "inlet"; break;
            case emitter_t::outlet: node["type"] = "outlet"; break;
            case emitter_t::oneTime: node["type"] = "once"; break;
            case emitter_t::velocitySource: node["type"] = "source"; break;
        }
        return node;
    };
    pm.addDecoder(typeid(fluidSource), [this](const YAML::Node& node){
        fluidSource var;
        var.emitterMin = callDecoder(vec, node["min"], var.emitterMin);
        var.emitterMax = callDecoder(vec, node["max"], var.emitterMax);
        var.emitterDensity = callDecoder(scalar, node["density"], var.emitterDensity);
        var.emitterVelocity = callDecoder(vec, node["velocity"], var.emitterVelocity);
        var.timeLimit = callDecoder(scalar, node["timeLimit"], var.timeLimit);
        var.emitterRampTime = callDecoder(scalar, node["ramp"], var.emitterRampTime);
        var.emitterRadius = callDecoder(scalar, node["radius"], var.emitterRadius);        
        var.shape = shape_t::rectangular;
        if(node["shape"]){
            std::string s = node["shape"].as<std::string>();
            if(s == "rectangular")
                var.shape = shape_t::rectangular;
            if(s == "spherical")
                var.shape = shape_t::spherical;
        }
        var.emitter = emitter_t::oneTime;
        if(node["type"]){
            std::string s = node["type"].as<std::string>();
            if(s == "once")
                var.emitter = emitter_t::oneTime;
            if(s == "inlet")
                var.emitter = emitter_t::inlet;
            if(s == "outlet")
                var.emitter = emitter_t::outlet;
            if(s == "source")
                var.emitter = emitter_t::velocitySource;
        }
        return detail::iAny(var);
    });
    vectorEncoder(fluidSource, fluidSourceEncoder);

    static auto triangleEncoder = [](const Triangle& var){
	    auto node  = YAML::Node();
        node["v0"] = encode(vec, var.v0);
        node["v1"] = encode(vec, var.v1);
        node["v2"] = encode(vec, var.v2);
        return node;
    };
    pm.addDecoder(typeid(Triangle), [this](const YAML::Node& node){
        Triangle var;
        var.v0 = callDecoder(vec, node["v0"], var.v0);
        var.v1 = callDecoder(vec, node["v1"], var.v1);
        var.v2 = callDecoder(vec, node["v2"], var.v2);
        return detail::iAny(var);
    });
    vectorEncoder(Triangle, triangleEncoder);


    static auto gravityEncoder = [](const gravitySource& var){
	    auto node  = YAML::Node();
        node["pointSource"] = encode(bool, var.pointSource);
        node["direction"] = encode(vec, var.direction);
        node["location"] = encode(vec, var.location);
        node["magnitude"] = encode(scalar, var.magnitude);
        return node;
    };
    pm.addDecoder(typeid(gravitySource), [this](const YAML::Node& node){
        gravitySource var;
        var.pointSource = callDecoder(bool, node["pointSource"], var.pointSource);
        var.direction = callDecoder(vec, node["direction"], var.direction);
        var.location = callDecoder(vec, node["location"], var.location);
        var.magnitude = callDecoder(scalar, node["magnitude"], var.magnitude);
        return detail::iAny(var);
    });
    vectorEncoder(gravitySource, gravityEncoder);


    pm.newParameter("fluids", std::vector<fluidSource>{},{.constant=false});
    pm.newParameter("triangles", std::vector<Triangle>{}, {.constant=true});
    pm.newParameter("gravity", std::vector<gravitySource>{},{.constant=false});

    pm.loadTree(config);

    if(pm.get<scalar>("video.fps") < 0.){
        pm.get<scalar>("video.fps") = 1. / pm.get<scalar>("sim.minDt");
    }

    //::initializeParameters();
}


void SPHSimulation::initializeSPH(){
    std::size_t n = pm.get<int32_t>("props.maxnumptcls");
    //baseSupport = 0.012512272327;
    //baseRadius = 0.22360679774997896 * baseSupport;

    fluidSources = pm.get<std::vector<fluidSource>>("fluids");
    if(config["gravity"]){
    }else{
        pm.get<std::vector<gravitySource>>("gravity").push_back(gravitySource{
            .pointSource = false,
            .direction = vec(0,-1.),
            .magnitude = 9.81
        });
    }
    gravitySources = pm.get<std::vector<gravitySource>>("gravity");


    // config["fluids"];
    // for(int32_t i = 0; i < fluidVolumes.size(); ++i){
    //     auto vol = fluidVolumes[i];
    //     // std::cout << fluidVolumes[i] << std::endl;
    //     fluidSource fs;
    //     // std::cout << vol["min"] << std::endl;
    //     fs.emitterMin = vol["min"].as<vec>();
    //     fs.emitterMax = vol["max"].as<vec>();
    //     fs.emitterDensity = vol["density"] ? vol["density"].as<double>() : 998.0;
    //     fs.emitterVelocity = vol["velocity"] ? vol["velocity"].as<vec>() : vec(0.,0.);
    //     fs.timeLimit = vol["timeLimit"] ? vol["timeLimit"].as<double>() : -1.;
    //     fs.emitterRampTime = vol["ramp"] ? vol["ramp"].as<double>() : -1.;
    //     fs.emitterRadius = vol["radius"].as<scalar>();
    //     fs.shape = shape_t::rectangular;
    //     if(vol["shape"]){
    //         std::string s = vol["shape"].as<std::string>();
    //         if(s == "rectangular")
    //             fs.shape = shape_t::rectangular;
    //         if(s == "spherical")
    //             fs.shape = shape_t::spherical;
    //     }
    //     fs.emitter = emitter_t::oneTime;
    //     if(vol["type"]){
    //         std::string s = vol["type"].as<std::string>();
    //         if(s == "once")
    //             fs.emitter = emitter_t::oneTime;
    //         if(s == "inlet")
    //             fs.emitter = emitter_t::inlet;
    //         if(s == "outlet")
    //             fs.emitter = emitter_t::outlet;
    //         if(s == "source")
    //             fs.emitter = emitter_t::velocitySource;
    //     }
    //     fluidSources.push_back(fs);
    // }
    pm.get<scalar>("props.baseRadius") = 0.;
    for(auto& src : fluidSources){
        pm.get<scalar>("props.baseRadius") = std::max(pm.get<scalar>("props.baseRadius"), src.emitterRadius);
    }
    pm.get<scalar>("props.baseSupport") = std::sqrt(pm.get<scalar>("props.baseRadius") * pm.get<scalar>("props.baseRadius") * targetNeighbors);




    fluidPosition.resize(n);
    fluidVelocity.resize(n);
    fluidAccel.resize(n);
    fluidPredVelocity.resize(n);
    fluidPredAccel.resize(n);
    fluidPredPosition.resize(n);
    fluidDensity.resize(n);
    fluidVorticity.resize(n);
    fluidAngularVelocity.resize(n);
    fluidArea.resize(n);
    fluidRestDensity.resize(n);
    fluidSupport.resize(n);
    fluidAlpha.resize(n);
    fluidActualArea.resize(n);
    fluidPressure1.resize(n);
    fluidPressure2.resize(n);
    fluidBoundaryPressure.resize(n);
    fluidSourceTerm.resize(n);
    fluidDpDt.resize(n);
    fluidDensityStar.resize(n);
    fluidPriorPressure.resize(n);
    fluidNeighborList.resize(n);
    fluidUID.resize(n);
    fluidTriangleNeighborList.resize(n);

    scalar baseArea = std::numbers::pi * pm.get<scalar>("props.baseRadius") * pm.get<scalar>("props.baseRadius");
    baseRadius = pm.get<scalar>("props.baseRadius");
    baseSupport = pm.get<scalar>("props.baseSupport");

    domainMin = pm.get<vec>("domain.min");
    domainMax = pm.get<vec>("domain.max");
    auto domainWidth = domainMax.x() - domainMin.x();
    auto domainHeight = domainMax.y() - domainMin.y();


    pm.get<int32_t>("props.cellsX") = std::size_t(std::ceil(domainWidth / baseSupport));
    pm.get<int32_t>("props.cellsY") = std::size_t(std::ceil(domainHeight / baseSupport));
    cellsX = pm.get<int32_t>("props.cellsX");
    cellsY = pm.get<int32_t>("props.cellsY");
    for(int32_t i = 0; i < pm.get<int32_t>("props.cellsX") * pm.get<int32_t>("props.cellsY"); ++i){
        std::vector<int32_t> cells;
        cellArray.push_back(cells);
        cellBoundaryArray.push_back(cells);
        cellTriangleArray.push_back(cells);
    }

    auto& boundaryTriangles = pm.get<std::vector<Triangle>>("triangles");

    //gravity = vec(0.0,0.0);
    float d = 1.5;
    auto addRect = [&](scalar xmin, scalar ymin, scalar xmax, scalar ymax) {
        boundaryTriangles.push_back(Triangle{ {xmin,ymin},{xmax,ymax},{xmin,ymax} });
        boundaryTriangles.push_back(Triangle{ {xmin,ymin},{xmax,ymin},{xmax,ymax} });

    };
    auto addRectSub = [&](scalar xmin, scalar ymin, scalar xmax, scalar ymax, int32_t nx, int32_t ny) {
        if (nx == ny && (nx == 1 || ny == 2)) {
            boundaryTriangles.push_back(Triangle{ {xmin,ymin},{xmax,ymax},{xmin,ymax} });
            boundaryTriangles.push_back(Triangle{ {xmin,ymin},{xmax,ymin},{xmax,ymax} });
            return;
        }
        auto dx = (xmax - xmin) / (scalar)(nx - 1);
        auto dy = (ymax - ymin) / (scalar)(ny - 1);
        for (int32_t i = 0; i < nx - 1; ++i) {
            auto xn = xmin + dx * (scalar)(i + 0);
            auto xp = xmin + dx * (scalar)(i + 1);
            for (int32_t j = 0; j < ny - 1; ++j) {
                auto yn = ymin + dy * (scalar)(j + 0);
                auto yp = ymin + dy * (scalar)(j + 1);
                addRect(xn, yn, xp, yp);
            }
        }

    };
    auto addRectSubh = [&](scalar xmin, scalar ymin, scalar xmax, scalar ymax, scalar lim) {
        auto nx = (int32_t) ::ceil((xmax - xmin) / lim);
        auto ny = (int32_t) ::ceil((ymax - ymin) / lim);
        addRectSub(xmin, ymin, xmax, ymax, nx + 1, ny+1);
    };


    // if(config["triangles"]){
    //     auto triangles = config["triangles"];
    //     for(int32_t i = 0; i < triangles.size(); ++i){
    //         auto node = triangles[i];
    //         boundaryTriangles.push_back(Triangle{
    //             .v0 = node["v0"].as<vec>(),
    //             .v1 = node["v1"].as<vec>(),
    //             .v2 = node["v2"].as<vec>()
    //         });
    //     }
    // }

        auto domainEpsilon = pm.get<scalar>("domain.epsilon");
        addRect(domainMin.x() + 0, 
                domainMin.y() + 0, 
                domainMin.x() + domainEpsilon, 
                domainMin.y() + domainHeight);
        addRect(domainMin.x() + domainEpsilon, 
                domainMin.y() + 0, 
                domainMin.x() + domainWidth - domainEpsilon, 
                domainMin.y() + domainEpsilon);
        addRect(domainMin.x() + domainWidth - domainEpsilon, 
                domainMin.y() + 0, 
                domainMin.x() + domainWidth, 
                domainMin.y() + domainHeight);
        addRect(domainMin.x() + domainEpsilon,
                domainMin.y() + domainHeight - domainEpsilon, 
                domainMin.x() + domainWidth - domainEpsilon, 
                domainMin.y() + domainHeight);


    
        scalar threshold = 0.39960767069134208 * baseSupport * 0.5;
        auto checkPosition = [&](vec pos){
            auto [ix, iy] = getCellIdx(pos.x(), pos.y());
            for (int32_t xi = -1; xi <= 1; ++xi) {
                if(xi + ix < 0) continue;
                for (int32_t yi = -1; yi <= 1; ++yi) {
                    if(yi + iy < 0) continue;
                    for(int32_t idx : getBoundaryCell(ix + xi, iy + yi)){
                        auto [pb,b] = boundaryParticles[idx];
                        auto d = (pb - pos).norm();
                        if(d <= threshold)
                            return false;
                    }
                }
            }
            return true;
        };
        auto checkedEmit = [this, checkPosition](vec pos){
            if(checkPosition(pos)){
                boundaryParticles.push_back(std::make_pair(pos,vec(0,0)));
                auto idx = boundaryParticles.size() - 1;
                auto [x,y] = getCellIdx(pos.x(), pos.y());
                getBoundaryCell(x,y).push_back(idx);
            }
        };
        auto triangleSign = [](vec v1, vec v2, vec v3){
            return (v1.x() - v3.x()) * (v2.y() - v3.y()) - (v2.x() - v3.x()) * (v1.y() - v3.y());
        };
        auto pointInTriangle = [triangleSign](vec pt, vec v1, vec v2, vec v3){
            auto d1 = triangleSign(pt, v1, v2);
            auto d2 = triangleSign(pt, v2, v3);
            auto d3 = triangleSign(pt, v3, v1);

            bool has_neg = (d1 < 0) || (d2 < 0) || (d3 < 0);
            bool has_pos = (d1 > 0) || (d2 > 0) || (d3 > 0);

            return !(has_neg && has_pos);
        };

        std::mutex m;

        int32_t xSlices = domainWidth / 0.39960767069134208 / baseSupport;
        int32_t ySlices = domainHeight / 0.39960767069134208 / baseSupport;
#pragma omp parallel for
        for(int32_t i = 0; i <= xSlices; ++i){
            for(int32_t j = 0; j <= ySlices; ++ j){
                auto x = ((double)i )/((double)xSlices) * domainWidth + domainMin.x();
                auto y = ((double)j )/((double)ySlices) * domainHeight + domainMin.y();
                vec p(x,y);
                for(int32_t ti = 0; ti < boundaryTriangles.size(); ++ti){
                    const auto& t = boundaryTriangles[ti];
                    auto [v0,v1,v2] = t;
                    if(pointInTriangle(p,v0,v1,v2))
                    {
                        std::lock_guard lg(m);
                        boundaryParticles.push_back(std::make_pair(p,vec(0,0)));
                        break;
                    }
                }
            }
        }
        for(int32_t i = 0; i < boundaryParticles.size(); ++i){
            auto [pos,normal] = boundaryParticles[i];
            auto [x,y] = getCellIdx(pos.x(), pos.y());
            getBoundaryCell(x,y).push_back(i);
        }
        for(int32_t i = 0; i < boundaryParticles.size(); ++i){
            auto [pos,normal] = boundaryParticles[i];
            
            auto [ix, iy] = getCellIdx(pos.x(), pos.y());
            auto d = 4;
            for (int32_t xi = -d; xi <= d; ++xi) {
                if(xi + ix < 0) continue;
                if(xi + ix > cellsX - 1) continue;
                for (int32_t yi = -d; yi <= d; ++yi) {
                    if(yi + iy < 0) continue;
                    if(yi + iy > cellsY - 1) continue;

                    for(int32_t idx : getBoundaryCell(ix + xi, iy + yi)){
                        auto [pb,nb] = boundaryParticles[idx];
                        // auto d = (pb - pos).norm();
                        // auto q = d / support;
                        
                        vec a = pos;
                        vec b = pb;
                        vec r = a - b;
                        auto rn = r.norm();
                        auto q = rn / (baseSupport * (scalar) d);
                        if (q < epsilon || q > scalar(1.0))
                            continue;
                            
                        vec grad = -r / rn * 7.0 / std::numbers::pi / (baseSupport * baseSupport * baseSupport) * (20.0 * q * power(1 - q, 3));

                        normal += grad;
                        // std::cout << "\t" << idx << " " << i << " : " << xi << " : " << yi << " => " << grad.x() << " " << grad.y() << " -> " << W(pos,pb) << " @ " << d << " -> " << q << std::endl;
                    }
                }
            }
            // if(normal.norm() > 0.)
            //normal = normal.normalized();
                // normal = normal / normal.norm();
            // std::cout << i << " @ " << pos.x() << " : " << pos.y() << " -> " << normal.x() << " : " << normal.y() << std::endl;
            normal = normal.normalized();
            boundaryParticles[i] = std::make_pair(pos, normal);

        }



        for (int32_t i = 0; i < boundaryTriangles.size(); ++i) {
            const auto& t = boundaryTriangles[i];
            scalar xmin = std::min(std::min(t.v0.x(), t.v1.x()), t.v2.x()) - baseSupport;
            scalar xmax = std::max(std::max(t.v0.x(), t.v1.x()), t.v2.x()) + baseSupport;
            scalar ymin = std::min(std::min(t.v0.y(), t.v1.y()), t.v2.y()) - baseSupport;
            scalar ymax = std::max(std::max(t.v0.y(), t.v1.y()), t.v2.y()) + baseSupport;

            auto [xmi, ymi] = getCellIdx(xmin, ymin);
            auto [xma, yma] = getCellIdx(xmax, ymax);

            // std::cout << xmi << " : " << ymi << " -> " << xma << " : " << yma << " @ " << cellsX << " : " << cellsY << std::endl;

            for (int32_t x = xmi; x <= xma; ++x) {
                for (int32_t y = ymi; y <= yma; ++y) {
                    //std::cout << x << " : " << y << std::endl;
                    getTriangleCell(x, y).push_back(i);
                }
            }


            //cellTriangleArray
        }

    auto& numPtcls = pm.get<int32_t>("props.numPtcls");
    for(const auto& source: fluidSources){
        if(source.emitter != emitter_t::oneTime){ continue;}
        auto ptcls = source.genParticles();

        auto support = std::sqrt(source.emitterRadius * source.emitterRadius * targetNeighbors);

        for(int32_t i = 0; i < ptcls.size(); ++ i){
        auto [ix, iy] = getCellIdx(ptcls[i].x(), ptcls[i].y());
        bool emit = true;
        for (int32_t xi = -1; xi <= 1; ++xi) {
            for (int32_t yi = -1; yi <= 1; ++yi) {
                const auto& cell = getCell(ix + xi, iy + yi);
                for (auto j : cell) {
                    auto& pj = fluidPosition[j];
                    vec r = pj - ptcls[i];
                    if (r.squaredNorm() <= /*1.5 * 2.0 * 2.0 **/ 
                    0.95 * 0.39960767069134208 * 0.39960767069134208 * support * support) {
                        emit=false;
                    }
                }
            }
        }
        if(!emit)continue;
            fluidPosition[numPtcls]         = ptcls[i];
            fluidVelocity[numPtcls]         = source.emitterVelocity;
            fluidArea[numPtcls]             = source.emitterRadius * source.emitterRadius * std::numbers::pi;
            fluidRestDensity[numPtcls]      = source.emitterDensity;
            fluidSupport[numPtcls]          = std::sqrt(fluidArea[numPtcls] * targetNeighbors / std::numbers::pi);
            fluidPriorPressure[numPtcls]    = 0.;
            fluidVorticity[numPtcls] = 0.;
            fluidAngularVelocity[numPtcls] = 0.;
            fluidUID[numPtcls] = fluidCounter++;

            getCell(ptcls[i].x(), ptcls[i].y()).push_back(numPtcls);

            numPtcls += 1;
        }
    }
    this->boundaryTriangles = boundaryTriangles;
}
