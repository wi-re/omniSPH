#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"
#include <glad/glad.h> 
#include "glui.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include <tools/stb_image_write.h>
#define STB_IMAGE_READ_IMPLEMENTATION
#include <tools/stb_image.h>
#include <simulation/SPH.h>
#include <tools/ParameterManager.h>


struct v4 { float x, y, z, w; };
GLuint create1DTexture(v4* colorMap, int32_t elements);
GLuint createProgram(std::string vertexSource, std::string fragmentSource);



GLuint create1DTexture(v4* colorMap, int32_t elements) {
    GLuint textureId_;

    // generate the specified number of texture objects
    glGenTextures(1, &textureId_);
     assert(glGetError() == GL_NO_ERROR);

    // bind texture
    glBindTexture(GL_TEXTURE_1D, textureId_);
     assert(glGetError() == GL_NO_ERROR);

    // tells OpenGL how the data that is going to be uploaded is aligned
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
     assert(glGetError() == GL_NO_ERROR);

    glTexImage1D(
        GL_TEXTURE_1D, // Specifies the target texture. Must be GL_TEXTURE_1D or GL_PROXY_TEXTURE_1D.
        0, // Specifies the level-of-detail number. Level 0 is the base image level. Level n is the
           // nth mipmap reduction image.
        GL_RGBA32F, elements,
        0, // border: This value must be 0.
        GL_RGBA, GL_FLOAT, colorMap);
     assert(glGetError() == GL_NO_ERROR);

    // texture sampling/filtering operation.
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
     assert(glGetError() == GL_NO_ERROR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
     assert(glGetError() == GL_NO_ERROR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
     assert(glGetError() == GL_NO_ERROR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
     assert(glGetError() == GL_NO_ERROR);

    glBindTexture(GL_TEXTURE_1D, 0);
    // assert(glGetError() == GL_NO_ERROR);

    return textureId_;
}
#define STB_IMAGE_IMPLEMENTATION
#include <tools/stb_image.h>

GLuint createProgram(std::string vertexSource, std::string fragmentSource) {
    // Create an empty vertex shader handle
    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);

    // Send the vertex shader source code to GL
    // Note that std::string's .c_str is NULL character terminated.
    const GLchar* source = (const GLchar*)vertexSource.c_str();
    glShaderSource(vertexShader, 1, &source, 0);

    // Compile the vertex shader
    glCompileShader(vertexShader);

    GLint isCompiled = 0;
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &isCompiled);
    if (isCompiled == GL_FALSE) {
        GLint maxLength = 0;
        glGetShaderiv(vertexShader, GL_INFO_LOG_LENGTH, &maxLength);

        // The maxLength includes the NULL character
        std::vector<GLchar> infoLog(maxLength);
        glGetShaderInfoLog(vertexShader, maxLength, &maxLength, &infoLog[0]);

        // We don't need the shader anymore.
        glDeleteShader(vertexShader);

        // Use the infoLog as you see fit.
        std::cerr << infoLog.data() << std::endl;

        // In this simple program, we'll just leave
        return -1;
    }

    // Create an empty fragment shader handle
    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

    // Send the fragment shader source code to GL
    // Note that std::string's .c_str is NULL character terminated.
    source = (const GLchar*)fragmentSource.c_str();
    glShaderSource(fragmentShader, 1, &source, 0);

    // Compile the fragment shader
    glCompileShader(fragmentShader);

    glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &isCompiled);
    if (isCompiled == GL_FALSE) {
        GLint maxLength = 0;
        glGetShaderiv(fragmentShader, GL_INFO_LOG_LENGTH, &maxLength);

        // The maxLength includes the NULL character
        std::vector<GLchar> infoLog(maxLength);
        glGetShaderInfoLog(fragmentShader, maxLength, &maxLength, &infoLog[0]);

        // We don't need the shader anymore.
        glDeleteShader(fragmentShader);
        // Either of them. Don't leak shaders.
        glDeleteShader(vertexShader);
        std::cerr << infoLog.data() << std::endl;

        // Use the infoLog as you see fit.

        // In this simple program, we'll just leave
        return -1;
    }

    // Vertex and fragment shaders are successfully compiled.
    // Now time to link them together into a program.
    // Get a program object.
    GLuint program = glCreateProgram();

    // Attach our shaders to our program
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);

    // Link our program
    glLinkProgram(program);

    // Note the different functions here: glGetProgram* instead of glGetShader*.
    GLint isLinked = 0;
    glGetProgramiv(program, GL_LINK_STATUS, (int*)&isLinked);
    if (isLinked == GL_FALSE) {
        GLint maxLength = 0;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &maxLength);

        // The maxLength includes the NULL character
        std::vector<GLchar> infoLog(maxLength);
        glGetProgramInfoLog(program, maxLength, &maxLength, &infoLog[0]);
        std::cerr << infoLog.data() << std::endl;

        // We don't need the program anymore.
        glDeleteProgram(program);
        // Don't leak shaders either.
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);

        // Use the infoLog as you see fit.

        // In this simple program, we'll just leave
        return -1;
    }

    // Always detach shaders after a successful link.
    glDetachShader(program, vertexShader);
    glDetachShader(program, fragmentShader);
    return program;
}


//scalar velmap(const Particle& p) {
//    return p.vel.norm();
//}
//scalar rhomap(const Particle& p) {
//    return p.rho;
//}
//scalar neighmap(const Particle& p) {
//    return (scalar)p.neighbors.size();
//}
//scalar vorticitymap(const Particle& p) {
//    return (scalar)p.angularVelocity;
//}
#include <random>

enum struct buffer_t {
    velocity, angularVelocity, acceleration, density, neighbors, UID, alpha, area, pressure1, pressure2, pressureBoundary, source, dpdt, rhoStar, predictedVelocity, pressureAcceleration, vorticity
};
enum struct vecmode_t {
    magnitude, x, y
};



void render(int32_t screenWidth, int32_t screenHeight) {
    glPointSize(1.f);
    // used to find the proper scaling for color mapping
    //auto velmap = [](const Particle p) -> scalar { return p.vel.norm(); };
    //auto rhomap = [](const Particle p) -> scalar { return p.rho; };
    //auto neighmap = [](const Particle p) -> scalar {
    //    return (scalar) p.neighbors.size(); 
    //};

    // update the window title
    std::stringstream sstream;
    //sstream << "Simulation time: " << simulationTime << "Frame time: " << frameTime << "ms, cmin: " << min << ", cmax: " << max
    //    << ", dt: " << dt << std::endl;
    //glfwSetWindowTitle(window, sstream.str().c_str());

    // reset the openGL state
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glOrtho(
        simulationState.pm.get<vec>("domain.min").x(), simulationState.pm.get<vec>("domain.max").x(), 
        simulationState.pm.get<vec>("domain.min").y(), simulationState.pm.get<vec>("domain.max").y(), 
    0, 1);

    #define MAP_PROPERTY( x)\
        if(arg == #x) prop = SPHSimulation::property_t::x;
    auto arg = simulationState.pm.get<std::string>("colorMap.buffer");
    SPHSimulation::property_t prop = SPHSimulation::property_t::density;

    MAP_PROPERTY(position);
    MAP_PROPERTY(velocity);
    MAP_PROPERTY(accel);
    MAP_PROPERTY(predVelocity);
    MAP_PROPERTY(predAccel);
    MAP_PROPERTY(predPosition);

    MAP_PROPERTY(density);
    MAP_PROPERTY(vorticity);
    MAP_PROPERTY(angularVelocity);
    MAP_PROPERTY(area);
    MAP_PROPERTY(restDensity);
    MAP_PROPERTY(support);

    MAP_PROPERTY(alpha);
    MAP_PROPERTY(actualArea);
    MAP_PROPERTY(pressure1);
    MAP_PROPERTY(pressure2);
    MAP_PROPERTY(boundaryPressure);
    MAP_PROPERTY(sourceTerm);
    MAP_PROPERTY(dpdt);
    MAP_PROPERTY(rhoStar);
    MAP_PROPERTY(priorPressure);

    MAP_PROPERTY(UID);
    MAP_PROPERTY(neighbors);


  static std::string colorMap = "";
static v4 *color_map = nullptr;
static int32_t image_width = 1024;
static int32_t image_height = 1024;
  // std::cout << colorMap << " : " <<
  // ParameterManager::instance().get<std::string>("colorMap.colorMap") <<
  // std::endl;
  if (simulationState.pm.get<std::string>("colorMap.cmap") != colorMap) {
    colorMap = simulationState.pm.get<std::string>("colorMap.cmap");
    
    v4 *img = new v4[1024];
    for (int32_t it = 0; it < 1024; ++it)
      img[it] = v4{(float)it / (float)1024 * 255.f, (float)it / (float)1024 * 255.f, (float)it / (float)1024 * 255.f, 255.f};
    // std::string file_name = std::string("cfg/") + colorMap + ".png";
    try {
      std::string file_name = simulationState.resolveFile(std::string("cfg/") + simulationState.pm.get<std::string>("colorMap.cmap") + ".png").string();

      if (fs::exists(file_name)) {
        // std::cout << "Loading " << file_name << std::endl;
        int32_t comp;
        unsigned char *image_data = stbi_load(file_name.c_str(), &image_width, &image_height, &comp, 4);
        // std::cout << image_width << " " << image_height << " " << comp << std::endl;
        delete[] img;
        img = new v4[image_width];
        for (int32_t it = 0; it < image_width; ++it) {
          img[it] = v4{(float)image_data[it * 4 + 0 + image_width * (image_height / 2)], (float)image_data[it * 4 + 1 + image_width * (image_height / 2)],
                       (float)image_data[it * 4 + 2 + image_width * (image_height / 2)], 255.f};
          // std::cout << it << ": [ " << img[it].x << " " << img[it].y << " "
          // << img[it].z <<
          //     " ]" << std::endl;
        }
        // img = QImage(QString::fromStdString(file_name));
        // img.load(QString(file_name.c_str()));
        // std::cout << image_width << " : " << image_height << std::endl;
      }else{
        std::cout << "Color Map " << (std::string("cfg/") + simulationState.pm.get<std::string>("colorMap.cmap") + ".png") << "does not exist" << std::endl;
      }
    } catch (...) {
        std::cout << "Color Map " << (std::string("cfg/") + simulationState.pm.get<std::string>("colorMap.cmap") + ".png") << "does not exist" << std::endl;
    }
    color_map = (v4 *)realloc(color_map, sizeof(v4) * (image_width));
    for (int32_t it = 0; it < image_width; ++it) {
      color_map[it] = v4{(float)(img[it].x) / 256.f, (float)(img[it].y) / 256.f, (float)(img[it].z) / 256.f, 1.f};
      // std::cout << color_map[it].x() << " : " << color_map[it].y() << " : "
      // << color_map[it].z() << std::endl; if(it == img.width()-1)
      //	color_map[it + 1] = QVector4D{ (float)(col.red()) / 256.f,
      //(float)(col.green()) / 256.f, 	(float)(col.blue()) / 256.f, 1.f };
    }
    auto color_map_elements = image_width;
    delete[] img;

    // std::cout << "color_map_elements " << color_map_elements << std::endl;
    // texUnit = create1DTexture(color_map, color_map_elements);
  }

    {
        // draws the particles in the simulation as points with a color value indicating velocity
        auto& particlePositions = simulationState.getPositions();
        glColor4f((GLfloat)1.f, (GLfloat)0.f, (GLfloat)0.f, 1);
        auto& numPtcls = simulationState.pm.get<int32_t>("props.numPtcls");
        auto [min, max, vs] = simulationState.colorMap(prop, simulationState.pm.get<bool>("colorMap.auto"), simulationState.pm.get<scalar>("colorMap.min"), simulationState.pm.get<scalar>("colorMap.max"));
        simulationState.pm.get<scalar>("colorMap.min") = min;
        simulationState.pm.get<scalar>("colorMap.max") = max;
        auto flipped = simulationState.pm.get<bool>("colorMap.flip");

            // auto area = simulationState.fluidArea[i];
            auto radius =simulationState.pm.get<scalar>("props.baseRadius");
            auto ratio = screenWidth / (simulationState.pm.get<vec>("domain.max").x() - simulationState.pm.get<vec>("domain.min").x());
            glPointSize(ratio * radius);
            // std::cout << radius << " -> " << ratio << " => " << ratio * radius << std::endl;

        glBegin(GL_POINTS);
        for(int32_t i = 0; i< numPtcls; ++ i){
            auto s = vs[i];
            s = std::clamp(s, 0., 1.);
            if(flipped) s = 1. - s;
            auto ss = s * (scalar) (image_width - 1);
            int32_t il = std::clamp(int32_t(::floor(ss)),0, image_width - 1);
            int32_t ir = std::clamp(int32_t(::ceil(ss)),0, image_width - 1);
            auto frac = std::clamp(ss - ::floor(ss),0.,1.);
            auto cl = color_map[il];
            auto cr = color_map[ir];
            
            // auto col = tinycolormap::GetColor(vs[i], tinycolormap::ColormapType::Viridis);
            glColor4f(  (GLfloat)(cl.x * (1. - frac) + cr.x * frac), 
                        (GLfloat)(cl.y * (1. - frac) + cr.y * frac), 
                        (GLfloat)(cl.z * (1. - frac) + cr.z * frac), 1);
            glVertex2f((GLfloat)particlePositions[i](0), (GLfloat)particlePositions[i](1));
        }
        glEnd();

    }
    {
        glBegin(GL_POINTS);
        for (int32_t i = 0; i < simulationState.boundaryParticles.size(); ++i) {
            auto [p,n] = simulationState.boundaryParticles[i];
            glColor4f(n.x() * 0.5 + 0.5,n.y() * 0.5 + 0.5,0.5f, 1);
            glVertex2f((GLfloat)p.x(), (GLfloat)p.y());
        }
        glEnd();

    }
    glColor4f(1.f, 0.0f, 0.0f, 1);
    glLineWidth(.25f);
    glBegin(GL_LINES);
    for (auto& tri : simulationState.boundaryTriangles) {
        glVertex2f(tri.v0.x(), tri.v0.y());
        glVertex2f(tri.v1.x(), tri.v1.y());
        glVertex2f(tri.v1.x(), tri.v1.y());
        glVertex2f(tri.v2.x(), tri.v2.y());
        glVertex2f(tri.v2.x(), tri.v2.y());
        glVertex2f(tri.v0.x(), tri.v0.y());
    }
    glEnd();
}

GLuint shader_programme;
GLuint vao = 0;
void GUI::renderFunctions() {
    glPushAttrib(GL_ALL_ATTRIB_BITS);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glEnable(GL_MULTISAMPLE);
    glClearColor(0.2f, 0.2f, 0.2f, 1.f);
    // glClearColor(1.0f, 1.0f, 1.0f, 1.f);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHT0);
    glDisable(GL_LIGHTING);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glDisable(GL_COLOR_MATERIAL);
    glFlush();
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
   
        auto b_init = simulationState.pm.get<std::string>("colorMap.buffer");
        auto mi_init = simulationState.pm.get<scalar>("colorMap.min");
        auto ma_init = simulationState.pm.get<scalar>("colorMap.max");
        render(screenWidth, screenHeight);
    static int32_t i = 0;
    if (simulationRunning && ((i++) % 5 == 0)){
        simulationState.timestep();
    }

    // ACTUAL RENDERING CALL

    auto time = simulationState.pm.get<scalar>("sim.time");
    if (time > 30.0)
        shouldStop = true;


    glBindVertexArray(0);
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glPopAttrib();
}
void GUI::renderLoop() {
    show_demo_window = true;
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    while (!shouldStop && !glfwWindowShouldClose(window)) {
        glfwPollEvents();
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
        glClear(GL_COLOR_BUFFER_BIT);

        renderFunctions();
        auto dt = 1.0 / simulationState.pm.get<scalar>("video.fps");
        static auto lastTime = simulationState.pm.get<double>("sim.time") - dt;
        auto currTime = simulationState.pm.get<double>("sim.time");
        if (currTime > lastTime + dt && m_ffmpegPipe != nullptr) {
            lastTime = currTime;
            static int32_t* buffer = new int32_t[display_w * display_h];
            glReadPixels(0, 0, display_w, display_h, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
            fwrite(buffer, sizeof(int) * display_w * display_h, 1, m_ffmpegPipe);

        }

        ImGui::NewFrame();
        //ImGui::ShowDemoWindow(&show_demo_window);
        uiFunctions();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable)
        {
            GLFWwindow* backup_current_context = glfwGetCurrentContext();
            ImGui::UpdatePlatformWindows();
            ImGui::RenderPlatformWindowsDefault();
            glfwMakeContextCurrent(backup_current_context);
        }
        glfwSwapBuffers(window);

        if(simulationState.pm.get<bool>("export.active") &&
            simulationState.pm.get<double>("export.limit") > 0. &&
            currTime > simulationState.pm.get<double>("export.limit"))
            shouldStop = true;

    }
#ifdef WIN32
    if (m_ffmpegPipe != nullptr)
        _pclose(m_ffmpegPipe);
#else
    if (m_ffmpegPipe != nullptr)
        pclose(m_ffmpegPipe);
#endif
    quit();
}
