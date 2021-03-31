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


void initRender() {

    glClearColor(1.f, 1.f, 1.f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_POINT_SMOOTH);
    glPointSize(3.f * scale  * pointScale/** 5.f*/);
    glLineWidth(5.f);
    glMatrixMode(GL_PROJECTION);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glClear(GL_COLOR_BUFFER_BIT);
}

scalar velmap(const Particle& p) {
    return p.vel.norm();
}
scalar rhomap(const Particle& p) {
    return p.rho;
}
scalar neighmap(const Particle& p) {
    return (scalar)p.neighbors.size();
}
scalar vorticitymap(const Particle& p) {
    return (scalar)p.angularVelocity;
}
#include <random>
#include "2DMath.h"


enum struct buffer_t {
    velocity, angularVelocity, acceleration, density, neighbors, UID, alpha, area, pressure1, pressure2, pressureBoundary, source, dpdt, rhoStar, predictedVelocity, pressureAcceleration
};
enum struct mode_t {
    magnitude, x, y
};


std::pair<float, float> updateField(float* data) {
    //static std::random_device rd;
    //static std::default_random_engine generator(rd());
    //static std::uniform_real_distribution<double> distribution(0.0, 1.0);
    //std::cout << "Regenerating Data! " << std::endl;
    //for (int32_t i = 0; i < screenWidth * screenHeight; ++i) {
    //    data[i] = distribution(generator);
    //}
    auto bufferstr = ParameterManager::instance().get<std::string>("colorMap.buffer");
    auto modestr = ParameterManager::instance().get<std::string>("colorMap.vectorMode");

    buffer_t buffer;
#define MAPBUFFER(x) if(bufferstr == #x) buffer = buffer_t::x;
    MAPBUFFER(velocity);
    MAPBUFFER(angularVelocity);
    MAPBUFFER(acceleration);
    MAPBUFFER(density);
    MAPBUFFER(neighbors);
    MAPBUFFER(UID);
    MAPBUFFER(alpha);
    MAPBUFFER(area);
    MAPBUFFER(pressure1);
    MAPBUFFER(pressure2);
    MAPBUFFER(pressureBoundary);
    MAPBUFFER(source);
    MAPBUFFER(dpdt);
    MAPBUFFER(rhoStar);
    MAPBUFFER(predictedVelocity);
    MAPBUFFER(pressureAcceleration);
    mode_t mode;
    if (modestr == "magnitude") mode = mode_t::magnitude;
    if (modestr == "x") mode = mode_t::x;
    if (modestr == "y") mode = mode_t::y;





    #pragma omp parallel for
    for (int32_t y = 0; y < screenHeight; ++y) {
        for (int32_t x = 0; x < screenWidth; ++x) {
            scalar xPos = 0.0 + (scalar)x / (scalar)screenWidth * (scalar)domainWidth;
            scalar yPos = 0.0 + (scalar)y / (scalar)screenHeight * (scalar)domainHeight;
            vec pos(xPos, yPos);
            auto [xi, yi] = getCellIdx(xPos, yPos);
            scalar value = 0.0;
            vec xAvg(0.0, 0.0);
            scalar wSum = 0.0;
            //boundaryFunc(pos, [&value](auto bpos, auto d, auto k, auto gk, auto triangle) { value += k; });
            for (int32_t ix = -1; ix <= 1; ++ix) {
                for (int32_t iy = -1; iy <= 1; ++iy) {
                    const auto& ptcls = getCell(std::clamp(xi + ix,0,(int32_t) cellsX-1), std::clamp(yi + iy, 0, (int32_t)cellsY-1));
                    for (auto j : ptcls) {
                        const Particle& p = particles[j];
                        const dfsphState& dp = particlesDFSPH[j];

                        scalar v = 0.0;
                        switch (buffer) {
                        case buffer_t::velocity:
                            v = (mode == mode_t::magnitude
                                ? p.vel.norm()
                                : (mode == mode_t::x ? p.vel.x() : p.vel.y()));
                            break;
                        case buffer_t::angularVelocity: v = p.angularVelocity;  break;
                        case buffer_t::acceleration:
                            v = (mode == mode_t::magnitude
                                ? p.accel.norm()
                                : (mode == mode_t::x ? p.accel.x() : p.accel.y())); break;
                        case buffer_t::density: v = p.rho; break;
                        case buffer_t::neighbors:v = p.neighbors.size(); break;
                        case buffer_t::UID:v = p.uid; break;
                        case buffer_t::alpha:v = dp.alpha; break;
                        case buffer_t::area:v = dp.area; break;
                        case buffer_t::pressure1:v = dp.pressure1; break;
                        case buffer_t::pressure2:v = dp.pressure2; break;
                        case buffer_t::pressureBoundary:v = dp.pressureBoundary; break;
                        case buffer_t::source:v = dp.source; break;
                        case buffer_t::dpdt:v = dp.dpdt; break;
                        case buffer_t::rhoStar:v = dp.rhoStar; break;
                        case buffer_t::predictedVelocity:
                            v = (mode == mode_t::magnitude
                                ? dp.vel.norm()
                                : (mode == mode_t::x ? dp.vel.x() : dp.vel.y())); break;
                        case buffer_t::pressureAcceleration:
                            v = (mode == mode_t::magnitude
                                ? dp.accel.norm()
                                : (mode == mode_t::x ? dp.accel.x() : dp.accel.y())); break;
                        }



                        value += area * v * W(pos, particles[j].pos);
                        //float wi = renderKernel((pos - particles[j].pos).norm() / (support));
                        //xAvg += wi * particles[j].pos;
                        //wSum += wi;
                    }                    
                }
            }
            //if (wSum > 0.0) {
            //    //return 0.0;
            //    xAvg /= wSum;
            //    value = (pos - xAvg).norm();
            //}
            //else
            //    value =  1.0;
            //value = gridPoint(pos, support, 1.0);
            
            data[y * screenWidth + x] = (float) value;
        }
    }

    auto min = *std::min_element(data, data + screenWidth * screenHeight);
    auto max = *std::max_element(data, data + screenWidth * screenHeight);
#pragma omp parallel for
    for (int32_t y = 0; y < screenHeight; ++y) {
        for (int32_t x = 0; x < screenWidth; ++x) {
            auto v = data[y * screenWidth + x];
            v = (v - min) / (max - min);
            v = std::clamp(v, 0.f, 1.f);
            data[y * screenWidth + x] = v;
        }
    }
    return std::make_pair(min, max);

}

void renderField() {
    static GLuint textureID, vao, texID, program, minUni, maxUni, texUnit;
    static bool once = true;
    static float* data = new float[screenWidth * screenHeight];
    if(once){
        std::default_random_engine generator;
        std::uniform_real_distribution<double> distribution(0.0, 1.0);
        for (int32_t i = 0; i < screenWidth * screenHeight; ++i) {
            data[i] = 0.f;
        }
    auto vtxShader = R"(#version 150

// Input vertex data, different for all executions of this shader.
in vec3 vertexPosition_modelspace;

// Output data ; will be interpolated for each fragment.
out vec2 UV;

void main(){
	gl_Position =  vec4(vertexPosition_modelspace,1);
	UV = (vertexPosition_modelspace.xy+vec2(1,1))/2.0;
})";
    auto frgShader = R"(#version 150

in vec2 UV;

out vec3 color;

uniform sampler2D renderedTexture;
uniform sampler1D colorRamp;
uniform float min = 0.0;
uniform float max = 1.0;

void main(){
float v = texture( renderedTexture, UV ).r;
float rel = v;//(v - min) / (max - min);
float clamped = clamp(rel, 0.0, 1.0);
vec3 col = texture(colorRamp, clamped).xyz;
color = col;
	//color =  vec3(v,v,v);
	//color =  vec3(UV.xy,1.0);
})";
    GLuint quad_VertexArrayID;
    glGenVertexArrays(1, &quad_VertexArrayID);
    glBindVertexArray(quad_VertexArrayID);

    static const GLfloat g_quad_vertex_buffer_data[] = {
        -1.0f, -1.0f, 0.0f,
        1.0f, -1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f,
        -1.0f,  1.0f, 0.0f,
        1.0f, -1.0f, 0.0f,
        1.0f,  1.0f, 0.0f,
    };
    // Create one OpenGL texture
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);
    glGenTextures(1, &textureID);

    // "Bind" the newly created texture : all future texture functions will modify this texture
    glBindTexture(GL_TEXTURE_2D, textureID);

    // Give the image to OpenGL
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, screenWidth, screenHeight, 0, GL_RED, GL_FLOAT, data);

    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    GLuint quad_vertexbuffer;
    glGenBuffers(1, &quad_vertexbuffer);
    glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(g_quad_vertex_buffer_data), g_quad_vertex_buffer_data, GL_STATIC_DRAW);

    // Create and compile our GLSL program from the shaders
    program = createProgram(vtxShader, frgShader);

    glUseProgram(program);
    texID = glGetUniformLocation(program, "renderedTexture");
    minUni = glGetUniformLocation(program, "min");
    maxUni = glGetUniformLocation(program, "max");
    //std::cout << "texID: " << texID << std::endl;
    //std::cout << "minUni: " << minUni << std::endl;
    //std::cout << "maxUni: " << maxUni << std::endl;

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, textureID);
    // Set our "myTextureSampler" sampler to use Texture Unit 0
    glUniform1i(texID, 1);
    glEnableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, quad_vertexbuffer);
    glVertexAttribPointer(
        0,                  // attribute 0. No particular reason for 0, but must match the layout in the shader.
        3,                  // size
        GL_FLOAT,           // type
        GL_FALSE,           // normalized?
        0,                  // stride
        (void*)0            // array buffer offset
    );
    once = false;
    int32_t image_width = 1024;
    int32_t image_height = 1024;
    v4* img = new v4[1024];
    v4* color_map = nullptr;
    for (int32_t it = 0; it < 1024; ++it)
        img[it] = v4{ (float)it / (float)1024 * 255.f,(float)it / (float)1024 * 255.f,(float)it / (float)1024 * 255.f, 255.f };
    std::string file_name = "viridis.png";
    if (std::filesystem::exists(file_name)) {
        //std::cout << "Loading " << file_name << std::endl;
        unsigned char* image_data = stbi_load(file_name.c_str(), &image_width, &image_height, NULL, 4);
        delete[] img;
        img = new v4[image_width];
        for (int32_t it = 0; it < image_width; ++it) {
            img[it] = v4{
            (float)image_data[it * 4 + 0],
            (float)image_data[it * 4 + 1],
            (float)image_data[it * 4 + 2], 255.f };
            //std::cout << it << ": [ " << img[it].x << " " << img[it].y << " " << img[it].z <<
            //    " ]" << std::endl;
        }
        //img = QImage(QString::fromStdString(file_name));
        //img.load(QString(file_name.c_str()));
        //std::cout << image_width << " : " << image_height << std::endl;
    }
    //catch (...) {}
    color_map = (v4*)realloc(color_map, sizeof(float4) * (image_width));
    for (int32_t it = 0; it < image_width; ++it) {
        color_map[it] = v4{ (float)(img[it].x) / 256.f, (float)(img[it].y) / 256.f,
                                  (float)(img[it].z) / 256.f, 1.f };
        //std::cout << color_map[it].x() << " : " << color_map[it].y() << " : " << color_map[it].z() << std::endl;
          //if(it == img.width()-1)
          //	color_map[it + 1] = QVector4D{ (float)(col.red()) / 256.f, (float)(col.green()) / 256.f,
          //	(float)(col.blue()) / 256.f, 1.f };
    }
    auto color_map_elements = image_width;

    //std::cout << "color_map_elements " << color_map_elements << std::endl;
    texUnit = create1DTexture(color_map, color_map_elements);
    GLuint samplerLocation = glGetUniformLocation(program, "colorRamp");
    //std::cout << "colorRamp: " << samplerLocation << std::endl;
    glUseProgram(program);
    glUniform1i(samplerLocation, 0);
    glActiveTexture(GL_TEXTURE0 + 0);
    glBindTexture(GL_TEXTURE_1D, texUnit);
    //std::cout << " unit " << texUnit << " ->  " << samplerLocation << " : " << 0 << std::endl;

    glUseProgram(0);
    glBindVertexArray(0);
}
    static  auto& frameR = ParameterManager::instance().get<int32_t>("sim.frame");
    static int32_t frame = -1;
    static auto buff = ParameterManager::instance().get<std::string>("colorMap.buffer");
    static auto vMode = ParameterManager::instance().get<std::string>("colorMap.buffer");
    if (frame != frameR || (buff != ParameterManager::instance().get<std::string>("colorMap.buffer"))||(vMode != ParameterManager::instance().get<std::string>("colorMap.vectorMode"))) {
        buff = ParameterManager::instance().get<std::string>("colorMap.buffer");
        vMode = ParameterManager::instance().get<std::string>("colorMap.vectorMode");
        auto [min,max] = updateField(data);
        //std::cout << buff << " -> " << vMode << std::endl;
        //std::cout << "New min and max: " << min << " x " << max << std::endl;
        ParameterManager::instance().get<scalar>("field.min") = min;
        ParameterManager::instance().get<scalar>("field.max") = max;
        glUseProgram(program);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, textureID);
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, screenWidth, screenHeight, GL_RED, GL_FLOAT, data);
        glUseProgram(0);
        frame = frameR;
    }
    glUseProgram(program);
    //std::cout << textureID << " -> " << texUnit << std::endl;
    glBindVertexArray(vao);
    glUniform1f(minUni, ParameterManager::instance().get<scalar>("field.min"));
    glUniform1f(maxUni, ParameterManager::instance().get<scalar>("field.max"));
    //GLuint samplerLocation = glGetUniformLocation(program, "colorRamp");
    //glUniform1i(samplerLocation, 0);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, textureID);
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_1D, texUnit);
    glDrawArrays(GL_TRIANGLES, 0, 6); // 2*3 indices starting at 0 -> 2 triangles
    glBindVertexArray(0);
    glUseProgram(0);
}

void render() {
    glPointSize(3.f * scale * pointScale/** 5.f*/);
    // used to find the proper scaling for color mapping
    //auto velmap = [](const Particle p) -> scalar { return p.vel.norm(); };
    //auto rhomap = [](const Particle p) -> scalar { return p.rho; };
    //auto neighmap = [](const Particle p) -> scalar {
    //    return (scalar) p.neighbors.size(); 
    //};
    std::function<scalar(const Particle&)> map = neighmap;
    auto bufferstr = ParameterManager::instance().get<std::string>("colorMap.buffer");
    auto modestr = ParameterManager::instance().get<std::string>("colorMap.vectorMode");

    buffer_t buffer;
//#define MAPBUFFER(x) if(bufferstr == #x) buffer = buffer_t::x;
    MAPBUFFER(velocity);
    MAPBUFFER(angularVelocity);
    MAPBUFFER(acceleration);
    MAPBUFFER(density);
    MAPBUFFER(neighbors);
    MAPBUFFER(UID);
    MAPBUFFER(alpha);
    MAPBUFFER(area);
    MAPBUFFER(pressure1);
    MAPBUFFER(pressure2);
    MAPBUFFER(pressureBoundary);
    MAPBUFFER(source);
    MAPBUFFER(dpdt);
    MAPBUFFER(rhoStar);
    MAPBUFFER(predictedVelocity);
    MAPBUFFER(pressureAcceleration);
    mode_t mode;
    if (modestr == "magnitude") mode = mode_t::magnitude;
    if (modestr == "x") mode = mode_t::x;
    if (modestr == "y") mode = mode_t::y;

    //auto min = std::numeric_limits<scalar>::max();
    //auto max = -std::numeric_limits<scalar>::max();
    auto& grid = ParameterManager::instance().get<bool>("render.showGrid");
    auto& marching = ParameterManager::instance().get<bool>("marching.render");
    auto& field = ParameterManager::instance().get<bool>("field.render");
    auto& ptcls = ParameterManager::instance().get<bool>("ptcl.render");
    auto& min = ParameterManager::instance().get<scalar>("colorMap.min");
    auto& max = ParameterManager::instance().get<scalar>("colorMap.max");
    auto& autoScaling = ParameterManager::instance().get<bool>("colorMap.auto");
    auto& limitScope = ParameterManager::instance().get<bool>("colorMap.limit");


    auto m = (domainHeight - 10.0 / domainScale) / 2.0 + domainEpsilon;
    auto t = 12.5 / domainScale;

    vec begin(domainEpsilon + 10/ domainScale, m - t);
    vec end(domainWidth - 15/ domainScale, m + t);
    vec mid = (end + begin) * 0.5;
    mid.x() -= 5.0/ domainScale;

    auto thickness = 1.0/ domainScale;
    auto c = 1.0 / domainScale;

    if (autoScaling) {
        min = std::numeric_limits<scalar>::max();
        max = -std::numeric_limits<scalar>::max();
        int32_t i = 0;
        for (int32_t i = 0; i < particles.size();++i) {
            const Particle& p = particles[i];
            const dfsphState& dp = particlesDFSPH[i];

            if (limitScope&&(p.pos.x() < begin.x() || p.pos.y() < begin.y() || p.pos.x() > end.x() || p.pos.y() > end.y()))
                continue;
            scalar v = 0.0;
            switch (buffer) {
              case buffer_t::velocity:
                v = (mode == mode_t::magnitude
                         ? p.vel.norm()
                         : (mode == mode_t::x ? p.vel.x() : p.vel.y()));
                break;
              case buffer_t::angularVelocity: v = p.angularVelocity;  break;
                case buffer_t::acceleration :
                    v = (mode == mode_t::magnitude
                        ? p.accel.norm()
                        : (mode == mode_t::x ? p.accel.x() : p.accel.y())); break;
                case buffer_t::density: v = p.rho; break;
                case buffer_t::neighbors :v = p.neighbors.size(); break;
                case buffer_t::UID :v = p.uid; break;
                case buffer_t::alpha:v = dp.alpha; break;
                case buffer_t::area:v = dp.area; break;
                case buffer_t::pressure1 :v = dp.pressure1; break;
                case buffer_t::pressure2 :v = dp.pressure2; break;
                case buffer_t::pressureBoundary :v = dp.pressureBoundary; break;
                case buffer_t::source :v = dp.source; break;
                case buffer_t::dpdt :v = dp.dpdt; break;
                case buffer_t::rhoStar :v = dp.rhoStar; break;
                case buffer_t::predictedVelocity :
                    v = (mode == mode_t::magnitude
                        ? dp.vel.norm()
                        : (mode == mode_t::x ? dp.vel.x() : dp.vel.y())); break;
                case buffer_t::pressureAcceleration :
                    v = (mode == mode_t::magnitude
                        ? dp.accel.norm()
                        : (mode == mode_t::x ? dp.accel.x() : dp.accel.y())); break;
            }

            min = std::min(v, min);
            max = std::max(v, max);
        }
    }
    //min = 0.0;
    //max = 25.0;

    // update the window title
    std::stringstream sstream;
    //sstream << "Simulation time: " << simulationTime << "Frame time: " << frameTime << "ms, cmin: " << min << ", cmax: " << max
    //    << ", dt: " << dt << std::endl;
    //glfwSetWindowTitle(window, sstream.str().c_str());

    // reset the openGL state
    //glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glLoadIdentity();
    glOrtho(0, domainWidth, 0, domainHeight, 0, 1);
    //glOrtho(40, 60, 3, 13, 0, 1);

    if (field) {
        renderField();
    }
    if (marching) {
        renderMarching();
    }
     //draw a background grid that shows the cell structure and as a visual aid
    if (grid) {
        glLineWidth(0.5f);
        glColor4f(0.3f, 0.3f, 0.3f, 1);
        glBegin(GL_LINES);
        for (int32_t x = 1; x < (int32_t)cellsX; x += 1) {
            for (int32_t y = 1; y < (int32_t)cellsY; y += 1) {
                if ((y % 5 == 0)) {
                    glLineWidth(1.f);
                    glColor4f(0.5f, 0.5f, 0.5f, 1);
                }
                glVertex2f(0.0, (float)y / (float)cellsY * (float)domainHeight);
                glVertex2f(domainWidth, (float)y / (float)cellsY * (float)domainHeight);
                if ((y % 5 == 0)) {
                    glLineWidth(0.5f);
                    glColor4f(0.3f, 0.3f, 0.3f, 1);
                }
            }
            if ((x % 5 == 0)) {
                glLineWidth(1.f);
                glColor4f(0.5f, 0.5f, 0.5f, 1);
            }
            glVertex2f((float)x / (float)cellsX * (float)domainWidth, 0.0);
            glVertex2f((float)x / (float)cellsX * (float)domainWidth, domainHeight);

            if ((x % 5 == 0)) {
                glLineWidth(0.5f);
                glColor4f(0.3f, 0.3f, 0.3f, 1);
            }
        }
        glEnd();
    }
    // draws the simulation boundary limits as thicker lines
    //glColor4f(1.f, 1.f, 1.f, 1);
    glColor4f(0.f, 0.f, 0.f, 1);
    glLineWidth(5.f);
    constexpr auto renderEps = 0.0;
    glBegin(GL_LINES);
    glVertex2f(domainEpsilon - renderEps, domainEpsilon - renderEps);
    glVertex2f(domainWidth - domainEpsilon + renderEps, domainEpsilon - renderEps);
    glVertex2f(domainWidth - domainEpsilon + renderEps, domainEpsilon - renderEps);
    glVertex2f(domainWidth - domainEpsilon + renderEps, domainHeight - domainEpsilon + renderEps);
    glVertex2f(domainWidth - domainEpsilon + renderEps, domainHeight - domainEpsilon + renderEps);
    glVertex2f(domainEpsilon - renderEps, domainHeight - domainEpsilon + renderEps);
    glVertex2f(domainEpsilon - renderEps, domainHeight - domainEpsilon + renderEps);
    glVertex2f(domainEpsilon - renderEps, domainEpsilon - renderEps);
    glEnd();

    if(ptcls) {

        // draws the particles in the simulation as points with a color value indicating velocity
        glBegin(GL_POINTS);
        for (int32_t i = 0; i < particles.size(); ++i) {
            const Particle& p = particles[i];
            const dfsphState& dp = particlesDFSPH[i];

            scalar v = 0.0;
            switch (buffer) {
            case buffer_t::velocity:
                v = (mode == mode_t::magnitude
                    ? p.vel.norm()
                    : (mode == mode_t::x ? p.vel.x() : p.vel.y()));
                break;
            case buffer_t::angularVelocity: v = p.angularVelocity;  break;
            case buffer_t::acceleration:
                v = (mode == mode_t::magnitude
                    ? p.accel.norm()
                    : (mode == mode_t::x ? p.accel.x() : p.accel.y())); break;
            case buffer_t::density: v = p.rho; break;
            case buffer_t::neighbors:v = p.neighbors.size(); break;
            case buffer_t::UID:v = p.uid; break;
            case buffer_t::alpha:v = dp.alpha; break;
            case buffer_t::area:v = dp.area; break;
            case buffer_t::pressure1:v = dp.pressure1; break;
            case buffer_t::pressure2:v = dp.pressure2; break;
            case buffer_t::pressureBoundary:v = dp.pressureBoundary; break;
            case buffer_t::source:v = dp.source; break;
            case buffer_t::dpdt:v = dp.dpdt; break;
            case buffer_t::rhoStar:v = dp.rhoStar; break;
            case buffer_t::predictedVelocity:
                v = (mode == mode_t::magnitude
                    ? dp.vel.norm()
                    : (mode == mode_t::x ? dp.vel.x() : dp.vel.y())); break;
            case buffer_t::pressureAcceleration:
                v = (mode == mode_t::magnitude
                    ? dp.accel.norm()
                    : (mode == mode_t::x ? dp.accel.x() : dp.accel.y())); break;
            }


            auto col = tinycolormap::GetColor((v - min) / (max - min), tinycolormap::ColormapType::Viridis);
            glColor4f((GLfloat)col.r(), (GLfloat)col.g(), (GLfloat)col.b(), 1);
            glVertex2f((GLfloat)p.pos(0), (GLfloat)p.pos(1));
        }
        glEnd();

    }
    glColor4f(1.f, 0.0f, 0.0f, 1);
    glLineWidth(2.5f);
    glBegin(GL_LINES);
    for (auto& tri : triangles) {
        glVertex2f(tri.v0.x(), tri.v0.y());
        glVertex2f(tri.v1.x(), tri.v1.y());
        glVertex2f(tri.v1.x(), tri.v1.y());
        glVertex2f(tri.v2.x(), tri.v2.y());
        glVertex2f(tri.v2.x(), tri.v2.y());
        glVertex2f(tri.v0.x(), tri.v0.y());
    }
    glEnd();


    renderRays();

    //glColor4f(0.f, 0.f, 0.f, 1);
    //glLineWidth(3.5f);
    //glBegin(GL_LINES);
    //glVertex2f(5, 5);
    //glVertex2f(95, 5);
    //glVertex2f(95, 5);
    //glVertex2f(95, 95);
    //glEnd();

    //glColor4f(0.f, 0.f, 0.f, 1);
    //glLineWidth(3.5f);
    //glBegin(GL_LINES);
    //glVertex2f(5, 25);
    //glVertex2f(50, 25);
    //glVertex2f(50, 25);
    //if (simulationCase == cornerAngle::acute) {
    //  glVertex2f(25, 45);
    //  glVertex2f(25, 45);
    //}
    //if (simulationCase == cornerAngle::ortho) {
    //  glVertex2f(50, 45);
    //  glVertex2f(50, 45);
    //}  
    //if (simulationCase == cornerAngle::obtuse) {
    //  glVertex2f(75, 45);
    //  glVertex2f(75, 45);
    //}

    //if (simulationCase == cornerAngle::nobtuse) {
    //  glVertex2f(75, 5);

    //  glVertex2f(75, 5);
    //  glVertex2f(95, 5);

    //  glVertex2f(95, 5);
    //  glVertex2f(95, 45);

    //  glVertex2f(95, 45);
    //}
    //if (simulationCase == cornerAngle::northo) {
    //  glVertex2f(50, 5);

    //  glVertex2f(50, 5);
    //  glVertex2f(95, 5);

    //  glVertex2f(95, 5);
    //  glVertex2f(95, 45);

    //  glVertex2f(95, 45);
    //}
    //if (simulationCase == cornerAngle::nacute) {
    //  glVertex2f(25, 5);

    //  glVertex2f(25, 5);
    //  glVertex2f(95, 5);

    //  glVertex2f(95, 5);
    //  glVertex2f(95, 45);

    //  glVertex2f(95, 45);
    //}

    glVertex2f(5, 45);
    glVertex2f(5, 45);
    glVertex2f(5, 25);
    glEnd();

    glColor4f(0.8f, 0.f, 0.f, 1);
    glLineWidth(1.5f);
    glBegin(GL_LINES);
    vec l(5, 25);
    vec u(20, 40);

    //glVertex2f(l.x(), l.y());
    //glVertex2f(u.x(), l.y());

    //glVertex2f(u.x(), l.y());
    //glVertex2f(u.x(), u.y());

    //glVertex2f(u.x(), u.y());
    //glVertex2f(l.x(), u.y());

    //glVertex2f(l.x(), u.y());
    //glVertex2f(l.x(), l.y());
    glEnd();

}