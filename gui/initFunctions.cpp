#include "imgui/imgui_impl_glfw.h"
#include "imgui/imgui_impl_opengl3.h"
#include <glad/glad.h> 
#include "glui.h"
#pragma warning( push )
#pragma warning(disable: 4251 4275 )
#include <yaml-cpp/yaml.h>
#pragma warning(pop)

static void glfw_error_callback(int error, const char* description)
{
    fprintf(stderr, "Glfw Error %d: %s\n", error, description);
}
void GUI::initGL(int argc, char* argv[]) {
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
        throw std::runtime_error("Could not setup glfw context");
    const char* glsl_version = "#version 150";
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_COMPAT_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    auto domainMin = simulationState.pm.get<vec>("domain.min");
    auto domainMax = simulationState.pm.get<vec>("domain.max");
    auto domainWidth = domainMax.x() - domainMin.x();
    auto domainHeight = domainMax.y() - domainMin.y();
    if(domainWidth > domainHeight){
        screenWidth = 2048;
        screenHeight =  2048. * domainHeight / domainWidth;
    }
    else{
        screenWidth = 2048 * domainWidth / domainHeight;
        screenHeight =  2048.;
    }
    window = glfwCreateWindow(screenWidth, screenHeight, "OmniFlow", NULL, NULL);
    // window = glfwCreateWindow(1024, 1024, "OmniFlow", NULL, NULL);
    if (window == NULL)
        throw std::runtime_error("Could not setup window");
    glfwSetWindowPos(window, 5, 35);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    bool err = gladLoadGL() == 0;
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);
    if (err) {
        fprintf(stderr, "Failed to initialize OpenGL loader!\n");
        throw std::runtime_error("Failed to initialize OpenGL loader");
    }
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    io.ConfigFlags |= ImGuiConfigFlags_ViewportsEnable;
    ImGui::StyleColorsLight();
    ImGuiStyle& style = ImGui::GetStyle();
    if (io.ConfigFlags & ImGuiConfigFlags_ViewportsEnable) {
        style.WindowRounding = 0.0f;
        style.Colors[ImGuiCol_WindowBg].w = 1.0f;
    }
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init(glsl_version);
    ImGui::PushStyleVar(ImGuiStyleVar_ChildRounding, 5.0f);
    io.Fonts->AddFontFromFileTTF("Cascadia.ttf", 16.0f);
    io.Fonts->AddFontFromFileTTF("NotoMono-Regular.ttf", 16.0f);
    clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
    //LoggingWindow::instance().AddLog("Starting log...");

    glfwSetKeyCallback(window, [](GLFWwindow* window, int key, int scancode, int action, int mods) {GUI::instance().keyCallback(window, key, scancode, action, mods); });
    glfwSetCursorPosCallback(window, [](GLFWwindow* window, double xpos, double ypos) {GUI::instance().cursorPositionCallback(window, xpos, ypos); });
    glfwSetMouseButtonCallback(window, [](GLFWwindow* window, int button, int action, int mods) {GUI::instance().mouseButtonCallback(window, button, action, mods); });
    glfwSetScrollCallback(window, [](GLFWwindow* window, double xpos, double ypos) {GUI::instance().scrollCallback(window, xpos, ypos); });
    glfwSwapInterval(0);
    

    glClearColor(1.f, 1.f, 1.f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT);
    glEnable(GL_POINT_SMOOTH);
    glLineWidth(5.f);
    glMatrixMode(GL_PROJECTION);
    glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
    glClear(GL_COLOR_BUFFER_BIT);

    
    if (simulationState.pm.get<bool>("video.active")) {
   namespace fs = std::filesystem;
        
   fs::path basePath = simulationState.pm.get<std::string>("internal.working_directory");
std::cout << basePath << std::endl;

       const std::chrono::time_point<std::chrono::system_clock> now =
           std::chrono::system_clock::now();
       const std::time_t t_c = std::chrono::system_clock::to_time_t(now);
       auto time = std::time(nullptr);
       std::stringstream ss;
       char buf[256];
       ss << std::put_time(std::localtime(&t_c), "%F_%T"); // ISO 8601 without timezone information.
       auto s = ss.str();
       std::replace(s.begin(), s.end(), ':', '-');
       fs::path actualPath = basePath / (s + ".mp4");
    //    actualPath.replace_extension(".mp4");


        std::stringstream sstream;
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        sstream << "ffmpeg -r " << 60.0 // ParameterManager::instance().get<scalar>("sim.maxDt")
            //<< (!arguments::cmd::instance().vm.count("verbose") ? " -hide_banner -nostats -loglevel 0" : "")
            //<< " -hide_banner -nostats -loglevel 0"
            << " -f rawvideo -pix_fmt rgba -s "<< screenWidth << "x" << screenHeight <<  " -i - "
            "-threads 0 -preset ultrafast -y -vf vflip -c:v libx264 "
            "-pix_fmt yuv420p -b:v 100M " 
            << actualPath.string();
            std::cout << sstream.str() << std::endl;
#ifdef WIN32
        m_ffmpegPipe = _popen(sstream.str().c_str(), "wb");
#else
        m_ffmpegPipe = popen(sstream.str().c_str(), "w");
#endif
    }
    else
        m_ffmpegPipe = nullptr;

}

void GUI::initSimulation() {
    // simulationState.initializeParameters();
    // simulationState.initializeSPH();
}

std::filesystem::path resolveFileLocal(std::string fileName, std::vector<std::string> search_paths) {
    namespace fs = std::filesystem;
    fs::path expanded = fs::path(fileName);

    fs::path base_path = "";
    if (fs::exists(fs::path(fileName)))
        return fs::path(fileName);
    for (const auto& path : search_paths) {
        auto p = (fs::path(path));
        if (fs::exists(p / fileName))
            return p.string() + std::string("/") + fileName;
    }

    if (fs::exists(fileName)) return fs::path(fileName);
    if (fs::exists(expanded))
        return expanded;

    for (const auto& pathi : search_paths) {
        auto path = (fs::path(pathi));
    }

    std::stringstream sstream;
    sstream << "File '" << fileName << "' could not be found in any provided search path" << std::endl;
    throw std::runtime_error(sstream.str().c_str());
}

void GUI::initParameters(int argc, char* argv[]) {
    // auto& pm = ParameterManager::instance();
    std::string stylePath;
    std::string fileName = "cfg/style.css";
    namespace fs = std::filesystem;
    fs::path working_dir = fs::absolute(fs::path(argv[0])).remove_filename();
    fs::path binary_dir = fs::current_path();
    fs::path expanded = (fs::path(fileName));
    //stylePath = resolveFileLocal(fileName, { working_dir.string(),binary_dir.string()}).string();
    //pm.newParameter("stylePath", stylePath, { .hidden = true });
    //pm.init();
    auto binary_directory = fs::absolute(fs::path(argv[0])).remove_filename();
    auto working_directory = fs::current_path();

    // Initialize the simulation based on command line parameters.
}


void GUI::initGVDB() {
} 