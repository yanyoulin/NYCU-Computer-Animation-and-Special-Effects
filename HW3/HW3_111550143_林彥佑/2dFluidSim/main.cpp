
#define GLAD_GL_IMPLEMENTATION
#include "glad/gl.h"
#undef GLAD_GL_IMPLEMENTATION

#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>
#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

#include "constants.h"
#include "fluid.h"

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <Eigen/Dense>
#include <vector>
#include <iostream>
#include <algorithm>

int windowW = 1600, windowH = 900;
int particle_num = 8000;
int iterations = 3;                     // number of iterations for solvers

float   particle_radius = 4.0f,
        obstacle_radius = 50.0f,
        flip_ratio = 0.9f,              // in the range of [0, 1]
        stiffness_coefficient = 1.0f,
        cell_dim = 10.0f,               // might want to make this divide windowW and windowH
        relaxation_cell_dim = 10.0f,    // might want to make this divide windowW and windowH
        dt = 0.075f,
        grav = 10.0f,
        over_relaxation = 1.9f;

bool density_correction = true;

GLuint point_program, obstacle_program;
Fluid* fluid;

int render_option = 2;
int lastMouseX = windowW, lastMouseY = 0, dx = 0, dy = 0;
bool lmbDown = false;

static GLuint createShaderProgram(const char* v_shader, const char* f_shader)
{
    GLuint vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertex_shader, 1, &v_shader, NULL);
    glCompileShader(vertex_shader);

    GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragment_shader, 1, &f_shader, NULL);
    glCompileShader(fragment_shader);

    GLuint program = glCreateProgram();
    glAttachShader(program, vertex_shader);
    glAttachShader(program, fragment_shader);
    glLinkProgram(program);

    GLint success;
    GLchar infoLog[512];

    // Check vertex shader compilation
    glGetShaderiv(vertex_shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(vertex_shader, 512, NULL, infoLog);
        std::cerr << "Vertex Shader Compilation Error:\n" << infoLog << std::endl;
    }

    // Check fragment shader compilation
    glGetShaderiv(fragment_shader, GL_COMPILE_STATUS, &success);
    if (!success)
    {
        glGetShaderInfoLog(fragment_shader, 512, NULL, infoLog);
        std::cerr << "Fragment Shader Compilation Error:\n" << infoLog << std::endl;
    }

    // Check program linking
    glGetProgramiv(program, GL_LINK_STATUS, &success);
    if (!success)
    {
        glGetProgramInfoLog(program, 512, NULL, infoLog);
        std::cerr << "Shader Program Linking Error:\n" << infoLog << std::endl;
    }


    return program;
}

static void error_callback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}

static void framebuffer_size_callback(GLFWwindow* window, int width, int height) 
{
    glViewport(0, 0, width, height);
    windowW = width;
    windowH = height;
}

static void mouse_button_callback(GLFWwindow* window, int button, int action, int mods) 
{
    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureMouse)
        return;

    if (button == GLFW_MOUSE_BUTTON_LEFT) 
    {
        if (action == GLFW_PRESS) 
        {
            lmbDown = true;
            double xpos, ypos;
            glfwGetCursorPos(window, &xpos, &ypos);
            lastMouseX = static_cast<int>(xpos);
            lastMouseY = static_cast<int>(ypos);
        }
        else if (action == GLFW_RELEASE) 
        {
            lmbDown = false;
        }
    }
}

static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos) 
{
    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureMouse) 
        return;

    if (lmbDown) 
    {
        int x = static_cast<int>(xpos);
        int y = static_cast<int>(ypos);
        dx = x - lastMouseX;
        dy = y - lastMouseY;
        lastMouseX = x;
        lastMouseY = y;
    }
}

static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods) 
{
    ImGuiIO& io = ImGui::GetIO();
    if (io.WantCaptureKeyboard)
        return;

    if (action == GLFW_PRESS) 
    {
        switch (key) 
        {
        case GLFW_KEY_ESCAPE:
            glfwSetWindowShouldClose(window, GLFW_TRUE);
            break;
        case GLFW_KEY_1:
            render_option = 2;
            break;
        case GLFW_KEY_2:
            render_option = 0;
            break;
        case GLFW_KEY_3:
            render_option = 1;
            break;
        }
    }
}

int main(int argc, char** argv) {
    // Initialize GLFW
    if (!glfwInit())
        return -1;

    // Configure GLFW
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_RESIZABLE, GL_FALSE);

    // Create window
    GLFWwindow* window = glfwCreateWindow(windowW, windowH, "Fluid Simulation | ID: 111550143", NULL, NULL);
    if (!window) 
    {
        glfwTerminate();
        return -1;
    }

    // Make OpenGL context current
    glfwMakeContextCurrent(window);
    glfwSwapInterval(0); // Disable vsync

    // Set up callbacks
    glfwSetErrorCallback(error_callback);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
    glfwSetMouseButtonCallback(window, mouse_button_callback);
    glfwSetCursorPosCallback(window, cursor_position_callback);
    glfwSetKeyCallback(window, key_callback);
    
    // Load OpenGL functions
    if (!gladLoadGL(glfwGetProcAddress))
    {
        glfwTerminate();
        return -1;
    }

    // Set up shaders
    point_program = createShaderProgram(point_vs_src.c_str(), point_fs_src.c_str());
    obstacle_program = createShaderProgram(mono_color_vs_src.c_str(), mono_color_fs_src.c_str());

    // Enable OpenGL features
    glEnable(GL_PROGRAM_POINT_SIZE);

    // Initialize ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL3_Init("#version 330");


    // Create fluid simulation
    fluid = new Fluid(particle_num, particle_radius, obstacle_radius, flip_ratio, 
        cell_dim, relaxation_cell_dim, iterations, windowW, windowH, dt, grav, stiffness_coefficient, density_correction, over_relaxation);
    fluid->setupRendering(point_program, obstacle_program);
    
    double frameDuration = 0;
    // Main loop
    while (!glfwWindowShouldClose(window)) 
    {
        double frameStart = glfwGetTime();

        glfwPollEvents();

        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        glClear(GL_COLOR_BUFFER_BIT);

        ImGui::SetNextWindowSize(ImVec2(400, 300), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Once);
        ImGui::Begin("Control Panel");
        ImGui::Text("Current FPS: %d", static_cast<int>(1 / frameDuration));
        ImGui::Text("Press 1, 2, 3 to switch render option.");

        switch (render_option)
        {
        case 2:
            ImGui::Text("Current Option: Particles");
            break;
        case 0:
            ImGui::Text("Current Option: Cells");
            break;
        case 1:
            ImGui::Text("Current Option: Particles + Cells");
            break;
        default:
            ImGui::Text("oops ig somethings wrong");
            break;
        }
        
        ImGui::SliderFloat("FLIP Ratio", &flip_ratio, 0.0f, 1.0f);
        ImGui::InputFloat("Stiffness Coefficient", &stiffness_coefficient, 0.1f);
        ImGui::InputFloat("dt", &dt, 0.01f);
        ImGui::InputFloat("Gravity", &grav, 5.0f);
        ImGui::InputFloat("Over Relaxation", &over_relaxation, 0.1f);
        over_relaxation = std::max(std::min(over_relaxation, 2.0f), 1.0f);
        ImGui::InputInt("Num Iterations", &iterations, 1);
        iterations = std::max(std::min(iterations, 10), 1);

        ImGui::Checkbox("Use Density Correction", &density_correction);
        ImGui::End();

        // fluid
        fluid->setGravity(grav);
        fluid->setDt(dt);
        fluid->setFlipRatio(flip_ratio);
        fluid->setStiffnessCoeff(stiffness_coefficient);
        fluid->setDensityCorrection(density_correction);
        fluid->setObstacle(lastMouseX, lastMouseY, dx / dt, dy / dt);
        fluid->setOverRelaxation(over_relaxation);
        fluid->setIterations(iterations);
        fluid->update(render_option, lmbDown);

        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(window);

        frameDuration = glfwGetTime() - frameStart;
    }

    // Cleanup
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    delete fluid;
    glfwDestroyWindow(window);
    glfwTerminate();

    return 0;
}