#include "graphics/cube.h"

#include "glad/gl.h"
#include "graphics/shader.h"

namespace graphics {

bool Cube::isInitialized = false;
std::weak_ptr<Buffer<1, GL_ARRAY_BUFFER>> Cube::vbo_weak;
std::weak_ptr<Buffer<1, GL_ELEMENT_ARRAY_BUFFER>> Cube::ebo_weak;

Cube::Cube() noexcept {
    if (!(vbo = Cube::vbo_weak.lock())) {
        vbo = std::make_shared<Buffer<1, GL_ARRAY_BUFFER>>();
        Cube::vbo_weak = vbo;
    }
    if (!(ebo = Cube::ebo_weak.lock())) {
        ebo = std::make_shared<Buffer<1, GL_ELEMENT_ARRAY_BUFFER>>();
        Cube::ebo_weak = ebo;
    }
    initialize();
    glBindVertexArray(vao);
    ebo->bind();
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

Cube::Cube(const Cube& other) noexcept : Rigidbody(other), vbo(other.vbo), ebo(other.ebo) {
    initialize();
    glBindVertexArray(vao);
    ebo->bind();
    glBindVertexArray(0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

Cube::Cube(Cube&& other) noexcept
    : Rigidbody(std::forward<Cube>(other)), vbo(std::move(other.vbo)), ebo(std::move(other.ebo)) {}

Cube& Cube::operator=(const Cube& other) noexcept {
    if (this != &other) {
        Rigidbody::operator=(other);
        vbo = other.vbo;
        ebo = other.ebo;
        initialize();
        glBindVertexArray(vao);
        ebo->bind();
        glBindVertexArray(0);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }
    return *this;
}

Cube& Cube::operator=(Cube&& other) noexcept {
    if (this != &other) {
        Rigidbody::operator=(std::forward<Cube>(other));
        vbo = std::move(other.vbo);
        ebo = std::move(other.ebo);
    }
    return *this;
}

void Cube::bindVBO() const { vbo->bind(); }

void Cube::generateVertices() {
    if (!Cube::isInitialized) {
        GLfloat vertices[] = {
            // position        // normal         // texcoord
            // front
            -1, -1, 1, 0, 0, 1, 0, 0, 1, -1, 1, 0, 0, 1, 1, 0, 1, 1, 1, 0, 0, 1, 1, 1, -1, 1, 1, 0, 0, 1, 0, 1,
            // back
            -1, -1, -1, 0, 0, -1, 1, 0, -1, 1, -1, 0, 0, -1, 1, 1, 1, 1, -1, 0, 0, -1, 0, 1, 1, -1, -1, 0, 0, -1, 0, 0,
            // left
            -1, -1, -1, -1, 0, 0, 0, 0, -1, -1, 1, -1, 0, 0, 1, 0, -1, 1, 1, -1, 0, 0, 1, 1, -1, 1, -1, -1, 0, 0, 0, 1,
            // right
            1, -1, -1, 1, 0, 0, 1, 0, 1, 1, -1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, -1, 1, 1, 0, 0, 0, 0,
            // top
            -1, 1, -1, 0, 1, 0, 0, 1, -1, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, -1, 0, 1, 0, 1, 1,
            // bottom
            -1, -1, -1, 0, -1, 0, 0, 0, 1, -1, -1, 0, -1, 0, 1, 0, 1, -1, 1, 0, -1, 0, 1, 1, -1, -1, 1, 0, -1, 0, 0, 1};

        GLuint indices[] = {0,  1,  2,  2,  3,  0,  4,  5,  6,  6,  7,  4,  8,  9,  10, 10, 11, 8,
                            12, 13, 14, 14, 15, 12, 16, 17, 18, 18, 19, 16, 20, 21, 22, 22, 23, 20};

        vbo->bind();
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
        ebo->bind();
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

        Cube::isInitialized = true;
    }
}

void Cube::render(Program* shaderProgram) {
    if (texture) {
        shaderProgram->setUniform("useTexture", 1);
        shaderProgram->setUniform("diffuseTexture", texture->getIndex());
    } else {
        shaderProgram->setUniform("useTexture", 0);
        shaderProgram->setUniform("baseColor", baseColor);
    }
    shaderProgram->setUniform("model", modelMatrix);
    shaderProgram->setUniform("invtransmodel", inverseTransposeModel);

    glBindVertexArray(vao);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, nullptr);
    glBindVertexArray(0);
}

}  // namespace graphics
