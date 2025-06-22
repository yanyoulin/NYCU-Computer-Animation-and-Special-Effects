#pragma once
#include <memory>

#include "buffer.h"
#include "rigidbody.h"

namespace graphics {
class Cube final : public Rigidbody {
 public:
    Cube() noexcept;
    Cube(const Cube&) noexcept;
    Cube(Cube&&) noexcept;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Cube& operator=(const Cube&) noexcept;
    Cube& operator=(Cube&&) noexcept;

    void bindVBO() const override;
    void generateVertices() override;
    void render(Program* shaderProgram) override;

 private:
    std::shared_ptr<Buffer<1, GL_ARRAY_BUFFER>> vbo = nullptr;
    std::shared_ptr<Buffer<1, GL_ELEMENT_ARRAY_BUFFER>> ebo = nullptr;
    static std::weak_ptr<Buffer<1, GL_ARRAY_BUFFER>> vbo_weak;
    static std::weak_ptr<Buffer<1, GL_ELEMENT_ARRAY_BUFFER>> ebo_weak;
    static bool isInitialized;
};
}  // namespace graphics
