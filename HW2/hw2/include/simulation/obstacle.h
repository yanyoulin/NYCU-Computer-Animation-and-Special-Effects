#pragma once

#include <memory>

#include "Eigen/Core"
#include "graphics/cube.h"

namespace graphics {
class Program;
}

namespace kinematics {
class Obstacle {
 public:
    Obstacle() noexcept;
    // no copy constructor
    Obstacle(const Obstacle&) = delete;
    Obstacle(Obstacle&&) noexcept;
    // You need this for alignment otherwise it may crash
    // Ref: https://eigen.tuxfamily.org/dox/group__TopicStructHavingEigenMembers.html
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Obstacle& operator=(const Obstacle&) = delete;
    Obstacle& operator=(Obstacle&&) noexcept;
    void render(graphics::Program* program);
    Eigen::Vector4d& getCurrentPosition();
    // Calculate and set position
    void setModelMatrix();
    void setCurrentPosition(const Eigen::Vector4d& pos);

 private:
    Eigen::Vector4d current_position = Eigen::Vector4d::Zero();
    std::unique_ptr<graphics::Cube> graphics;
};
}  // namespace kinematics
