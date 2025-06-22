#include "simulation/obstacle.h"

#include <utility>

namespace kinematics {
Obstacle::Obstacle() noexcept
    : current_position(0.0f, 1.0f, 0.0f, 0.0), graphics(std::make_unique<graphics::Cube>()) {
    auto texture_folder = util::PathFinder::find("Texture");
    auto slime = std::make_shared<graphics::Texture>(texture_folder / "slime.png");
    graphics->setTexture(std::move(slime));
}

Obstacle::Obstacle(Obstacle&& other) noexcept
    : current_position(std::move(other.current_position)), graphics(std::move(other.graphics)) {}

Obstacle& Obstacle::operator=(Obstacle&& other) noexcept {
    if (this != &other) {
        current_position = std::move(other.current_position);
        graphics = std::move(other.graphics);
    }
    return *this;
}

void Obstacle::render(graphics::Program* program) { graphics->render(program); }

Eigen::Vector4d& Obstacle::getCurrentPosition() { return current_position; }

void Obstacle::setModelMatrix() {
    Eigen::Affine3d trans = Eigen::Affine3d::Identity();
    trans.translate(current_position.head<3>());
    trans.scale(Eigen::Vector3d(1, 1, 1));
    graphics->setModelMatrix(trans.cast<float>());
}

void Obstacle::setCurrentPosition(const Eigen::Vector4d& pos) { current_position = pos; }
}  // namespace kinematics
