#pragma once
#include <memory>

#include "Eigen/Dense"

#include "jelly.h"
#include "../util/helper.h"

namespace simulation {
class Terrain;
enum class TerrainType : char { Plane, Elevator };

class TerrainFactory final {
 public:
    // no instance, only static usage
    TerrainFactory() = delete;
    static std::unique_ptr<Terrain> CreateTerrain(TerrainType type);
};

// a virtual class
class Terrain {
 public:
    virtual ~Terrain() = default;
    Eigen::Matrix4f getModelMatrix();
    //getter
    float getMass() const;
    Eigen::Vector3f getPosition() const;
    Eigen::Vector3f getVelocity() const;
    Eigen::Vector3f getAcceleration() const;
    Eigen::Vector3f getForce() const;
    //  setter
    void setMass(const float _mass);
    void setPosition(const Eigen::Vector3f &_position);
    void setVelocity(const Eigen::Vector3f &_velocity);
    void setAcceleration(const Eigen::Vector3f &_acceleration);
    void setForce(const Eigen::Vector3f &_force);
    //  method
    void addPosition(const Eigen::Vector3f &_position);
    void addVelocity(const Eigen::Vector3f &_velocity);
    void addAcceleration(const Eigen::Vector3f &_acceleration);
    void addForce(const Eigen::Vector3f &_force);

    virtual void reset() = 0;
    virtual TerrainType getType() = 0;
    virtual void handleCollision(const float delta_T, Jelly& jelly) = 0;

 protected:
    Eigen::Matrix4f modelMatrix = Eigen::Matrix4f::Identity();

    float mass = 1.0f;
    Eigen::Vector3f initialPosition = Eigen::Vector3f::Zero();
    Eigen::Vector3f position = Eigen::Vector3f::Zero();
    Eigen::Vector3f velocity = Eigen::Vector3f::Zero();
    Eigen::Vector3f force = Eigen::Vector3f::Zero();

};

class PlaneTerrain final : public Terrain {
 public:
    friend class TerrainFactory;

    PlaneTerrain();

    void reset();

    TerrainType getType() override;
    void handleCollision(const float delta_T, Jelly& jelly) override;

 private:
    Eigen::Vector3f normal = Eigen::Vector3f(sin(util::PI<float>() / 9), cos(util::PI<float>() / 9), 0);
};

class ElevatorTerrain final : public Terrain {
 public:
    friend class TerrainFactory;

    ElevatorTerrain();

    TerrainType getType() override;
    void handleCollision(const float delta_T, Jelly& jelly) override;
    void reset();
 private:
    Eigen::Vector3f normal = Eigen::Vector3f(0.0f, 1.0f, 0.0f);
};

}  // namespace simulation
