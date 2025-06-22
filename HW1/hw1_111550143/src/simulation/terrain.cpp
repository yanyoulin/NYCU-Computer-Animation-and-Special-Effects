#include "terrain.h"

#include <stdexcept>

#include "../util/helper.h"
#include <iostream>

#define PI 3.14159265

namespace simulation {
// Factory
std::unique_ptr<Terrain> TerrainFactory::CreateTerrain(TerrainType type) {
    switch (type) {
        case simulation::TerrainType::Plane:
            return std::make_unique<PlaneTerrain>();
        case simulation::TerrainType::Elevator:
            return std::make_unique<ElevatorTerrain>();

        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}
// Terrain

Eigen::Matrix4f Terrain::getModelMatrix() { return modelMatrix; }

float Terrain::getMass() const { return mass; }

Eigen::Vector3f Terrain::getPosition() const { return position; }

Eigen::Vector3f Terrain::getVelocity() const { return velocity; }

Eigen::Vector3f Terrain::getAcceleration() const { return force / mass; }

Eigen::Vector3f Terrain::getForce() const { return force; }

void Terrain::setMass(const float _mass) { mass = _mass; }

void Terrain::setPosition(const Eigen::Vector3f& _position) { 
    modelMatrix = util::translate(_position - position) * modelMatrix;
    position = _position;
}

void Terrain::setVelocity(const Eigen::Vector3f& _velocity) { velocity = _velocity; }

void Terrain::setAcceleration(const Eigen::Vector3f& _acceleration) { force = _acceleration * mass; }

void Terrain::setForce(const Eigen::Vector3f& _force) { force = _force; }

void Terrain::addPosition(const Eigen::Vector3f& _position) { 
    position += _position;
    modelMatrix = util::translate(_position) * modelMatrix;
}

void Terrain::addVelocity(const Eigen::Vector3f& _velocity) { velocity += _velocity; }

void Terrain::addAcceleration(const Eigen::Vector3f& _acceleration) { force += _acceleration * mass; }

void Terrain::addForce(const Eigen::Vector3f& _force) { force += _force; }

// Note:
// You should update each particles' velocity (base on the equation in
// slide) and force (contact force : resist + friction) in handleCollision function

// PlaneTerrain //

PlaneTerrain::PlaneTerrain() { reset(); }

void PlaneTerrain::reset() { 
    modelMatrix = util::translate(0.0f, position[1], 0.0f) * util::rotateDegree(0, 0, -20) * util::scale(30, 1, 30);
}


TerrainType PlaneTerrain::getType() { return TerrainType::Plane; }

void PlaneTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.5f;
    // TODO#3-1: Handle collision when a particle collide with the plane terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data.
    //   3. The plane spans 30x30 units in the XZ plane and is rotated -20 degrees around the Z-axis.
    // Hint:
    //   1. Review "particles.pptx" from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.

    Eigen::Vector3f normal = Eigen::Vector3f(sinf(PI/9), cosf(PI/9), 0);
    Eigen::Vector3f PointOfPlane = Eigen::Vector3f(0, 0, 0);
    Eigen::Vector3f n = normal.normalized();

    for (int i = 0; i < jelly.getParticleNum(); i++) {
        Particle& particle = jelly.getParticle(i);
        Eigen::Vector3f x = particle.getPosition();
        Eigen::Vector3f v = particle.getVelocity();
        Eigen::Vector3f v_n = v.dot(n) * n;
        Eigen::Vector3f v_t = v - v_n;
        Eigen::Vector3f force = particle.getForce();

        // collision detection
        if(n.dot(x - PointOfPlane) < eEPSILON && v.dot(n) < 0)
        {
            Eigen::Vector3f v_n_new = -1 * coefResist * v_n;
            Eigen::Vector3f v_t_new = v_t;
            particle.setVelocity(v_n_new + v_t_new);

            // contact force
            if(std::abs(n.dot(x - PointOfPlane)) < eEPSILON || std::abs(v.dot(n)) < eEPSILON)  
            {
                if(n.dot(force) < 0){
                    Eigen::Vector3f frictionForce = -1 * coefFriction * (-n.dot(force)) * v_t.normalized();
                    Eigen::Vector3f contactForce = -1 * (n.dot(force)) * n;
                    particle.addForce(frictionForce + contactForce);
                }
            }
        }
    }
}

ElevatorTerrain::ElevatorTerrain() { 
    reset();
}

void ElevatorTerrain::reset() {
    modelMatrix = util::translate(0.0f, 1.0f, 0.0f) * util::rotateDegree(0, 0, 0) * util::scale(5, 1, 5);
    position = Eigen::Vector3f(0.0f, 1.0f, 0.0f);
    velocity = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
}

TerrainType ElevatorTerrain::getType() { return TerrainType::Elevator; }

void ElevatorTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.5f;

    // TODO#3-2: Implement the collision handling between the jelly and the elevator
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data.
    //   3. The elevator plane spans 5x5 units in the XZ plane.
    // Hint:
    //   1. Review "particles.pptx" from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.

    Eigen::Vector3f normal = Eigen::Vector3f(0, 1, 0);
    Eigen::Vector3f PointOfPlane = Terrain::getPosition();
    Eigen::Vector3f n = normal.normalized();
    Eigen::Vector3f elevator_v = Terrain::getVelocity();
    float elevator_m = Terrain::getMass();

    for (int i = 0; i < jelly.getParticleNum(); i++) {
        Particle& particle = jelly.getParticle(i);
        Eigen::Vector3f x = particle.getPosition();
        Eigen::Vector3f v = particle.getVelocity();
        Eigen::Vector3f force = particle.getForce();
        float particle_m = particle.getMass();

        // collision detection
        if(n.dot(x - PointOfPlane) < eEPSILON && (v - elevator_v).dot(n) < 0)
        {
            float u_a = v.dot(n);
            float u_b = elevator_v.dot(n);
            Eigen::Vector3f v_t = v - u_a * n;

            float v_a = (particle_m * u_a + elevator_m * u_b + elevator_m * coefResist * (u_b - u_a)) / (particle_m + elevator_m);

            particle.setVelocity(v_t + v_a * n);

            // contact force
            if(std::abs(n.dot(x - PointOfPlane)) < eEPSILON || std::abs(v.dot(n)) < eEPSILON)  
            {
                if(n.dot(force) < 0){
                    Eigen::Vector3f frictionForce = -1 * coefFriction * (-n.dot(force)) * v_t.normalized();
                    Eigen::Vector3f contactForce = -1 * (n.dot(force)) * n;
                    particle.addForce(frictionForce + contactForce);
                }
            }
        }
        

    }
}

}  // namespace simulation
