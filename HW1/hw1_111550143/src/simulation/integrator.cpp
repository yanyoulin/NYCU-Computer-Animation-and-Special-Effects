#include "integrator.h"
#include <iostream>

#include <vector>
namespace simulation {
// Factory
std::unique_ptr<Integrator> IntegratorFactory::CreateIntegrator(IntegratorType type) {
    switch (type) {
        case simulation::IntegratorType::ExplicitEuler:
            return std::make_unique<ExplicitEulerIntegrator>();
        case simulation::IntegratorType::ImplicitEuler:
            return std::make_unique<ImplicitEulerIntegrator>();
        case simulation::IntegratorType::MidpointEuler:
            return std::make_unique<MidpointEulerIntegrator>();
        case simulation::IntegratorType::RungeKuttaFourth:
            return std::make_unique<RungeKuttaFourthIntegrator>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}

//
// ExplicitEulerIntegrator
//

IntegratorType ExplicitEulerIntegrator::getType() { return IntegratorType::ExplicitEuler; }

void ExplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-1: Integrate position and velocity
    //   1. Integrate position using current velocity.
    //   2. Integrate velocity using current acceleration.
    //   3. Clear force
    // Note:
    //   1. You should do this first because it is very simple. Then you can check whether your collision is correct or
    //   not.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Review "ODE_basics.pptx" from p.15 - p.16
    float dt = particleSystem.deltaTime;

    for (int i = 0; i < particleSystem.getJellyCount(); i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            Eigen::Vector3f acceleration = particle.getAcceleration();
            Eigen::Vector3f velocity = particle.getVelocity();
            Eigen::Vector3f position = particle.getPosition();

            
            position = position + velocity * dt;
            velocity = velocity + acceleration * dt;

            particle.setPosition(position);
            particle.setVelocity(velocity);
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }

    if (particleSystem.sceneIdx == 1 && particleSystem.elevatorTerrain != nullptr) {
        Terrain& elevator = *(particleSystem.elevatorTerrain);

        Eigen::Vector3f a = elevator.getAcceleration();
        Eigen::Vector3f v = elevator.getVelocity();
        Eigen::Vector3f x = elevator.getPosition();

        Eigen::Vector3f v_new = v + dt * a;
        Eigen::Vector3f x_new = x + dt * v;

        elevator.setVelocity(v_new);
        elevator.setPosition(x_new);
        elevator.setForce(Eigen::Vector3f::Zero());
    }

}

//
// ImplicitEulerIntegrator
//

IntegratorType ImplicitEulerIntegrator::getType() { return IntegratorType::ImplicitEuler; }

void ImplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-2: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 and Vn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce`and `MassSpringSystem::computeElevatorForce`
    //      with modified position and velocity to get Xn+1.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Remember to reset particleSystem.elevatorCounter back to the original state
    //   5. Review "ODE_implicit.pptx" from p.18 - p.19
    float dt = particleSystem.deltaTime;
    int jellyCount = particleSystem.getJellyCount();
    int elevatorCounterBU = particleSystem.elevatorCounter;

    std::vector<std::vector<Eigen::Vector3f>> oldPos(jellyCount);
    std::vector<std::vector<Eigen::Vector3f>> oldVel(jellyCount);
    Eigen::Vector3f ElevatorPos = particleSystem.elevatorTerrain->getPosition();
    Eigen::Vector3f ElevatorVel = particleSystem.elevatorTerrain->getVelocity();

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        int particleNum = jelly->getParticleNum();
        oldPos[i].resize(particleNum);
        oldVel[i].resize(particleNum);
        for (int j = 0; j < particleNum; j++) {
            Particle& particle = jelly->getParticle(j);
            oldPos[i][j] = particle.getPosition();
            oldVel[i][j] = particle.getVelocity();
            particle.setPosition(particle.getPosition());
            particle.setVelocity(particle.getVelocity());
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        particleSystem.computeJellyForce(*jelly);
    }

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for(int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            Eigen::Vector3f acceleration = particle.getAcceleration();
            Eigen::Vector3f v_new = oldVel[i][j] + acceleration * dt;
            Eigen::Vector3f x_new = oldPos[i][j] + particle.getVelocity() * dt;

            particle.setPosition(x_new);
            particle.setVelocity(v_new);
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }

    if (particleSystem.sceneIdx == 1 && particleSystem.elevatorTerrain != nullptr) {
        Terrain& elevator = *(particleSystem.elevatorTerrain);
    
        elevator.setVelocity(ElevatorVel);
        elevator.setPosition(ElevatorPos);
        elevator.setForce(Eigen::Vector3f::Zero());
    
        particleSystem.computeElevatorForce();
    
        Eigen::Vector3f v_new = ElevatorVel + elevator.getAcceleration() * dt;
        Eigen::Vector3f x_new = ElevatorPos + elevator.getVelocity() * dt;
    
        particleSystem.elevatorCounter = elevatorCounterBU;

        elevator.setVelocity(v_new);
        elevator.setPosition(x_new);
        elevator.setForce(Eigen::Vector3f::Zero());
    }


}

//
// MidpointEulerIntegrator
//

IntegratorType MidpointEulerIntegrator::getType() { return IntegratorType::MidpointEuler; }

void MidpointEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-3: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce`and `MassSpringSystem::computeElevatorForce`
    //      with modified position and velocity to get Xn+1.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Remember to reset particleSystem.elevatorCounter back to the original state
    //   5. Review "ODE_basics.pptx" from p.18 - p.19
    float dt = particleSystem.deltaTime;
    int jellyCount = particleSystem.getJellyCount();
    int elevatorCounterBU = particleSystem.elevatorCounter;

    std::vector<std::vector<Eigen::Vector3f>> oldPos(jellyCount);
    std::vector<std::vector<Eigen::Vector3f>> oldVel(jellyCount);
    Eigen::Vector3f ElevatorPos = particleSystem.elevatorTerrain->getPosition();
    Eigen::Vector3f ElevatorVel = particleSystem.elevatorTerrain->getVelocity();

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        int particleNum = jelly->getParticleNum();
        oldPos[i].resize(particleNum);
        oldVel[i].resize(particleNum);
        for (int j = 0; j < particleNum; j++) {
            Particle& particle = jelly->getParticle(j);
            oldPos[i][j] = particle.getPosition();
            oldVel[i][j] = particle.getVelocity();
            Eigen::Vector3f acceleration = particle.getAcceleration();
            Eigen::Vector3f v_mid = particle.getVelocity() + acceleration * dt / 2.0f;
            Eigen::Vector3f x_mid = particle.getPosition() + particle.getVelocity() * dt / 2.0f;

            particle.setPosition(x_mid);
            particle.setVelocity(v_mid);
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        particleSystem.computeJellyForce(*jelly);
    }

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for(int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            Eigen::Vector3f acceleration = particle.getAcceleration();
            Eigen::Vector3f v_new = oldVel[i][j] + acceleration * dt;
            Eigen::Vector3f x_new = oldPos[i][j] + particle.getVelocity() * dt;

            particle.setPosition(x_new);
            particle.setVelocity(v_new);
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }

    if (particleSystem.sceneIdx == 1 && particleSystem.elevatorTerrain != nullptr) {
        Terrain& elevator = *(particleSystem.elevatorTerrain);
    
        Eigen::Vector3f a0 = elevator.getAcceleration();
        Eigen::Vector3f v_mid = ElevatorVel + a0 * (dt * 0.5f);
        Eigen::Vector3f x_mid = ElevatorPos + ElevatorVel * (dt * 0.5f);
    
        elevator.setVelocity(v_mid);
        elevator.setPosition(x_mid);
        elevator.setForce(Eigen::Vector3f::Zero());
    
        particleSystem.computeElevatorForce();
    
        Eigen::Vector3f a_mid = elevator.getAcceleration();
        Eigen::Vector3f v_new = ElevatorVel + a_mid * dt;
        Eigen::Vector3f x_new = ElevatorPos + v_mid * dt;
    
        particleSystem.elevatorCounter = elevatorCounterBU;
        elevator.setVelocity(v_new);
        elevator.setPosition(x_new);
        elevator.setForce(Eigen::Vector3f::Zero());
    }


}

//
// RungeKuttaFourthIntegrator
//

IntegratorType RungeKuttaFourthIntegrator::getType() { return IntegratorType::RungeKuttaFourth; }

void RungeKuttaFourthIntegrator::integrate(MassSpringSystem& particleSystem) {
    struct StateStep {
        Eigen::Vector3f deltaVel;
        Eigen::Vector3f deltaPos;
    };
    // TODO#4-4: Integrate velocity and acceleration
    //   1. Backup original particles' data.
    //   2. Compute k1, k2, k3, k4
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce`and `MassSpringSystem::computeElevatorForce`
    //      with modified position and velocity to get Xn+1.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Remember to reset particleSystem.elevatorCounter back to the original state
    //   5. StateStep struct is just a hint, you can use whatever you want.
    //   6. Review "ODE_basics.pptx" from p.21
    float dt = particleSystem.deltaTime;
    int jellyCount = particleSystem.getJellyCount();
    int elevatorCounterBU = particleSystem.elevatorCounter;

    std::vector<std::vector<Eigen::Vector3f>> oldPos(jellyCount);
    std::vector<std::vector<Eigen::Vector3f>> oldVel(jellyCount);
    std::vector<std::vector<StateStep>> k1(jellyCount);
    std::vector<std::vector<StateStep>> k2(jellyCount);
    std::vector<std::vector<StateStep>> k3(jellyCount);
    std::vector<std::vector<StateStep>> k4(jellyCount);
    StateStep k1_elevator, k2_elevator, k3_elevator, k4_elevator;

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        int particleNum = jelly->getParticleNum();
        oldPos[i].resize(particleNum);
        oldVel[i].resize(particleNum);
        k1[i].resize(particleNum);
        k2[i].resize(particleNum);
        k3[i].resize(particleNum);
        k4[i].resize(particleNum);
        for (int j = 0; j < particleNum; j++) {
            Particle& particle = jelly->getParticle(j);
            oldPos[i][j] = particle.getPosition();
            oldVel[i][j] = particle.getVelocity();
        }
    }

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            Eigen::Vector3f acceleration = particle.getAcceleration();
            Eigen::Vector3f velocity = particle.getVelocity();
            Eigen::Vector3f position = particle.getPosition();

            k1[i][j].deltaVel = acceleration * dt;
            k1[i][j].deltaPos = velocity * dt;

            particle.setPosition(oldPos[i][j] + k1[i][j].deltaPos * 0.5f);
            particle.setVelocity(oldVel[i][j] + k1[i][j].deltaVel * 0.5f);
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        particleSystem.computeJellyForce(*jelly);
    }

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            Eigen::Vector3f acceleration = particle.getAcceleration();
            Eigen::Vector3f velocity = particle.getVelocity();
            Eigen::Vector3f position = particle.getPosition();

            k2[i][j].deltaVel = acceleration * dt;
            k2[i][j].deltaPos = velocity * dt;

            particle.setPosition(oldPos[i][j] + k2[i][j].deltaPos * 0.5f);
            particle.setVelocity(oldVel[i][j] + k2[i][j].deltaVel * 0.5f);
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        particleSystem.computeJellyForce(*jelly);
    }

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            Eigen::Vector3f acceleration = particle.getAcceleration();
            Eigen::Vector3f velocity = particle.getVelocity();
            Eigen::Vector3f position = particle.getPosition();

            k3[i][j].deltaVel = acceleration * dt;
            k3[i][j].deltaPos = velocity * dt;

            particle.setPosition(oldPos[i][j] + k3[i][j].deltaPos);
            particle.setVelocity(oldVel[i][j] + k3[i][j].deltaVel);
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        particleSystem.computeJellyForce(*jelly);
    }

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);
            Eigen::Vector3f acceleration = particle.getAcceleration();
            Eigen::Vector3f velocity = particle.getVelocity();
            Eigen::Vector3f position = particle.getPosition();

            k4[i][j].deltaVel = acceleration * dt;
            k4[i][j].deltaPos = velocity * dt;
        }
    }

    for(int i = 0; i < jellyCount; i++) {
        Jelly* jelly = particleSystem.getJellyPointer(i);
        for (int j = 0; j < jelly->getParticleNum(); j++) {
            Particle& particle = jelly->getParticle(j);

            Eigen::Vector3f x_new = oldPos[i][j] + (k1[i][j].deltaPos + 2 * k2[i][j].deltaPos + 2 * k3[i][j].deltaPos + k4[i][j].deltaPos) / 6.0f;
            Eigen::Vector3f v_new = oldVel[i][j] + (k1[i][j].deltaVel + 2 * k2[i][j].deltaVel + 2 * k3[i][j].deltaVel + k4[i][j].deltaVel) / 6.0f;

            particle.setPosition(x_new);
            particle.setVelocity(v_new);
            particle.setForce(Eigen::Vector3f::Zero());
        }
    }

    if (particleSystem.sceneIdx == 1 && particleSystem.elevatorTerrain != nullptr) {
        Terrain& elevator = *(particleSystem.elevatorTerrain);
        
        Eigen::Vector3f a0 = elevator.getAcceleration();
        Eigen::Vector3f v0 = elevator.getVelocity();
        Eigen::Vector3f x0 = elevator.getPosition();

        k1_elevator.deltaVel = a0 * dt;
        k1_elevator.deltaPos = v0 * dt;

        elevator.setPosition(x0 + k1_elevator.deltaPos * 0.5f);
        elevator.setVelocity(v0 + k1_elevator.deltaVel * 0.5f);
        elevator.setForce(Eigen::Vector3f::Zero());

        particleSystem.computeElevatorForce();

        Eigen::Vector3f a1 = elevator.getAcceleration();
        Eigen::Vector3f v1 = elevator.getVelocity();

        k2_elevator.deltaVel = a1 * dt;
        k2_elevator.deltaPos = v1 * dt;

        elevator.setPosition(x0 + k2_elevator.deltaPos * 0.5f);
        elevator.setVelocity(v0 + k2_elevator.deltaVel * 0.5f);
        elevator.setForce(Eigen::Vector3f::Zero());

        particleSystem.computeElevatorForce();

        Eigen::Vector3f a2 = elevator.getAcceleration();
        Eigen::Vector3f v2 = elevator.getVelocity();
        
        k3_elevator.deltaVel = a2 * dt;
        k3_elevator.deltaPos = v2 * dt;

        elevator.setPosition(x0 + k3_elevator.deltaPos);
        elevator.setVelocity(v0 + k3_elevator.deltaVel);
        elevator.setForce(Eigen::Vector3f::Zero());

        particleSystem.computeElevatorForce();

        Eigen::Vector3f a3 = elevator.getAcceleration();
        Eigen::Vector3f v3 = elevator.getVelocity();

        k4_elevator.deltaVel = a3 * dt;
        k4_elevator.deltaPos = v3 * dt;

        particleSystem.elevatorCounter = elevatorCounterBU;
        elevator.setPosition(x0 + (k1_elevator.deltaPos + 2 * k2_elevator.deltaPos + 2 * k3_elevator.deltaPos + k4_elevator.deltaPos) / 6.0f);
        elevator.setVelocity(v0 + (k1_elevator.deltaVel + 2 * k2_elevator.deltaVel + 2 * k3_elevator.deltaVel + k4_elevator.deltaVel) / 6.0f);
        elevator.setForce(Eigen::Vector3f::Zero());
    }


}
}  // namespace simulation
