#include "simulation/kinematics.h"

#include <iostream>
#include "Eigen/Dense"
#include "acclaim/bone.h"
#include "util/helper.h"


namespace kinematics {
void dfs(acclaim::Bone* bone, const acclaim::Posture& posture, std::vector<bool>& visited) {
    if (visited[bone->idx]) {
        return;
    }

    visited[bone->idx] = true;

    if (bone->parent) {
        bone->start_position = bone->parent->end_position;
        bone->rotation = bone->parent->rotation * bone->rot_parent_current * util::rotateDegreeZYX(posture.bone_rotations[bone->idx]);
        bone->end_position = bone->start_position + bone->rotation * bone->dir.normalized() * bone->length;
    }
    
    if (bone->child && !visited[bone->child->idx]) {
        dfs(bone->child, posture, visited);
    }
    
    for (acclaim::Bone* sibling = bone->sibling; sibling != nullptr; sibling = sibling->sibling) {
        if (!visited[sibling->idx]) {
            dfs(sibling, posture, visited);
        }
    }
}

void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
    // TODO#1: Forward Kinematic
    // Hint:
    // - Traverse the skeleton tree from root to leaves.
    // - Compute each bone's global rotation and global position.
    // - Use local rotation (from posture) and bone hierarchy (parent rotation, offset, etc).
    // - Remember to update both bone->start_position and bone->end_position.
    // - Use bone->rotation to store global rotation (after combining parent, local, etc).
    std::vector<bool> visited(31, false);
    bone->start_position = posture.bone_translations[0];
    bone->rotation = bone->rot_parent_current * util::rotateDegreeZYX(posture.bone_rotations[bone->idx]);
    bone->end_position = bone->start_position + bone->rotation * bone->dir.normalized() * bone->length;
    visited[bone->idx] = true;
    if (bone->child) {
        dfs(bone->child, posture, visited);
    }
}

Eigen::VectorXd pseudoInverseLinearSolver(const Eigen::Matrix4Xd& Jacobian, const Eigen::Vector4d& target) {
    Eigen::VectorXd deltatheta;
    // TODO#2: Inverse linear solver (find x which min(| jacobian * x - target |))
    // Hint:
    //   1. Linear algebra - least squares solution
    //   2. https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Construction
    // Note:
    //   1. SVD or other pseudo-inverse method is useful
    //   2. Some of them have some limitation, if you use that method you should check it.
    Eigen::JacobiSVD<Eigen::Matrix4Xd> svd(Jacobian, Eigen::ComputeThinU | Eigen::ComputeThinV);
    deltatheta = svd.solve(target);

    return deltatheta;
}

/**
 * @brief Perform inverse kinematics (IK)
 *
 * @param target_pos The position where `end_bone` will move to.
 * @param obs_pos The position where the obstacle is at
 * @param obsActive Whether the obstacle is active or not
 * @param start_bone This bone is the last bone you can move while doing IK
 * @param end_bone This bone will try to reach `target_pos`
 * @param posture The original AMC motion's reference, you need to modify this
 *
 * @return True if IK is stable (HW2 bonus)
 */
bool inverseJacobianIKSolver(const Eigen::Vector4d& target_pos, const Eigen::Vector4d& obs_pos, bool obsActive,
                             acclaim::Bone* start_bone, acclaim::Bone* end_bone, acclaim::Posture& posture) {
    constexpr int max_iteration = 1000;
    constexpr double epsilon = 1E-10;
    constexpr double step = 0.1;
    constexpr double obsAvoidThreshold = 1.01;  // if bone is within 1 unit from obstacle 0.01

    // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is the
    // root.
    acclaim::Bone* root_bone = start_bone - start_bone->idx;

    // TODO#3:
    // Perform inverse kinematics (IK)
    // HINTs will tell you what should do in that area.
    // Of course you can ignore it (Any code below this line) and write your own code.
    acclaim::Posture original_posture(posture);

    size_t bone_num = 0;
    std::vector<acclaim::Bone*> boneList;
    acclaim::Bone* current = end_bone;
    // TODO#3-1:
    // Calculate number of bones need to move to perform IK, store in `bone_num`
    // (a.k.a. how may bones from end_bone to its parent than to start_bone (include both side))
    // Store the bones need to move to perform IK into boneList
    // Hint:
    //   1. Traverse from end_bone to start_bone is easier than start to end (since there is only 1 parent)
    //   2. If start bone is not reachable from end. Go to root first.
    // Note:
    //   1. Both start and end should be in the list

    // 3-1 start
    while (current != nullptr) {
        boneList.push_back(current);
        bone_num++;
        if (current == start_bone) break;
        current = current->parent;
    }
    // 3-1 end

    Eigen::Matrix4Xd Jacobian(4, 3 * bone_num);
    Jacobian.setZero();
    for (int iter = 0; iter < max_iteration; ++iter) {
        forwardSolver(posture, root_bone);
        Eigen::Vector4d desiredVector = target_pos - end_bone->end_position;

        if (desiredVector.norm() < epsilon) {
            break;
        }
        // TODO#3-2 (compute jacobian)
        //   1. Compute arm vectors
        //   2. Compute jacobian columns, store in `Jacobian`
        // Hint:
        //   1. You should not put rotation in jacobian if it doesn't have that DoF.
        //   2. jacobian.col(/* some column index */) = /* jacobian column */

        // 3-2 start
        Jacobian.setZero();
        for (int i = 0; i < bone_num; i++) {
            acclaim::Bone* bone = boneList[i];
            Eigen::Vector3d arm = end_bone->end_position.head<3>() - bone->start_position.head<3>();

            if (bone->dofrx) {
                Jacobian.col(i * 3).head<3>() = ((bone->rotation * Eigen::Vector4d(1, 0, 0, 0)).head<3>()).cross(arm);
                Jacobian.col(i * 3)(3) = 0;
            }
            if (bone->dofry) {
                Jacobian.col(i * 3 + 1).head<3>() = ((bone->rotation * Eigen::Vector4d(0, 1, 0, 0)).head<3>()).cross(arm);
                Jacobian.col(i * 3 + 1)(3) = 0;
            }
            if (bone->dofrz) {
                Jacobian.col(i * 3 + 2).head<3>() = ((bone->rotation * Eigen::Vector4d(0, 0, 1, 0)).head<3>()).cross(arm);
                Jacobian.col(i * 3 + 2)(3) = 0;
            }
        }
        // 3-2 end

        // TODO#3-3 (obstacle avoidance)
        //  1. Iterate through all bones in `boneList`.
        //  2. Compute the center of each bone (average of start and end positions).
        //  3. Calculate the vector from obstacle center to bone center.
        //  4. If distance is below threshold, compute repulsive vector.
        //  5. Add this repulsive vector to `desiredVector`.
        // Hint:
        // - Use a constant threshold distance to determine proximity.
        // - The repulsive vector should point away from the obstacle.
        // - Use `.head<3>().norm()` to compute 3D distance from a 4D vector.
        // - Normalize the repulsive vector and scale it based on how close it is.
        if (obsActive) {
            for (int i = 0; i < bone_num; i++) {
                acclaim::Bone* bone = boneList[i];
                Eigen::Vector4d center = (bone->start_position + bone->end_position) / 2;
                Eigen::Vector4d obs_vector = center - obs_pos;
                double dist = obs_vector.head<3>().norm();

                if (dist < obsAvoidThreshold && dist > 1e-6) {
                    Eigen::Vector3d d = obs_vector.head<3>().normalized();
                    desiredVector.head<3>() += d * (obsAvoidThreshold - dist);
                }
            }
        }

        Eigen::VectorXd deltatheta = step * pseudoInverseLinearSolver(Jacobian, desiredVector);
        // TODO#3-4 (update rotation)
        //   Update `posture.bone_rotation` (in euler angle / degrees) using deltaTheta
        // Hint:
        //   1. You can ignore rotation limit of the bone.
        // Bonus:
        //   1. You cannot ignore rotation limit of the bone.

        // 3-4 start
        double PI = 3.14159265358979323846;
        for (int i = 0; i < bone_num; i++) {
            acclaim::Bone* bone = boneList[i];
            if (bone->dofrx) {
                posture.bone_rotations[bone->idx][0] += deltatheta(i * 3) * 180 / PI;
            }
            if (bone->dofry) {
                posture.bone_rotations[bone->idx][1] += deltatheta(i * 3 + 1) * 180 / PI;
            }
            if (bone->dofrz) {
                posture.bone_rotations[bone->idx][2] += deltatheta(i * 3 + 2) * 180 / PI;
            }

            // Bonus
            if (bone->dofrx) {
                if (posture.bone_rotations[bone->idx][0] < bone->rxmin) {
                    posture.bone_rotations[bone->idx][0] = bone->rxmin;
                } else if (posture.bone_rotations[bone->idx][0] > bone->rxmax) {
                    posture.bone_rotations[bone->idx][0] = bone->rxmax;
                }
            }
            if (bone->dofry) {
                if (posture.bone_rotations[bone->idx][1] < bone->rymin) {
                    posture.bone_rotations[bone->idx][1] = bone->rymin;
                } else if (posture.bone_rotations[bone->idx][1] > bone->rymax) {
                    posture.bone_rotations[bone->idx][1] = bone->rymax;
                }
            }
            if (bone->dofrz) {
                if (posture.bone_rotations[bone->idx][2] < bone->rzmin) {
                    posture.bone_rotations[bone->idx][2] = bone->rzmin;
                } else if (posture.bone_rotations[bone->idx][2] > bone->rzmax) {
                    posture.bone_rotations[bone->idx][2] = bone->rzmax;
                }
            }
            
        }
        // 3-4 end
    }
    // TODO#3-5
    // Return whether IK is stable
    // i.e. whether the ball is reachable
    // Hint:
    //      1. comment out the line here and return whether the IK is stable or not
    //      2. if the ball is reachable,  swinging its hand in air

    // 3-5 start
    forwardSolver(posture, root_bone);
    Eigen::Vector4d finalDesiredVector = target_pos - end_bone->end_position;
    if (finalDesiredVector.head<3>().norm() < epsilon) {
        return true;
    } else {
        posture = original_posture;
        forwardSolver(posture, root_bone);
        return false;
    }
    // 3-5 end
}

}  // namespace kinematics
