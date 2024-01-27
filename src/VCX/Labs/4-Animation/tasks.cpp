#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <spdlog/spdlog.h>
#include <iostream>
#include "Labs/4-Animation/tasks.h"
#include "IKSystem.h"
#include "CustomFunc.inl"


namespace VCX::Labs::Animation 
{
    void ForwardKinematics(IKSystem & ik, int StartIndex)
    {
        if (StartIndex == 0)
        {
            ik.JointGlobalRotation[0] = ik.JointLocalRotation[0];
            ik.JointGlobalPosition[0] = ik.JointLocalOffset[0];
            StartIndex                = 1;
        }
        
        for (int i = StartIndex; i < ik.JointLocalOffset.size(); i++) 
        {
            //先计算缩放，然后计算旋转，最后计算平移
            ik.JointGlobalRotation[i] = ik.JointLocalRotation[i] * ik.JointGlobalRotation[i-1];

            glm::vec3 position = ik.JointGlobalRotation[i-1] * ik.JointLocalOffset[i];
            ik.JointGlobalPosition[i] = ik.JointGlobalPosition[i-1] + position;
        }
    }

    void InverseKinematicsCCD(IKSystem & ik, const glm::vec3 & EndPosition, int maxCCDIKIteration, float eps) 
    {
        ForwardKinematics(ik, 0);
        
        int size = ik.NumJoints();
        if (size == 0) 
        {
            return;
        }
        int last = size - 1;
        glm::vec3 goal = EndPosition;
        
        // These functions will be useful: glm::normalize, glm::rotation, glm::quat * glm::quat
        for (int CCDIKIteration = 0; CCDIKIteration < maxCCDIKIteration && glm::l2Norm(ik.EndEffectorPosition() - EndPosition) > eps; CCDIKIteration++) {
            
            //每一步迭代就是将当前节点与末端点之间的连线方向旋转到当前节点与目标连线的方向，从倒数第二个节点开始旋转，因为最后一个节点旋转不产生影响
            for (int j = (int)size - 2; j >= 0; --j)
            {
                //获得末端节点的世界坐标位置
                glm::vec3& effector = ik.JointGlobalPosition[last];

                //获得当前节点的世界变换
                glm::vec3 position = ik.JointGlobalPosition[j];
                glm::quat rotation = ik.JointGlobalRotation[j];

                glm::vec3 toEffector = glm::normalize(effector - position);    //由当前节点指向末端节点
                glm::vec3 toGoal = glm::normalize(goal - position);            //由当前节点指向目标节点

                //计算从toEffector方向到toGoal方向的旋转四元数
                glm::quat effectorToGoal;
                if (glm::length2(toGoal) > 0.00001f) 
                {
                    effectorToGoal = glm::rotation(toEffector, toGoal);
                }

                //更新当前节点的局部旋转
                glm::quat worldRotated = rotation * effectorToGoal;
                glm::quat localRotate = worldRotated * glm::inverse(rotation);
                ik.JointLocalRotation[j] = localRotate * ik.JointLocalRotation[j];
                
                //更新当前节点到末端节点的全局位置和旋转
                ForwardKinematics(ik, j);

                //判断末端节点和目标位置是否足够近
                effector = ik.JointGlobalPosition[last];
                if (glm::length2(goal - effector) < eps)
                {
                    printf("CCD IK的迭代次数 = %d\n", CCDIKIteration + 1);
                    return;
                }
            }
        }
    }

void WorldToIKChain(IKSystem & ik, const std::vector<glm::vec3>& worldPosition)
{
    int size = ik.NumJoints();
    if (size == 0)
    {
        return;
    }
    
    typedef glm::vec<3, double, glm::defaultp> vec3d;
    typedef glm::qua<double, glm::defaultp> quatd;

    for (int i = 0; i < size - 1; ++i)
    {
        vec3d position = ik.JointGlobalPosition[i];
        quatd rotation = ik.JointGlobalRotation[i];

        vec3d toNext = glm::normalize(vec3d(ik.JointGlobalPosition[i + 1]) - position);
        toNext = glm::normalize(glm::inverse(rotation) * toNext);

        vec3d toDesired = glm::normalize(vec3d(worldPosition[i + 1]) - position);
        toDesired = glm::normalize(glm::inverse(rotation) * toDesired);

        quatd delta = glm::rotation(toNext, toDesired);
        ik.JointLocalRotation[i] = delta * quatd(ik.JointLocalRotation[i]);
    }
}

    void InverseKinematicsFABR(IKSystem & ik, const glm::vec3 & EndPosition, int maxFABRIKIteration, float eps) 
    {
        int nJoints = ik.NumJoints();
        if (0 == nJoints)
        {
            return;
        }
        
        ForwardKinematics(ik, 0);
        
        glm::vec3 goal = EndPosition;
        glm::vec3 base = ik.JointGlobalPosition[0];
        
        //初始化存储中间状态的全局位置向量
        std::vector<glm::vec3> global_positions(nJoints);
        memcpy(global_positions.data(), ik.JointGlobalPosition.data(), sizeof(glm::vec3) * nJoints);
        
        int iterCount = 0;
        
        for (int IKIteration = 0; IKIteration < maxFABRIKIteration; IKIteration++)
        {
            // 如果末端节点和目标点非常接近，退出循环
            if (glm::l2Norm(global_positions[nJoints - 1] - goal) <= eps)
            {
                break;
            }
            
            // backward update，将末端节点拉到目标点的位置
            global_positions[nJoints - 1] = goal;
            for (int i = nJoints - 2; i >= 0; i--) 
            {
                //从子节点指向父节点的方向
                glm::vec3 direction = glm::normalize(global_positions[i] - global_positions[i + 1]);
                glm::vec3 offset = direction * ik.JointOffsetLength[i + 1];
                
                //将当前节点往子节点方向拉
                global_positions[i] = global_positions[i + 1] + offset;
            }

            // forward update，将起始节点拉到本身节点的位置
            global_positions[0] = base;
            for (int i = 1; i < nJoints; i++)
            {
                //从父节点向子节点的方向
                glm::vec3 direction = glm::normalize(global_positions[i] - global_positions[i - 1]);
                
                //沿着从父节点向子节点的方向的偏移向量
                glm::vec3 offset = direction * ik.JointOffsetLength[i];
                
                //父节点+偏移向量 即更新其节点的位置
                global_positions[i] = global_positions[i - 1] + offset;
            }
            
            iterCount ++;
        }
        
        printf("FABR IK 的迭代次数 = %d\n", iterCount);
        
        ik.JointGlobalPosition = global_positions;

        // 根据全局位置计算关节的旋转
        for (int i = 0; i < nJoints - 1; i++)
        {
            ik.JointGlobalRotation[i] = glm::rotation(glm::normalize(ik.JointLocalOffset[i + 1]), glm::normalize(ik.JointGlobalPosition[i + 1] - ik.JointGlobalPosition[i]));
        }
        
        //计算每个关节的局部旋转
        ik.JointLocalRotation[0] = ik.JointGlobalRotation[0];
        for (int i = 1; i < nJoints - 1; i++) 
        {
            ik.JointLocalRotation[i] = glm::inverse(ik.JointGlobalRotation[i - 1]) * ik.JointGlobalRotation[i];
        }
        
//        for (int i = 0; i < nJoints; i ++)
//        {
//            printf("%d ik point x = %f, y = %f, z = %f,  global point x = %f, y = %f, z= %f\n", i, 
//                   ik.JointGlobalPosition[i].x, ik.JointGlobalPosition[i].y, ik.JointGlobalPosition[i].z,
//                   global_positions[i].x, global_positions[i].y, global_positions[i].z);
//        }
//        
//        printf("target point x = %f, y = %f, z = %f\n",
//               goal.x, goal.y, goal.z);
    }

    void InverseKinematicsJacobianTranspose(IKSystem & ik, const glm::vec3 & EndPosition, int maxIteration, float eps)
    {
        ForwardKinematics(ik, 0);
        
        int size = ik.NumJoints();
        if (size == 0)
        {
            return;
        }
        
        int iteration = 0;
        // These functions will be useful: glm::normalize, glm::rotation, glm::quat * glm::quat
        for (; iteration < maxIteration && glm::l2Norm(ik.EndEffectorPosition() - EndPosition) > eps; iteration++)
        {
            //step 1 计算出delta x，即当前末端点减去目标点的向量
            glm::vec3 delta = ik.EndEffectorPosition() - EndPosition;
            
            Eigen::MatrixXf deltaX(3, 1);
            deltaX(0, 0) = delta.x;
            deltaX(1, 0) = delta.y;
            deltaX(2, 0) = delta.z;
            
            //step 2: 计算雅可比矩阵
            
            Eigen::MatrixXf jacobianTransposeMatrix(3 * ik.NumJoints(), 3);
            Eigen::MatrixXf deltaTheta(3 * ik.NumJoints(), 1);
            
            for (int i = 0; i < ik.NumJoints(); i ++)
            {
                // 1 获得父亲关节的全局旋转
                int parentIndex = i - 1;
                if (0 == i)
                {
                    parentIndex = 0;
                }
                const glm::quat& parentRotate = ik.JointGlobalRotation[parentIndex];
                
                // 2 获得当前节点的局部旋转，并转换为XYZ的转角顺序的欧拉角
                const glm::quat& localRatate = ik.JointLocalRotation[i];
                glm::vec3 eulerAngles = glm::eulerAngles(localRatate);
                
                // 3 计算当前节点在全局空间中的朝向
#if 1
                glm::quat rotateX = glm::angleAxis(eulerAngles.x, glm::vec3(1.0, 0.0, 0.0));
                glm::vec3 axisX = glm::normalize(parentRotate * glm::vec3(1.0, 0.0, 0.0));
                glm::vec3 axisY = glm::normalize(parentRotate * rotateX * glm::vec3(0.0, 1.0, 0.0));
                
                glm::quat rotateY = glm::angleAxis(eulerAngles.y, glm::vec3(0.0, 1.0, 0.0));
                glm::vec3 axisZ = glm::normalize(parentRotate * rotateX * rotateY * glm::vec3(0.0, 0.0, 1.0));
#else
                glm::quat rotateZ = glm::angleAxis(eulerAngles.z, glm::vec3(0.0, 0.0, 1.0));
                glm::vec3 axisZ = glm::normalize(parentRotate * glm::vec3(0.0, 0.0, 1.0));
                glm::vec3 axisY = glm::normalize(parentRotate * rotateZ * glm::vec3(0.0, 1.0, 0.0));
                glm::quat rotateY = glm::angleAxis(eulerAngles.y, glm::vec3(0.0, 1.0, 0.0));
                glm::vec3 axisX = glm::normalize(parentRotate * rotateZ * rotateY * glm::vec3(1.0, 0.0, 0.0));
#endif
                
                //计算当前节点到末端节点连线的向量
                glm::vec3 ri = ik.EndEffectorPosition() - ik.JointGlobalPosition[i];
                
                glm::vec3 dx = (glm::cross(axisX, ri));
                glm::vec3 dy = (glm::cross(axisY, ri));
                glm::vec3 dz = (glm::cross(axisZ, ri));
                
                jacobianTransposeMatrix(i * 3, 0) = dx.x;
                jacobianTransposeMatrix(i * 3, 1) = dx.y;
                jacobianTransposeMatrix(i * 3, 2) = dx.z;
                
                jacobianTransposeMatrix(i * 3 + 1, 0) = dy.x;
                jacobianTransposeMatrix(i * 3 + 1, 1) = dy.y;
                jacobianTransposeMatrix(i * 3 + 1, 2) = dy.z;
                
                jacobianTransposeMatrix(i * 3 + 2, 0) = dz.x;
                jacobianTransposeMatrix(i * 3 + 2, 1) = dz.y;
                jacobianTransposeMatrix(i * 3 + 2, 2) = dz.z;
            }
            
            // step3 计算角度增量和步长
            deltaTheta = jacobianTransposeMatrix * deltaX;
            
            Eigen::MatrixXf jacobianMatrix = jacobianTransposeMatrix.transpose();
            Eigen::Vector3f vec1 = jacobianMatrix * jacobianTransposeMatrix * deltaX;
            glm::vec3 vecTemp(vec1.x(), vec1.y(), vec1.z());
            
            // 计算stepsize
            float beta = glm::dot(delta, vecTemp) / glm::dot(vecTemp, vecTemp);
            
            //step4 应用上当前的角度增量
            for (int i = 0; i < ik.NumJoints(); i ++)
            {
                glm::vec3 eularAngle;
                eularAngle.x = deltaTheta(i * 3, 0);
                eularAngle.y = deltaTheta(i * 3 + 1, 0);
                eularAngle.z = deltaTheta(i * 3 + 2, 0);
                
                glm::vec3 currentAngle = glm::eulerAngles(ik.JointLocalRotation[i]);
                currentAngle -= eularAngle * beta;
                
                ik.JointLocalRotation[i] = glm::quat(currentAngle);
            }
            
            ForwardKinematics(ik, 0);
        }
        
//        for (int i = 0; i < ik.NumJoints(); i ++)
//        {
//            printf("%d ik point x = %f, y = %f, z = %f\n", i,
//                   ik.JointGlobalPosition[i].x, ik.JointGlobalPosition[i].y, ik.JointGlobalPosition[i].z);
//        }
//
//        printf("target point x = %f, y = %f, z = %f\n",
//               EndPosition.x, EndPosition.y, EndPosition.z);
        
        printf("jacobian 迭代次数 = %d\n", iteration);
    }

    void InverseKinematicsJacobianInverse(IKSystem & ik, const glm::vec3 & EndPosition, int maxIteration, float eps)
    {
        ForwardKinematics(ik, 0);
        
        int size = ik.NumJoints();
        if (size == 0)
        {
            return;
        }
        
        int iteration = 0;
        // These functions will be useful: glm::normalize, glm::rotation, glm::quat * glm::quat
        for (; iteration < maxIteration && glm::l2Norm(ik.EndEffectorPosition() - EndPosition) > eps; iteration++)
        {
            //step 1 计算出delta x，即当前末端点减去目标点的向量
            glm::vec3 delta = ik.EndEffectorPosition() - EndPosition;
            
            Eigen::MatrixXf deltaX(3, 1);
            deltaX(0, 0) = delta.x;
            deltaX(1, 0) = delta.y;
            deltaX(2, 0) = delta.z;
            
            //step 2: 计算雅可比矩阵
            
            Eigen::MatrixXf jacobianMatrix(3, 3 * ik.NumJoints());
            Eigen::MatrixXf deltaTheta(3 * ik.NumJoints(), 1);
            
            for (int i = 0; i < ik.NumJoints(); i ++)
            {
                // 1 获得父亲关节的全局旋转
                int parentIndex = i - 1;
                if (0 == i)
                {
                    parentIndex = 0;
                }
                const glm::quat& parentRotate = ik.JointGlobalRotation[parentIndex];
                
                // 2 获得当前节点的局部旋转，并转换为XYZ的转角顺序的欧拉角
                const glm::quat& localRatate = ik.JointLocalRotation[i];
                glm::vec3 eulerAngles = glm::eulerAngles(localRatate);
                
                // 3 计算当前节点在全局空间中的朝向
                glm::quat rotateX = glm::angleAxis(eulerAngles.x, glm::vec3(1.0, 0.0, 0.0));
                glm::vec3 axisX = glm::normalize(parentRotate * glm::vec3(1.0, 0.0, 0.0));
                glm::vec3 axisY = glm::normalize(parentRotate * rotateX * glm::vec3(0.0, 1.0, 0.0));
                
                glm::quat rotateY = glm::angleAxis(eulerAngles.y, glm::vec3(0.0, 1.0, 0.0));
                glm::vec3 axisZ = glm::normalize(parentRotate * rotateX * rotateY * glm::vec3(0.0, 0.0, 1.0));
                
                //计算当前节点到末端节点连线的向量
                glm::vec3 ri = ik.EndEffectorPosition() - ik.JointGlobalPosition[i];
                
                glm::vec3 dx = (glm::cross(axisX, ri));
                glm::vec3 dy = (glm::cross(axisY, ri));
                glm::vec3 dz = (glm::cross(axisZ, ri));
                
                jacobianMatrix(0, i * 3) = dx.x;
                jacobianMatrix(1, i * 3) = dx.y;
                jacobianMatrix(2, i * 3) = dx.z;
                
                jacobianMatrix(0, i * 3 + 1) = dy.x;
                jacobianMatrix(1, i * 3 + 1) = dy.y;
                jacobianMatrix(2, i * 3 + 1) = dy.z;
                
                jacobianMatrix(0, i * 3 + 2) = dz.x;
                jacobianMatrix(1, i * 3 + 2) = dz.y;
                jacobianMatrix(2, i * 3 + 2) = dz.z;
            }
            
            // step3 计算角度增量和步长
            //求雅可比矩阵的伪逆
            Eigen::MatrixXf jacobianTransposeMatrix = jacobianMatrix.transpose();
            Eigen::MatrixXf pseudoInverse = jacobianMatrix.completeOrthogonalDecomposition().pseudoInverse();;
            deltaTheta = pseudoInverse * deltaX;
            
//            Eigen::IOFormat CleanFmt(3, 0, ", ", "\n", "[", "]");
//            
//            std::cout << (jacobianMatrix * jacobianTransposeMatrix).inverse().format(CleanFmt) << std::endl;
//            std::cout << "下一次迭代" << std::endl;
            
            
//            Eigen::Vector3f vec1 = jacobianMatrix * jacobianTransposeMatrix * deltaX;
//            glm::vec3 vecTemp(vec1.x(), vec1.y(), vec1.z());
//            
//            // 计算stepsize
//            float beta = glm::dot(delta, vecTemp) / glm::dot(vecTemp, vecTemp);
            
            float beta = 1.0;
            
            //step4 应用上当前的角度增量
            for (int i = 0; i < ik.NumJoints(); i ++)
            {
                glm::vec3 eularAngle;
                eularAngle.x = deltaTheta(i * 3, 0);
                eularAngle.y = deltaTheta(i * 3 + 1, 0);
                eularAngle.z = deltaTheta(i * 3 + 2, 0);
                
                glm::vec3 currentAngle = glm::eulerAngles(ik.JointLocalRotation[i]);
                currentAngle -= eularAngle * beta;
                
                ik.JointLocalRotation[i] = glm::quat(currentAngle);
            }
            
            ForwardKinematics(ik, 0);
        }
        
//        for (int i = 0; i < ik.NumJoints(); i ++)
//        {
//            printf("%d ik point x = %f, y = %f, z = %f\n", i,
//                   ik.JointGlobalPosition[i].x, ik.JointGlobalPosition[i].y, ik.JointGlobalPosition[i].z);
//        }
//
//        printf("target point x = %f, y = %f, z = %f\n",
//               EndPosition.x, EndPosition.y, EndPosition.z);
        
        printf("jacobian inverse 迭代次数 = %d\n", iteration);
    }

    IKSystem::Vec3ArrPtr IKSystem::BuildCustomTargetPosition() {
        // get function from https://www.wolframalpha.com/input/?i=Albert+Einstein+curve
        int nums = 5000;
        using Vec3Arr = std::vector<glm::vec3>;
        std::shared_ptr<Vec3Arr> custom(new Vec3Arr(nums));
        int index = 0;
        for (int i = 0; i < nums; i++) {
            float x_val = 1.5e-3f * custom_x(92 * glm::pi<float>() * i / nums);
            float y_val = 1.5e-3f * custom_y(92 * glm::pi<float>() * i / nums);
            if (std::abs(x_val) < 1e-3 || std::abs(y_val) < 1e-3) continue;
            (*custom)[index++] = glm::vec3(1.6f - x_val, 0.0f, y_val - 0.2f);
        }
        custom->resize(index);
        return custom;
    }

    static Eigen::VectorXf glm2eigen(std::vector<glm::vec3> const & glm_v) {
        Eigen::VectorXf v = Eigen::Map<Eigen::VectorXf const, Eigen::Aligned>(reinterpret_cast<float const *>(glm_v.data()), static_cast<int>(glm_v.size() * 3));
        return v;
    }

    static std::vector<glm::vec3> eigen2glm(Eigen::VectorXf const & eigen_v) {
        return std::vector<glm::vec3>(
            reinterpret_cast<glm::vec3 const *>(eigen_v.data()),
            reinterpret_cast<glm::vec3 const *>(eigen_v.data() + eigen_v.size())
        );
    }

    static Eigen::SparseMatrix<float> CreateEigenSparseMatrix(std::size_t n, std::vector<Eigen::Triplet<float>> const & triplets) {
        Eigen::SparseMatrix<float> matLinearized(n, n);
        matLinearized.setFromTriplets(triplets.begin(), triplets.end());
        return matLinearized;
    }

    // solve Ax = b and return x
    static Eigen::VectorXf ComputeSimplicialLLT(
        Eigen::SparseMatrix<float> const & A,
        Eigen::VectorXf const & b) {
        auto solver = Eigen::SimplicialLLT<Eigen::SparseMatrix<float>>(A);
        return solver.solve(b);
    }

    void AdvanceMassSpringSystem(MassSpringSystem & system, float const dt) {
        // your code here: rewrite following code
        int const steps = 1000;
        float const ddt = dt / steps; 
        for (std::size_t s = 0; s < steps; s++) {
            std::vector<glm::vec3> forces(system.Positions.size(), glm::vec3(0));
            for (auto const spring : system.Springs) {
                auto const p0 = spring.AdjIdx.first;
                auto const p1 = spring.AdjIdx.second;
                glm::vec3 const x01 = system.Positions[p1] - system.Positions[p0];
                glm::vec3 const v01 = system.Velocities[p1] - system.Velocities[p0];
                glm::vec3 const e01 = glm::normalize(x01);
                glm::vec3 f = (system.Stiffness * (glm::length(x01) - spring.RestLength) + system.Damping * glm::dot(v01, e01)) * e01;
                forces[p0] += f;
                forces[p1] -= f;
            }
            for (std::size_t i = 0; i < system.Positions.size(); i++) {
                if (system.Fixed[i]) continue;
                system.Velocities[i] += (glm::vec3(0, -system.Gravity, 0) + forces[i] / system.Mass) * ddt;
                system.Positions[i] += system.Velocities[i] * ddt;
            }
        }
    }
}
