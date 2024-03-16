#pragma once

#include <utility>
#include <vector>

#include <glm/glm.hpp>

namespace VCX::Labs::Animation {

    //质点弹簧系统
    struct MassSpringSystem {
        struct Spring {
            std::pair<std::size_t, std::size_t> AdjIdx;        //弹簧边的索引
            float                               RestLength;   //弹簧原长度
        };

        std::vector<glm::vec3>   Positions;    //位置
        std::vector<glm::vec3>   Velocities;   //速度
        std::vector<int>         Fixed;       //固定节点
        float                    Mass      { 1 };

        std::vector<Spring>      Springs;           //所有的弹簧
        float                    Stiffness { 100 };   //弹性系数
        float                    Damping   { .2f };   //摩擦系数
        float                    Gravity   { .3f };   //重力系数

        void AddParticle(glm::vec3 const & position, glm::vec3 const & velocity = glm::vec3(0)) {
            Positions.push_back(position);
            Velocities.push_back(velocity);
            Fixed.push_back(false);
        }

        void AddSpring(std::size_t const adjIdx0, std::size_t const adjIdx1, float const restLength = -1) {
            Springs.push_back({
                .AdjIdx     { adjIdx0, adjIdx1 },
                .RestLength { restLength < 0 ? glm::length(Positions[adjIdx0] - Positions[adjIdx1]) : restLength } 
            });
        }
    };
}
