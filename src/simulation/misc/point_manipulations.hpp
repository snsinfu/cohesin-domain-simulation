#pragma once

#include <random>
#include <vector>

#include <md.hpp>

#include "../common_types.hpp"


std::vector<md::point> generate_random_walk(md::index n_steps, md::scalar step_size, random_engine& random);
std::vector<md::point> interpolate_points(md::array_view<md::point const> knots, md::index n_points);
void                   move_centroid(md::array_view<md::point> points, md::point target);
