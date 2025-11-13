#include <random>

#include <spline.hpp>

#include "point_manipulations.hpp"


std::vector<md::point>
generate_random_walk(md::index n_steps, md::scalar step_size, random_engine& random)
{
    std::vector<md::point> points;
    points.reserve(n_steps);

    auto direction = [](auto& random) {
        std::normal_distribution<md::scalar> normal;
        return md::vector{normal(random), normal(random), normal(random)}.normalize();
    };

    md::point pos;
    for (md::index i = 0; i < n_steps; i++) {
        points.push_back(pos);
        pos += step_size * direction(random);
    }

    return points;
}


std::vector<md::point>
interpolate_points(md::array_view<md::point const> knots, md::index n_points)
{
    std::vector<double> input_coords;
    for (md::index i = 0; i < knots.size(); i++) {
        input_coords.push_back(double(i) / double(knots.size() - 1));
    }

    std::vector<double> knots_x, knots_y, knots_z;
    for (auto const& knot : knots) {
        knots_x.push_back(knot.x);
        knots_y.push_back(knot.y);
        knots_z.push_back(knot.z);
    }

    cubic_spline spline_x(input_coords, knots_x, cubic_spline::natural);
    cubic_spline spline_y(input_coords, knots_y, cubic_spline::natural);
    cubic_spline spline_z(input_coords, knots_z, cubic_spline::natural);

    std::vector<md::point> points;
    for (md::index i = 0; i < n_points; i++) {
        auto const t = double(i) / double(n_points - 1);
        auto const point = md::point{spline_x(t), spline_y(t), spline_z(t)};
        points.push_back(point);
    }

    return points;
}


void
move_centroid(md::array_view<md::point> points, md::point target)
{
    md::vector offset;
    for (auto const& point : points) {
        offset += point - target;
    }
    offset /= md::scalar(points.size());

    for (auto& point : points) {
        point -= offset;
    }
}
