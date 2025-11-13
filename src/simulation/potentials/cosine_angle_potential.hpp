#pragma once

#include <cmath>
#include <tuple>

#include <md/basic_types.hpp>


namespace md
{
    // cosine_angle_potential implements the three-body bending potential
    // function for constraining the bond angle θ around a preferred value:
    //
    //     u(r,s) = e (1 - cos(θ - θ0))
    //
    // where e is the strength of the potential, θ0 is the preferred angle,
    // and θ is the angle between r and s:
    //
    //     cos(θ) = dot(r,s) / (|r| |s|) .
    //
    struct cosine_angle_potential
    {
        md::scalar bending_energy = 0;
        md::scalar preferred_angle = 0;

        md::scalar evaluate_energy(md::vector rij, md::vector rjk) const
        {
            md::scalar const dij_sq = rij.squared_norm();
            md::scalar const djk_sq = rjk.squared_norm();
            md::scalar const dij_djk = std::sqrt(dij_sq * djk_sq);

            if (dij_djk == 0) [[unlikely]] {
                return 0;
            }

            md::scalar const cos_0 = std::cos(preferred_angle);
            md::scalar const sin_0 = std::sin(preferred_angle);

            md::scalar const dot_ijk = md::dot(rij, rjk);
            md::scalar const cross_ijk = md::cross(rij, rjk).norm();
            md::scalar const cos_delta = (dot_ijk * cos_0 + cross_ijk * sin_0) / dij_djk;

            return bending_energy - bending_energy * cos_delta;
        }

        std::tuple<md::vector, md::vector, md::vector> evaluate_force(
            md::vector rij, md::vector rjk
        ) const
        {
            md::scalar const dij_sq = rij.squared_norm();
            md::scalar const djk_sq = rjk.squared_norm();
            md::scalar const dij_djk = std::sqrt(dij_sq * djk_sq);

            if (dij_djk == 0) [[unlikely]] {
                return std::make_tuple(md::vector{}, md::vector{}, md::vector{});
            }

            md::scalar const dot_ijk = md::dot(rij, rjk);
            md::scalar const cross_ijk = md::cross(rij, rjk).norm();

            if (cross_ijk == 0) [[unlikely]] {
                return std::make_tuple(md::vector{}, md::vector{}, md::vector{});
            }

            md::scalar const cos_0 = std::cos(preferred_angle);
            md::scalar const sin_0 = std::sin(preferred_angle);
            md::scalar const omega_ijk = cos_0 - sin_0 * dot_ijk / cross_ijk;

            md::scalar const e_div_dd = bending_energy * omega_ijk / dij_djk;
            md::vector const fij = e_div_dd * (rjk - dot_ijk / dij_sq * rij);
            md::vector const fjk = e_div_dd * (rij - dot_ijk / djk_sq * rjk);

            return std::make_tuple(fij, fjk - fij, -fjk);
        }
    };
}
