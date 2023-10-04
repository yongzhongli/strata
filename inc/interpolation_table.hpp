// Author: Damian Marek

// Copyright 2021 Shashwat Sharma and Piero Triverio

// This file is part of Strata.

// Strata is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// Strata is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with Strata.  If not, see <https://www.gnu.org/licenses/>.

/************************ interpolation_table.hpp ************************

 * Interpolation table for 3D datasets using piecewise polynomials.
 *
 * Author: Damian Marek
 * Created on: April 09, 2022

 ***************************************************************/

#ifndef STRATA_INTERPOLATION_TABLE_HPP
#define STRATA_INTERPOLATION_TABLE_HPP

#include <cassert>
#include <functional>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <complex>
#include <vector>
#include <algorithm>
#include <iostream>

using Int3D = std::array<int, 3>;
using Double3D = std::array<double, 3>;

inline Int3D operator-(const Int3D &lhs, const Int3D &rhs)
{
    Int3D diff;
    for (int i = 0; i < 3; i++)
        diff[i] = lhs[i] - rhs[i];
    return diff;
}

inline Int3D operator+(const Int3D &lhs, const Int3D &rhs)
{
    Int3D sum;
    for (int i = 0; i < 3; i++)
        sum[i] = lhs[i] + rhs[i];
    return sum;
}

inline Int3D operator+(const Int3D &lhs, const int rhs)
{
    Int3D sum;
    for (int i = 0; i < 3; i++)
        sum[i] = lhs[i] + rhs;
    return sum;
}
template <size_t N>
inline std::array<std::complex<double>, N> operator*(const std::array<std::complex<double>, N> &lhs,
                                                     const double rhs)
{
    std::array<std::complex<double>, N> result;
    for (int i = 0; i < N; i++)
        result[i] = lhs[i] * rhs;
    return result;
}
template <size_t N>
inline std::array<std::complex<double>, N> &
operator+=(std::array<std::complex<double>, N> &lhs, const std::array<std::complex<double>, N> &rhs)
{
    for (int i = 0; i < N; i++)
        lhs[i] += rhs[i];
    return lhs;
}

template <typename Iterator, typename T>
ptrdiff_t find_closest_point(T item, Iterator begin, Iterator end)
{
    assert(std::is_sorted(begin, end));
    assert(*begin <= item);
    assert(item <= *(end - 1));
    // check edge cases
    if (*begin == item)
        return 0;
    else if (item == *(end - 1))
        return end - 1 - begin;
    auto lb = std::lower_bound(begin, end, item);
    assert(lb != end);
    if (*lb == item)
    {
        // test point is coincident with a data point
        return lb - begin;
    }
    else // test < *lb
    {
        T other_b = *(lb - 1);
        return std::abs(other_b - item) < std::abs(*lb - item) ? lb - begin - 1 : lb - begin;
    }
}
/// Determine the interval where an item belongs.
/// Each interval is defined by a range [a, b) where "b" is not included
/// A collection of ranges are defined as {a, b, c, ...}, where the first interval is [a, b), and
/// the second is defined as [b, c), etc. The final interval will include its end point in the
/// interval.
template <typename Iterator, typename T>
ptrdiff_t find_interval(T item, Iterator begin, Iterator end)
{
    assert(std::is_sorted(begin, end));
    assert(*begin <= item);
    assert(item <= *(end - 1));
    // check edge cases
    if (*begin == item)
        return 0;
    else if (item == *(end - 1))
        return end - 1 - begin;

    auto ub = std::upper_bound(begin, end, item);
    assert(ub != end);
    ptrdiff_t interval = (ub - begin - 1);
    return interval;
}

inline std::vector<double> get_midpoints(const std::vector<double> &coords)
{
    std::vector<double> midpoints(coords.size() - 1);
    for (int i = 0; i < coords.size() - 1; i++)
    {
        midpoints[i] = 0.5 * (coords[i] + coords[i + 1]);
    }
    return midpoints;
}

/// Piecewise polynomial interpolation of 3D data
template <int N> class InterpolationTable
{
public:
    InterpolationTable() = default;
    /// Sets the grid and the associated data. Data is stored row-major, so the "z" coordinate
    /// changes the quickest.
    InterpolationTable(std::vector<double> _x, std::vector<double> _y, std::vector<double> _z,
                       const std::vector<std::array<std::complex<double>, N>> &_data)
            : xgrid(_x), ygrid(_y), zgrid(_z), data(_data)
    {
    }

    using InterFunction =
            std::function<std::array<std::complex<double>, N>(double x, double y, double z)>;
    /// Sets the grid and also sets data by a callable function oject
    InterpolationTable(std::vector<double> _x, std::vector<double> _y, std::vector<double> _z,
                       const InterFunction &function)
            : xgrid(_x), ygrid(_y), zgrid(_z)
    {
        // Set data using provided function
        data.resize(grid_size());
        int index = 0;
        for (int i = 0; i < xgrid.size(); i++)
        {
            double x = xgrid[i];
            for (int j = 0; j < ygrid.size(); j++)
            {
                double y = ygrid[j];
                for (int k = 0; k < zgrid.size(); k++)
                {
                    double z = zgrid[k];
                    data[index++] = function(x, y, z);
                }
            }
        }
    }

    /// Sets the grid and uses an existing table to set available values
    InterpolationTable(const InterpolationTable &existing_tbl, std::vector<double> _x,
                       std::vector<double> _y, std::vector<double> _z,
                       const InterFunction &function)
            : xgrid(_x), ygrid(_y), zgrid(_z)
    {
        assert(existing_tbl.xgrid.size() <= xgrid.size());
        assert(existing_tbl.ygrid.size() <= ygrid.size());
        assert(existing_tbl.zgrid.size() <= zgrid.size());
        this->stencil_size = existing_tbl.stencil_size;
        std::array<std::vector<double>, 3> existing_coords, new_coords;
        existing_coords[0] = (existing_tbl.xgrid);
        existing_coords[1] = (existing_tbl.ygrid);
        existing_coords[2] = (existing_tbl.zgrid);
        new_coords[0] = (xgrid);
        new_coords[1] = (ygrid);
        new_coords[2] = (zgrid);

        std::array<std::vector<double>, 3> common_coords, unique_coords;
        std::array<std::vector<std::pair<int, int>>, 3> common_coords_convert;
        std::array<std::vector<std::pair<double, int>>, 3> unique_coords_indices;
        std::array<std::unordered_map<int, int>, 3> map_to_old_data;
        // Find points in common and track the indices
        for (int i = 0; i < 3; i++)
        {
            common_coords[i].resize(existing_coords[i].size());
            map_to_old_data[i].reserve(existing_coords[i].size());
            auto end = std::set_intersection(existing_coords[i].cbegin(), existing_coords[i].cend(),
                                             new_coords[i].cbegin(), new_coords[i].cend(),
                                             common_coords[i].begin());
            common_coords[i].resize(end - common_coords[i].begin());
            common_coords_convert[i].resize(common_coords[i].size());
            for (int j = 0; j < common_coords[i].size(); j++)
            {
                auto it_new = std::lower_bound(new_coords[i].cbegin(), new_coords[i].cend(),
                                               common_coords[i][j]);
                assert(*it_new == common_coords[i][j]);
                auto it_existing = std::lower_bound(existing_coords[i].cbegin(),
                                                    existing_coords[i].cend(), common_coords[i][j]);
                assert(*it_existing == common_coords[i][j]);
                common_coords_convert[i][j] = {it_new - new_coords[i].cbegin(),
                                               it_existing - existing_coords[i].cbegin()};
                map_to_old_data[i].insert(common_coords_convert[i][j]);
            }
        }

        // Find points unique to this new grid
        for (int i = 0; i < 3; i++)
        {
            unique_coords[i].resize(new_coords[i].size() - common_coords[i].size());
            std::set_difference(new_coords[i].cbegin(), new_coords[i].cend(),
                                existing_coords[i].cbegin(), existing_coords[i].cend(),
                                unique_coords[i].begin());
            unique_coords_indices[i].resize(unique_coords[i].size());
            for (int j = 0; j < unique_coords[i].size(); j++)
            {
                auto it = std::lower_bound(new_coords[i].cbegin(), new_coords[i].cend(),
                                           unique_coords[i][j]);
                assert(*it == unique_coords[i][j]);
                unique_coords_indices[i][j] = {*it, it - new_coords[i].cbegin()};
            }
        }

        size_t num_coords = (common_coords[0].size() + unique_coords[0].size()) *
                            (common_coords[1].size() + unique_coords[1].size()) *
                            (common_coords[2].size() + unique_coords[2].size());

        data.resize(xgrid.size() * ygrid.size() * zgrid.size());
        assert(num_coords == data.size());

        size_t copied = 0, computed = 0;
        size_t index = 0;
        for (int i = 0; i < xgrid.size(); i++)
        {
            double x = xgrid[i];
            int x2 = -1;
            auto search = map_to_old_data[0].find(i);
            if (search != map_to_old_data[0].end())
                x2 = search->second;
            for (int j = 0; j < ygrid.size(); j++)
            {
                double y = ygrid[j];
                int y2 = -1;
                search = map_to_old_data[1].find(j);
                if (search != map_to_old_data[1].end())
                    y2 = search->second;
                for (int k = 0; k < zgrid.size(); k++)
                {
                    double z = zgrid[k];
                    int z2 = -1;
                    search = map_to_old_data[2].find(k);
                    if (search != map_to_old_data[2].end())
                        z2 = search->second;
                    if (x2 == -1 || y2 == -1 || z2 == -1)
                    {
                        computed++;
                        data[index++] = function(x, y, z);
                    }
                    else
                    {
                        copied++;
                        data[index++] = existing_tbl.data[existing_tbl.get_index({x2, y2, z2})];
                    }
                }
            }
        }
        std::cout << "Num Computed: " << computed << "\t Num Copied: " << copied << std::endl;
    }

    /// Compute the function value by interpolation at the specified point
    std::array<std::complex<double>, N> compute_at(double x, double y, double z) const
    {
        Int3D stencil_ijk = locate_stencil(x, y, z);
        std::array<std::complex<double>, N> result;
        result.fill(0.0);
        int LX = stencil_size[0];
        int LY = stencil_size[1];
        int LZ = stencil_size[2];

        Int3D Lijk = {0, 0, 0};
        Double3D position = {x, y, z};
        for (int i = 0; i < LX; i++)
        {
            Lijk[0] = i;
            for (int j = 0; j < LY; j++)
            {
                Lijk[1] = j;
                for (int k = 0; k < LZ; k++)
                {
                    Lijk[2] = k;
                    Int3D data_ijk = Lijk + stencil_ijk;
                    auto y_ijk = data[get_index(data_ijk)];
                    double l_ijk = compute_lagrange_polynomial(stencil_ijk, Lijk, position);
                    result += y_ijk * l_ijk;
                }
            }
        }
        return result;
    }

    /// Total number of grid points
    int grid_size() const { return xgrid.size() * ygrid.size() * zgrid.size(); }

    void get_grids(std::vector<double> &_x, std::vector<double> &_y, std::vector<double> &_z) const
    {
        _x = xgrid;
        _y = ygrid;
        _z = zgrid;
    }

    /// The number of stencils is directly related to the number of points.
    Int3D num_stencils() const
    {
        Int3D num_points = {xgrid.size(), ygrid.size(), zgrid.size()};
        return num_points - stencil_size + 1;
    }

    /// Helpful functions for converting from (i,j,k) indices to
    /// actual index in row major array, will check range in debug.
    int get_index(const Int3D &ijk) const
    {
        assert(ijk[0] >= 0 && ijk[1] >= 0 && ijk[2] >= 0);
        assert(ijk[0] < xgrid.size() && ijk[1] < ygrid.size() && ijk[2] < zgrid.size());
        const int ny = ygrid.size();
        const int nz = zgrid.size();
        return ijk[2] + nz * (ijk[1] + ny * ijk[0]);
    }

    /// Reverse of above, given an index returns the (i,j,k) indices, will check range in debug.
    Int3D get_ijk(const int index) const
    {
        assert(index >= 0);
        assert(index < grid_size());
        const int ny = ygrid.size();
        const int nz = zgrid.size();
        Int3D ijk;

        // Row Major - slightly optimized
        ijk[1] = index / nz;
        ijk[2] = index % nz;
        ijk[0] = ijk[1] / ny;
        ijk[1] = ijk[1] % ny;
        return ijk;
    }
    /// Number of interpolation points along each dimension. Equivalent to p+1, where p is the
    /// order of polynomial used for interpolation.
    Int3D stencil_size = {3, 3, 3};

private:
    Int3D locate_stencil(double x, double y, double z) const
    {
        Int3D stencil_ijk;
        stencil_ijk[0] = locate_stencil_1d(x, xgrid, stencil_size[0]);
        stencil_ijk[1] = locate_stencil_1d(y, ygrid, stencil_size[1]);
        stencil_ijk[2] = locate_stencil_1d(z, zgrid, stencil_size[2]);
        return stencil_ijk;
    }

    static int locate_stencil_1d(double coord, std::vector<double> grid, int stencil_size)
    {
        auto begin = grid.cbegin();
        auto end = grid.cend();
        int stencil_1d = -1;
        int max_stencil = grid.size() - stencil_size;
        // Even stencils are centered on intervals
        if (stencil_size % 2 == 0)
        {
            stencil_1d = find_interval(coord, begin, end);
            // stencil indices are defined by the first point that belongs to them
            // So the width of the stencil is subtracted from the interval location
            stencil_1d -= stencil_size / 2 - 1;
        }
        else // Odd stencils are centered on grid points
        {
            stencil_1d = find_closest_point(coord, begin, end);
            stencil_1d -= (stencil_size - 1) / 2;
        }
        // Take care of stencil indices that were computed past the boundary
        stencil_1d = std::max(0, stencil_1d);
        stencil_1d = std::min(max_stencil, stencil_1d);

        return stencil_1d;
    }
    /// Computes one of the lagrange polynomials, given by Lijk, at a position within a given
    /// stencil.
    double compute_lagrange_polynomial(Int3D stencil_ijk, Int3D Lijk, Double3D pos) const
    {
        double product = 1.0;
        int LX = stencil_size[0];
        int LY = stencil_size[1];
        int LZ = stencil_size[2];
        // Compute product along x dimension
        for (int i = 0; i < LX; ++i)
        {
            if (i != Lijk[0])
            {
                const double x_m = xgrid[i + stencil_ijk[0]];
                const double x_idx = xgrid[Lijk[0] + stencil_ijk[0]];

                product *= (pos[0] - x_m) / (x_idx - x_m);
            }
        }
        // Compute product along y dimension
        for (int j = 0; j < LY; ++j)
        {
            if (j != Lijk[1])
            {
                const double y_m = ygrid[j + stencil_ijk[1]];
                const double y_idx = ygrid[Lijk[1] + stencil_ijk[1]];

                product *= (pos[1] - y_m) / (y_idx - y_m);
            }
        }
        // Compute product along z dimension
        for (int k = 0; k < LZ; ++k)
        {
            if (k != Lijk[2])
            {
                const double z_m = zgrid[k + stencil_ijk[2]];
                const double z_idx = zgrid[Lijk[2] + stencil_ijk[2]];

                product *= (pos[2] - z_m) / (z_idx - z_m);
            }
        }
        return product;
    }
    /// The 3D structured grid.
    std::vector<double> xgrid{}, ygrid{}, zgrid{};
    /// Stores many components or fields.
    std::vector<std::array<std::complex<double>, N>> data;
};

template <int N>
double check_max_error(const std::array<std::complex<double>, N> &ref,
                       const std::array<std::complex<double>, N> &trial)
{
    double max_error = std::numeric_limits<double>::min();
    for (int i = 0; i < N; i++)
    {
        double error = std::abs(ref[i] - trial[i]) / std::abs(ref[i]);
        max_error = std::max(error, max_error);
    }
    return max_error;
}

template <int N>
bool check_interpolation_and_update_grid(
        InterpolationTable<N> &trial, InterpolationTable<N> &new_table,
        const std::function<std::array<std::complex<double>, N>(double x, double y, double z)>
        &function,
        double &max_err_observed, double mid_rtol = 1e-4)
{
    std::vector<double> xgrid, ygrid, zgrid;
    trial.get_grids(xgrid, ygrid, zgrid);
    std::set<double> xgrid_add, ygrid_add, zgrid_add;
    max_err_observed = std::numeric_limits<double>::min();
    bool table_changed = false;
    for (int i = 0; i < xgrid.size() - 1; i++)
    {
        double x = 0.5 * (xgrid[i] + xgrid[i + 1]);
        for (int j = 0; j < ygrid.size() - 1; j++)
        {
            double y = 0.5 * (ygrid[j] + ygrid[j + 1]);
            for (int k = 0; k < zgrid.size() - 1; k++)
            {
                double z = 0.5 * (zgrid[k] + zgrid[k + 1]);
                auto reference = function(x, y, z);
                auto result = trial.compute_at(x, y, z);
                auto max_error = check_max_error<N>(reference, result);
                max_err_observed = std::max(max_error, max_err_observed);
                if (max_error > mid_rtol)
                {
                    table_changed = true;
                    xgrid_add.insert(x);
                    ygrid_add.insert(y);
                    zgrid_add.insert(z);
                }
            }
        }
    }
    if (table_changed)
    {
        xgrid.insert(xgrid.end(), xgrid_add.cbegin(), xgrid_add.cend());
        ygrid.insert(ygrid.end(), ygrid_add.cbegin(), ygrid_add.cend());
        zgrid.insert(zgrid.end(), zgrid_add.cbegin(), zgrid_add.cend());
        std::sort(xgrid.begin(), xgrid.end());
        std::sort(ygrid.begin(), ygrid.end());
        std::sort(zgrid.begin(), zgrid.end());
        InterpolationTable<N> tmp_table(xgrid, ygrid, zgrid, function);
        tmp_table.stencil_size = new_table.stencil_size;
        std::swap(new_table, tmp_table);
    }
    return table_changed;
}

template <int N>
bool check_interpolation_one_dimension_and_update_grid(
        int dim, InterpolationTable<N> &trial, InterpolationTable<N> &new_table,
        const std::function<std::array<std::complex<double>, N>(double, double, double)> &function,
        double &max_err_observed, double mid_rtol = 1e-4)
{
    assert(dim >= 0 && dim < 3);
    std::array<std::vector<double>, 3> coords;
    trial.get_grids(coords[0], coords[1], coords[2]);

    max_err_observed = std::numeric_limits<double>::min();

    std::array<std::vector<double> *, 3> trial_coords = {&coords[0], &coords[1], &coords[2]};

    auto coord_midpoints = get_midpoints(coords[dim]);
    trial_coords[dim] = &coord_midpoints;
    std::unordered_set<double> coord_add(coord_midpoints.size());

    bool table_changed = false;
    std::array<double, 3> coord{};
    for (auto &x : *trial_coords[0])
    {
        coord[0] = x;
        for (auto &y : *trial_coords[1])
        {
            coord[1] = y;
            for (auto &z : *trial_coords[2])
            {
                coord[2] = z;
                if (coord_add.find(coord[dim]) != coord_add.cend())
                    continue;
                auto reference = function(x, y, z);
                auto result = trial.compute_at(x, y, z);
                auto max_error = check_max_error<N>(reference, result);
                max_err_observed = std::max(max_error, max_err_observed);
                if (max_error > mid_rtol)
                {
                    table_changed = true;
                    coord_add.insert(coord[dim]);
                }
            }
        }
    }
    if (table_changed)
    {
        coords[dim].insert(coords[dim].end(), coord_add.cbegin(), coord_add.cend());
        std::sort(coords[dim].begin(), coords[dim].end());
        InterpolationTable<N> tmp_table(trial, coords[0], coords[1], coords[2], function);
        // InterpolationTable<N> tmp_table(coords[0], coords[1], coords[2], function);
        tmp_table.stencil_size = trial.stencil_size;
        std::swap(new_table, tmp_table);
    }
    return table_changed;
}

template <int N>
bool check_interpolation_each_dimension_and_update_grid(
        InterpolationTable<N> &trial, InterpolationTable<N> &new_table,
        const std::function<std::array<std::complex<double>, N>(double, double, double)> &function,
        double &max_err_observed, double mid_rtol = 1e-4)
{
    double max_err2, max_err1, max_err0;
    bool changed2 = check_interpolation_one_dimension_and_update_grid<N>(
            2, trial, new_table, function, max_err2, mid_rtol);

    bool changed1 = check_interpolation_one_dimension_and_update_grid<N>(
            1, new_table, new_table, function, max_err1, mid_rtol);
    bool changed0 = check_interpolation_one_dimension_and_update_grid<N>(
            0, new_table, new_table, function, max_err0, mid_rtol);
    max_err_observed = std::max({max_err2, max_err1, max_err0});
    return changed2 || changed1 || changed0;
}

#endif // STRATA_INTERPOLATION_TABLE_HPP
