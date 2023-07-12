/*
 *
 * Copyright: 2023 Ilya Popov <ilya.popov@isteq.nl>, ISTEQ BV
 *
 * SPDX-License-Identifier: GPL3.0-or-later
 *
 */

#pragma once

#include "expressions.hpp"

#include "boost/mp11/algorithm.hpp"
#include "boost/mp11/tuple.hpp"

namespace Foam {
namespace fve {

template<typename Func, typename... Args>
struct map_expr {
    using value_type = decltype(std::declval<Func>()(std::declval<typename Args::value_type>()...));
    static constexpr bool all_args_vol = boost::mp11::mp_all_of<boost::mp11::mp_list<Args...>, is_vol_expr>::value;
    static constexpr bool all_args_surf = boost::mp11::mp_all_of<boost::mp11::mp_list<Args...>, is_surface_expr>::value;
    static_assert(all_args_vol || all_args_surf, "All arguments to map expr must be the same location, either all cell or all face.");
    static constexpr loc location = all_args_vol? loc::cell : loc::face;
    static constexpr bool has_surface_integrate = false; //FIXME

    Func func;
    std::tuple<Args&...> args;

    map_expr(Func&& func, Args&&... args)
        : func(std::forward<Func>(func))
        , args(std::forward_as_tuple(args...))
    {}

    value_type operator [](Foam::label i) const {
        return boost::mp11::tuple_apply(func, boost::mp11::tuple_transform(subscript{i}, args));
    }

    value_type on_boundary(Foam::label patchi, Foam::label facei) const {
        return boost::mp11::tuple_apply(func, boost::mp11::tuple_transform(boundary{patchi, facei}, args));
    }

    const fvMesh& mesh() const {
        return std::get<0>(args).mesh();
    }

    Foam::dimensionSet dimensions() const {
        // This is a hack just to represent unknown dimensions
        return Foam::dimless;
    }
};

template <typename Func, typename... Args>
auto map(Func&& func, Args&&... args) -> map_expr<Func, Args...> {
    return {std::forward<Func>(func), std::forward<Args>(args)...};
}

} // namespace fve
} // namespace Foam
