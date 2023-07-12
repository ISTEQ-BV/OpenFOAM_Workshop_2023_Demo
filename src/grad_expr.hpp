/*
 *
 * Copyright: 2023 Ilya Popov <ilya.popov@isteq.nl>, ISTEQ BV
 *
 * SPDX-License-Identifier: GPL3.0-or-later
 *
 */

#pragma once

#include "expressions.hpp"

#include "leastSquaresVectors.H"

namespace Foam {
namespace fve {

template<typename CellExpr>
struct grad_expr {
    using value_type = typename Foam::outerProduct<Foam::vector, typename CellExpr::value_type>::type;
    static_assert(CellExpr::location == loc::cell, "Argument to gradient must be cell expression");
    static constexpr loc location = loc::cell;
    static constexpr bool has_surface_integrate = false;

    CellExpr nested;
    const Foam::leastSquaresVectors& lsv;
    const Foam::labelUList& owner;
    const Foam::labelUList& neighbour;
    const Foam::cellList& cells;
    const Foam::surfaceVectorField& pVectors;
    const Foam::surfaceVectorField& nVectors;

    grad_expr(const CellExpr& arg)
        : nested(arg)
        , lsv(Foam::leastSquaresVectors::New(arg.mesh()))
        , owner(arg.mesh().owner())
        , neighbour(arg.mesh().neighbour())
        , cells(arg.mesh().cells())
        , pVectors(lsv.pVectors())
        , nVectors(lsv.nVectors())
    {}

    value_type operator [](Foam::label celli) const {
        const Foam::cell &cell = cells[celli];

        value_type grad{};
        auto val = nested[celli];

        for (Foam::label facei: cell) {
            if (facei < owner.size()) {
                auto own = owner[facei];
                auto nei = neighbour[facei];
                if (own == celli) {
                    grad += pVectors[facei] * (nested[nei] - val);
                }
                else {
                    grad += nVectors[facei] * (nested[own] - val);
                }
            }
            else {
                // TODO
            }
        }

        return grad;
    }

    value_type on_boundary(Foam::label patchi, Foam::label facei) const {
        // TODO
        return value_type::zero;
    }

    const fvMesh& mesh() const {
        return nested.mesh();
    }

    Foam::dimensionSet dimensions() const {
        // This is a hack just to represent unknown dimensions
        return nested.dimensions() / Foam::dimLength;
    }
};

template<typename Expr>
struct is_expression<grad_expr<Expr>> : std::true_type {};

template <typename CellExpr, typename std::enable_if<is_expression<CellExpr>::value, int>::type = 0>
auto grad(const CellExpr& arg) -> grad_expr<CellExpr> {
    return {arg};
}


} // namespace fve
} // namespace Foam
