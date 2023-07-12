/*
 *
 * Copyright: 2023 Ilya Popov <ilya.popov@isteq.nl>, ISTEQ BV
 *
 * SPDX-License-Identifier: GPL3.0-or-later
 *
 */

#pragma once

#include "expressions.hpp"

#include "microdomains.hpp"

#include "process_microdomains.hpp"

#include "leastSquaresVectors.H"

namespace Foam {
namespace fve {

template <typename Field, typename CellExpr>
struct grad_expr_2 {
    using value_type = typename Foam::outerProduct<Foam::vector, typename CellExpr::value_type>::type;
    static_assert(CellExpr::location == loc::cell, "Argument to gradient must be cell expression");
    static constexpr loc location = loc::cell;
    static constexpr bool has_surface_integrate = true;

    Field& field;
    CellExpr nested;
    const Foam::leastSquaresVectors& lsv;
    const Foam::labelUList& owner;
    const Foam::labelUList& neighbour;
    const Foam::cellList& cells;
    const Foam::surfaceVectorField& pVectors;
    const Foam::surfaceVectorField& nVectors;

    grad_expr_2(Field& f, const CellExpr& arg)
        : field(f)
        , nested(arg)
        , lsv(Foam::leastSquaresVectors::New(arg.mesh()))
        , owner(arg.mesh().owner())
        , neighbour(arg.mesh().neighbour())
        , cells(arg.mesh().cells())
        , pVectors(lsv.pVectors())
        , nVectors(lsv.nVectors())
    {}

    value_type operator[](Foam::label celli) const {
        return field[celli];
    }

    value_type on_boundary(Foam::label patchi, Foam::label facei) const {
        return field.boundaryField()[patchi][facei];
    }

    void process_face(Foam::label facei) const {
        auto own = owner[facei];
        auto nei = neighbour[facei];
        auto delta_v = nested[nei] - nested[own];
        field[own] += pVectors[facei] * delta_v;
        field[nei] -= nVectors[facei] * delta_v;
    }

    void process_microdomain(const microdomain& md) const {
        //std::cout << "Processing microdomain " << md.internal_faces << " + " << md.own_boundary_faces << std::endl;
        for (Foam::label facei: md.internal_faces) {
            process_face(facei);
        }

        for (Foam::label facei: md.own_boundary_faces) {
            process_face(facei);
        }
    }

    void process_boundary() const {
        // TODO:
    }

    const fvMesh& mesh() const {
        return nested.mesh();
    }

    Foam::dimensionSet dimensions() const {
        // check dimensions
        field.dimensions() = nested.dimensions() / Foam::dimLength;
        return nested.dimensions() / Foam::dimLength;
    }
};

template<typename Field, typename Expr>
struct is_expression<grad_expr_2<Field, Expr>> : std::true_type {};

template <typename Field, typename CellExpr, typename std::enable_if<is_expression<CellExpr>::value, int>::type = 0>
auto grad(Field& f, const CellExpr& arg) -> grad_expr_2<Field, CellExpr> {
    return {f, arg};
}

template <typename Field, typename CellExpr>
void process_microdomain(const grad_expr_2<Field, CellExpr>& expr, const microdomain& md)
{
    expr.process_microdomain(md);
}

} // namespace fve
} // namespace Foam
