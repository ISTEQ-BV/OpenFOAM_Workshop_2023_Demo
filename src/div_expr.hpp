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

namespace Foam {
namespace fve {

enum class where {
    internal_and_boundary_faces,
    cells,
    internal_faces,
};

template<typename Field, typename FaceExpr>
struct surface_integrate_expr {
    using value_type = typename FaceExpr::value_type;
    static_assert(FaceExpr::location == loc::face, "Argument to gradient must be cell expression");
    static_assert(!FaceExpr::has_surface_integrate, "More than one level of surface integrate expressions is not handled");
    static constexpr loc location = loc::cell;
    static constexpr bool has_surface_integrate = true;

    Field& field;
    FaceExpr nested;
    const Foam::labelUList& owner;
    const Foam::labelUList& neighbour;
    const DimensionedField<scalar, volMesh>& V;

    surface_integrate_expr(const Field& field, const FaceExpr& arg)
        : field(field)
        , nested(arg)
        , owner(arg.mesh().owner())
        , neighbour(arg.mesh().neighbour())
        , V(nested.mesh().V())
    {
        const fvMesh& mesh = arg.mesh();
        for (label patchi = 0; patchi < mesh.boundary().size(); ++patchi) {
            const auto& patch = mesh.boundary()[patchi];
            const auto& pFaceCells = patch.faceCells();

            for (label facei = 0; facei < patch.size(); facei++) {
                label celli = pFaceCells[facei];
                field[celli] += nested.on_boundary(patchi, facei) / V[celli];
            }
        }
    }

    void process_face(label facei) {
        label own = owner[facei];
        label nei = neighbour[facei];
        auto face_val = nested[facei];

        field[own] += face_val / V[own];
        field[nei] -= face_val / V[nei];
        //TODO: try moving division by V into a separate loop
    }

    void process_microdomain(const microdomain& md) {
        for (auto facei: md.internal_faces) {
            process_face(facei);
        }
        for (auto facei: md.own_boundary_faces) {
            process_face(facei);
        }
    }

    const value_type operator[](Foam::label celli) const {
        return field[celli];
    }

    const value_type on_boundary(Foam::label patchi, Foam::label celli) const {
        return field.boundaryField()[patchi][celli];
    }

    const fvMesh& mesh() const {
        return nested.mesh();
    }

    Foam::dimensionSet dimensions() const {
        return nested.dimensions() / Foam::dimVolume;
    }
};

template<typename Field, typename Expr>
struct is_expression<surface_integrate_expr<Field, Expr>> : std::true_type {};

template <typename Field, typename FaceExpr, typename std::enable_if<is_expression<FaceExpr>::value, int>::type = 0>
auto surfaceIntegrate(Field& f, const FaceExpr& arg) -> surface_integrate_expr<Field, FaceExpr> {
    return {f, arg};
}

template <typename Field, typename Expr,
         typename std::enable_if<is_expression<Expr>::value && Expr::location == loc::face, int>::type = 0>
auto div(Field& f, const Expr& arg) -> surface_integrate_expr<Field, Expr> {
    return surfaceIntegrate(arg);
}

template <typename Field, typename Expr,
         typename std::enable_if<is_expression<Expr>::value && Expr::location == loc::cell, int>::type = 0>
auto div(Field& f, const Expr& arg) -> surface_integrate_expr<Field, dot_expr<field_expr<surfaceVectorField>, linear_interpolate_expr<Expr>>> {
    return surfaceIntegrate(fve::read(arg.mesh().Sf()) & fve::interpolate(arg));
}

template <typename Field, typename FaceExpr, typename std::enable_if<is_expression<FaceExpr>::value, int>::type = 0>
auto gauss_grad(Field& f, const FaceExpr& arg) -> surface_integrate_expr<Field, mul_expr<field_expr<surfaceVectorField>, linear_interpolate_expr<FaceExpr>>> {
    return surfaceIntegrate(fve::read(arg.mesh().Sf()) * fve::interpolate(arg));
}

} // namespace fve
} // namespace Foam
