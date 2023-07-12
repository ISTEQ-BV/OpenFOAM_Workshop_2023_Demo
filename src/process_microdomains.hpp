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

#include "boost/mp11/function.hpp"

namespace Foam {
namespace fve {

template <typename Type, template<class> class PatchField, typename Expression,
         typename std::enable_if<Expression::has_surface_integrate, int>::type = 0>
[[gnu::noinline]]
void operator<<=(Foam::GeometricField<Type, PatchField, Foam::volMesh>& f, Expression e)
{
    static_assert(Expression::location == loc::cell,
                  "Expression must have same location (cell or face) as the target field");
    const fvMesh& mesh = e.mesh();

    if (!e.dimensions().dimensionless()) {
        // Check dimensions
        f.dimensions() = e.dimensions();
    }

    const auto mds = microdomains::New(mesh);

    for (const auto&  md: mds.domains) {
        process_microdomain(e, md);

        for (auto celli: md.cells) {
            f[celli] = e[celli];
        }
    }

    Foam::label nPatches = mesh.boundary().size();
    for (Foam::label patchi = 0; patchi < nPatches; patchi++) {
        PatchField<Type>& patchField = f.boundaryFieldRef()[patchi];
        Foam::label nFaces = patchField.size();

        for (Foam::label facei = 0; facei < nFaces; facei++) {
            patchField[facei] = e.on_boundary(patchi, facei);
        }
    }
}

template <typename Type, template<class> class PatchField, typename Expression,
         typename std::enable_if<Expression::has_surface_integrate, int>::type = 0>
[[gnu::noinline]]
void operator<<=(Foam::GeometricField<Type, PatchField, Foam::surfaceMesh>& f, Expression e)
{
    static_assert(Expression::location == loc::face,
                  "Expression must have same location (cell or face) as the target field");
    const fvMesh& mesh = e.mesh();

    if (!e.dimensions().dimensionless()) {
        // Check dimensions
        f.dimensions() = e.dimensions();
    }

    const auto mds = microdomains::New(mesh);

    for (const auto&  md: mds.domains) {
        process_microdomain(e, md);

        // after we precomputed a domain, we can compute values in internal faces
        for (auto facei: md.internal_faces) {
            f[facei] = e[facei];
        }
    }

    // after we computed all domains, we can compute boundary faces between domains
    for (const auto&  md: mds.domains) {
        for (auto facei: md.own_boundary_faces) {
            f[facei] = e[facei];
        }
    }

    Foam::label nPatches = mesh.boundary().size();
    for (Foam::label patchi = 0; patchi < nPatches; patchi++) {
        PatchField<Type>& patchField = f.boundaryFieldRef()[patchi];
        Foam::label nFaces = patchField.size();

        for (Foam::label facei = 0; facei < nFaces; facei++) {
            patchField[facei] = e.on_boundary(patchi, facei);
        }
    }
}

template<typename, typename = void>
struct has_nested : std::false_type {};

template<typename T>
struct has_nested<T, typename boost::mp11::mp_void<decltype(std::declval<T>().nested)> > : std::true_type {};

template<typename, typename = void>
struct has_lhs_rhs : std::false_type {};

template<typename T>
struct has_lhs_rhs<T, typename boost::mp11::mp_void<decltype(std::declval<T>().lhs), decltype(std::declval<T>().rhs)> > : std::true_type {};

template <typename Expr, typename std::enable_if<has_nested<Expr>::value, int>::type = 0>
void process_microdomain(const Expr& e, const microdomain& md) {
    process_microdomain(e.nested, md);
}

template <typename Expr, typename std::enable_if<has_lhs_rhs<Expr>::value, int>::type = 0>
void process_microdomain(const Expr& e, const microdomain& md) {
    process_microdomain(e.lhs, md);
    process_microdomain(e.rhs, md);
}

template <typename Field>
void process_microdomain(const field_expr<Field>& e, const microdomain& md) {
    // do nothing
}

} // namespace fve
} // namespace Foam
