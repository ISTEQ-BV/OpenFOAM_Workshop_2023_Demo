/*
 *
 * Copyright: 2023 Ilya Popov <ilya.popov@isteq.nl>, ISTEQ BV
 *
 * SPDX-License-Identifier: GPL3.0-or-later
 *
 */

#pragma once

#include "fvMesh.H"
#include "GeometricField.H"
#include "surfaceMesh.H"
#include "volMesh.H"

namespace Foam {
namespace fve {

enum class loc {
    cell,
    face,
};

//struct expression {
//    using value_type = ...;
//    static constexpr loc location = ...;
//    static constexpr bool has_surface_integrate = ...;
//    value_type operator[](Foam::label facei) const;
//    value_type on_boundary(Foam::label patchi, Foam::label facei) const;
//    void precompute_microdomain()
//    const Foam::fvMesh& mesh() const;
//    Foam::dimensionSet dimensions() const;
//};

template <typename T>
struct is_expression : std::false_type {};

///////////////////////////////////////////////////////////////////////////////

template<typename Field>
struct field_expr;

template <typename Type, template<class> class PatchField, typename GeoMesh>
struct field_expr<Foam::GeometricField<Type, PatchField, GeoMesh>> {

    using value_type = Type;
    static constexpr loc location = Foam::isVolMesh<GeoMesh>::value ? loc::cell : loc::face;
    static constexpr bool has_surface_integrate = false;
    using field_type = Foam::GeometricField<Type, PatchField, GeoMesh>;

    const field_type& field;

    field_expr(const field_type& f)
        : field(f)
    {}

    value_type operator [](Foam::label facei) const {
        return field[facei];
    }

    value_type on_boundary(Foam::label patchi, Foam::label facei) const {
        return field.boundaryField()[patchi][facei];
    }

    const Foam::fvMesh& mesh() const {
        return field.mesh();
    }

    Foam::dimensionSet dimensions() const {
        return field.dimensions();
    }
};

template <typename Field>
struct is_expression<field_expr<Field>> : std::true_type {};

template <typename Type, template<class> class PatchField, typename GeoMesh>
auto read(const Foam::GeometricField<Type, PatchField, GeoMesh>& f) -> field_expr<Foam::GeometricField<Type, PatchField, GeoMesh>> {
    return {f};
}

///////////////////////////////////////////////////////////////////////////////

template<typename CellExpression>
struct linear_interpolate_expr {
    using value_type = typename CellExpression::value_type;
    static_assert(CellExpression::location == loc::cell, "Argument ot interpolation must be cell expression");
    static constexpr loc location = loc::face;
    static constexpr bool has_surface_integrate = CellExpression::has_surface_integrate;

    CellExpression nested;
    const Foam::surfaceScalarField& weight;
    const Foam::labelUList& owner;
    const Foam::labelUList& neighbour;

    linear_interpolate_expr(CellExpression expr)
        : nested(expr)
        , weight(nested.mesh().weights())
        , owner(nested.mesh().owner())
        , neighbour(nested.mesh().neighbour())
    {}

    value_type operator[](Foam::label facei) const {
        auto w = weight[facei];
        auto own = owner[facei];
        auto nei = neighbour[facei];
        return w*nested[own] + (1-w)*nested[nei];
    }

    value_type on_boundary(Foam::label patchi, Foam::label facei) const {
        return nested.on_boundary(patchi, facei);
    }

    const Foam::fvMesh& mesh() const {
        return nested.mesh();
    }

    Foam::dimensionSet dimensions() const {
        return nested.dimensions();
    }
};

template <typename Expr>
struct is_expression<linear_interpolate_expr<Expr>> : std::true_type {};

template <typename Expression, typename std::enable_if<is_expression<Expression>::value, int>::type = 0>
auto interpolate(Expression e) -> linear_interpolate_expr<Expression> {
    return {e};
}

///////////////////////////////////////////////////////////////////////////////

struct subscript {
    Foam::label i;

    template<typename Expr>
    auto operator()(const Expr& e) const -> typename Expr::value_type {
        return e[i];
    }
};

struct boundary {
    Foam::label patchi;
    Foam::label facei;

    template<typename Expr>
    auto operator()(const Expr& e) const -> typename Expr::value_type {
        return e.on_boundary(patchi, facei);
    }
};

template <typename Expr>
struct is_vol_expr {
    static constexpr bool value = (Expr::location == loc::cell);
};

template <typename Expr>
struct is_surface_expr {
    static constexpr bool value = (Expr::location == loc::face);
};

///////////////////////////////////////////////////////////////////////////////

#define FVE_UNARY_FUNCTION(name, dimFunc)                             \
template<typename Expression>                                \
struct name##_expr {                                         \
    using value_type = decltype(name(std::declval<typename Expression::value_type>())); \
    static constexpr loc location = Expression::location;    \
    static constexpr bool has_surface_integrate = Expression::has_surface_integrate; \
    Expression nested;                                       \
    name##_expr(Expression e)                                \
        : nested(e)                                          \
    {}                                                       \
    value_type operator [](Foam::label i) const {            \
        return name(nested[i]);                              \
    }                                                        \
    value_type on_boundary(Foam::label patchi, Foam::label facei) const { \
        return name(nested.on_boundary(patchi, facei));                   \
    }                                                                     \
    const fvMesh& mesh() const {                             \
        return nested.mesh();                                \
    }                                                        \
    Foam::dimensionSet dimensions() const {                  \
        return dimFunc(nested.dimensions());                    \
    }                                                        \
};                                                           \
template <typename Expression, typename std::enable_if<is_expression<Expression>::value, int>::type = 0> \
auto name(Expression e) -> name##_expr<Expression> {         \
    return {e};                                              \
}                                                            \
template <typename Expression>                               \
struct is_expression<name##_expr<Expression>> : std::true_type {};

FVE_UNARY_FUNCTION(twoSymm, transform)
FVE_UNARY_FUNCTION(dev, transform)
FVE_UNARY_FUNCTION(sqr, sqr)
FVE_UNARY_FUNCTION(mag, mag)
FVE_UNARY_FUNCTION(magSq, magSqr)
FVE_UNARY_FUNCTION(pow2, pow2)
FVE_UNARY_FUNCTION(pow3, pow3)
FVE_UNARY_FUNCTION(pow4, pow4)
FVE_UNARY_FUNCTION(pow5, pow5)
FVE_UNARY_FUNCTION(pow6, pow6)
FVE_UNARY_FUNCTION(pow025, pow025)
FVE_UNARY_FUNCTION(sqrt, sqrt)
FVE_UNARY_FUNCTION(cbrt, cbrt)

///////////////////////////////////////////////////////////////////////////////

#define FVE_UNARY_OPERATOR(name, op)                         \
template<typename Expression>                                \
struct name##_expr {                                         \
    using value_type = decltype(op std::declval<typename Expression::value_type>()); \
    static constexpr loc location = Expression::location;    \
    static constexpr bool has_surface_integrate = Expression::has_surface_integrate; \
    Expression nested;                                       \
    name##_expr(Expression e)                                \
    : nested(e)                                              \
    {}                                                       \
    value_type operator [](Foam::label i) const {            \
        return op nested[i];                                 \
    }                                                        \
    value_type on_boundary(Foam::label patchi, Foam::label facei) const { \
        return op nested.on_boundary(patchi, facei);                      \
    }                                                                     \
    const fvMesh& mesh() const {                             \
        return nested.mesh();                                \
    }                                                        \
    Foam::dimensionSet dimensions() const {                  \
        return op nested.dimensions();                       \
    }                                                        \
};                                                           \
template <typename Expression, typename std::enable_if<is_expression<Expression>::value, int>::type = 0> \
auto operator op(Expression e) -> name##_expr<Expression> {  \
    return {e};                                              \
}                                                            \
template <typename Expression>                               \
struct is_expression<name##_expr<Expression>> : std::true_type {};

FVE_UNARY_OPERATOR(neg, -)

///////////////////////////////////////////////////////////////////////////////

#define FVE_BINARY_OPERATOR(name, op)                                                    \
template<typename Expression1, typename Expression2>                                     \
struct name##_expr {                                                                     \
    using value_type = decltype(std::declval<typename Expression1::value_type>() op std::declval<typename Expression2::value_type>()); \
    static_assert(Expression1::location == Expression2::location, "Operands to binary expression must have the same location"); \
    static constexpr loc location = Expression1::location;                               \
    static constexpr bool has_surface_integrate = Expression1::has_surface_integrate || Expression2::has_surface_integrate; \
    Expression1 lhs;                                                                     \
    Expression2 rhs;                                                                     \
        name##_expr(Expression1 e1, Expression2 e2)                                      \
        : lhs(e1)                                                                        \
        , rhs(e2)                                                                        \
    {}                                                                                   \
    value_type operator [](Foam::label i) const {                                        \
            return lhs[i] op rhs[i];                                                     \
    }                                                                                    \
    value_type on_boundary(Foam::label patchi, Foam::label facei) const {                \
        return lhs.on_boundary(patchi, facei) op rhs.on_boundary(patchi, facei);         \
    }                                                                                    \
    const fvMesh& mesh() const {                                                         \
        return lhs.mesh();                                                               \
    }                                                                                    \
    Foam::dimensionSet dimensions() const {                                              \
        return lhs.dimensions() op rhs.dimensions();                                     \
    }                                                                                    \
};                                                                                       \
template <typename Expr1, typename Expr2, \
typename std::enable_if<is_expression<Expr1>::value && is_expression<Expr2>::value , int>::type = 0> \
auto operator op(Expr1 e1, Expr2 e2) -> name##_expr<Expr1, Expr2> {                       \
        return {e1, e2};                                                                 \
}                                                                                        \
template <typename Expr1, typename Expr2>                                                \
struct is_expression<name##_expr<Expr1, Expr2>> : std::true_type {};


FVE_BINARY_OPERATOR(add, +)
FVE_BINARY_OPERATOR(sub, -)
FVE_BINARY_OPERATOR(mul, *)
FVE_BINARY_OPERATOR(divide, /)
FVE_BINARY_OPERATOR(dot, &)
FVE_BINARY_OPERATOR(cross, ^)
FVE_BINARY_OPERATOR(dotdot, &&)

///////////////////////////////////////////////////////////////////////////////

template <typename Type, template<class> class PatchField, typename GeoMesh, typename Expression,
         typename std::enable_if<!Expression::has_surface_integrate, int>::type = 0>
[[gnu::noinline]]
void operator<<=(Foam::GeometricField<Type, PatchField, GeoMesh>& f, Expression e)
{
    static_assert(Expression::location == (Foam::isVolMesh<GeoMesh>::value ? loc::cell : loc::face),
                  "Expression must have same location (cell or face) as the target field");
    const fvMesh& mesh = e.mesh();

    if (!e.dimensions().dimensionless()) {
        // Check dimensions
        f.dimensions() = e.dimensions();
    }

    Foam::label nInternalElems = f.internalField().size();
    for (Foam::label i = 0; i < nInternalElems; i++) {
        f[i] = e[i];
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

} // namespace fve
} // namespace Foam
