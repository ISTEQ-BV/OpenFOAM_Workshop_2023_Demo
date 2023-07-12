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
#include "surfaceFields.H"
#include "volFields.H"

#include <tuple>
#include <functional>
#include "boost/mp11/tuple.hpp"

template <typename Field>
struct value_type;

template <typename Type>
struct value_type<Foam::fvsPatchField<Type>> {
    using type = Type;
};

template <typename Type>
struct value_type<Foam::fvPatchField<Type>> {
    using type = Type;
};

template <typename Type, template<class> class PatchField, typename GeoMesh>
struct value_type<Foam::GeometricField<Type, PatchField, GeoMesh>> {
    using type = Type;
};

template <typename Field>
struct patch_type;

template <typename Type, template<class> class PatchField, typename GeoMesh>
struct patch_type<Foam::GeometricField<Type, PatchField, GeoMesh>> {
    using type = PatchField<Type>;
};



struct get_patch
{
    int patchi;

    get_patch(int patchi)
        :patchi(patchi)
    {}

//    template <typename Type, template<class> class PatchField, typename GeoMesh>
//    PatchField<Type>& operator()(Foam::GeometricField<Type, PatchField, GeoMesh>&& f) const {
//        return f.boundaryField()[patchi];
//    }
    template <typename Field>
    auto operator()(const Field& f) const -> const typename patch_type<Field>::type& {
        return f.boundaryField()[patchi];
    }
    template <typename Field>
    auto operator()(Field& f) const -> typename patch_type<Field>::type& {
        return f.boundaryFieldRef()[patchi];
    }
};

struct get_index
{
    int i;

    get_index(int i)
        : i(i)
    {}

    template <typename Patch>
    const typename value_type<Patch>::type& operator()(const Patch& f) const {
        return f[i];
    }
    template <typename Patch>
    typename value_type<Patch>::type& operator()(Patch& f) const {
        return f[i];
    }
};

struct get_face
{
    int facei;

    get_face(int i)
        : facei(i)
    {}

    template <typename FieldLike>
    auto operator()(const FieldLike& f) const -> const typename value_type<FieldLike>::type {
        return f[facei];
    }

    //template <typename FieldLike>
    //auto operator()(FieldLike& f) const -> typename value_type<FieldLike>::type& {
    //    return f[facei];
    //}

    template <typename Type, template<class> class PatchField>
    auto operator()(const Foam::GeometricField<Type, PatchField, Foam::surfaceMesh>& f) const -> const Type& {
        return f[facei];
    }
    template <typename Type, template<class> class PatchField>
    auto operator()(Foam::GeometricField<Type, PatchField, Foam::surfaceMesh>& f) const -> Type& {
        return f[facei];
    }
};

struct get_face_interp
{
    int facei;
    int own;
    int nei;
    Foam::scalar w;

    get_face_interp(int i, const Foam::labelUList& owner, const Foam::labelUList& neighbour, const Foam::surfaceScalarField& weights)
        : facei(i)
        , own(owner[facei])
        , nei(neighbour[facei])
        , w(weights[facei])
    {}

    template <typename Type, template<class> class PatchField>
    auto operator()(const Foam::GeometricField<Type, PatchField, Foam::surfaceMesh>& f) const -> const Type& {
        return f[facei];
    }
    template <typename Type, template<class> class PatchField>
    auto operator()(Foam::GeometricField<Type, PatchField, Foam::surfaceMesh>& f) const -> Type& {
        return f[facei];
    }
    template <typename Type, template<class> class PatchField>
    auto operator()(const Foam::GeometricField<Type, PatchField, Foam::volMesh>& f) const -> Type {
        return w*f[own] + (1.0 - w)*f[nei];
    }
    template <typename Type, template<class> class PatchField>
    auto operator()(Foam::GeometricField<Type, PatchField, Foam::volMesh>& f) const -> Type {
        return w*f[own] + (1.0 - w)*f[nei];
    }
};

template <typename Callable, typename... Fields>
[[gnu::noinline]]
void for_each_face_interp(Callable&& func, Fields&& ...fs) {
    auto fields = std::forward_as_tuple(fs...);
    const Foam::fvMesh& mesh = std::get<0>(fields).mesh();
    const auto& weights = mesh.weights();
    const auto& owner = mesh.owner();
    const auto& neighbour = mesh.neighbour();

    auto nInternalFaces = mesh.nInternalFaces();
    for (int facei = 0; facei < nInternalFaces; ++facei) {
        auto g = get_face_interp(facei, owner, neighbour, weights);
        func(g(fs)...);
    }

    for (int patchi = 0; patchi < mesh.boundary().size(); ++patchi) {
        const Foam::fvPatch& patch = mesh.boundary()[patchi];
        auto patchfields = boost::mp11::tuple_transform(get_patch(patchi), fields);
        int nFaces = patch.size();
        for (int facei = 0; facei < nFaces; ++facei) {
            auto values = boost::mp11::tuple_transform(get_index(facei), patchfields);
            boost::mp11::tuple_apply(func, values);
        }
    }
};

//template <typename Callable, typename... Fields>
//[[gnu::noinline]]
//void for_each_face(Callable&& func, Fields&& ...fs) {
//    auto fields = std::forward_as_tuple(fs...);
//    const Foam::fvMesh& mesh = std::get<0>(fields).mesh();

//    auto nInternalFaces = mesh.nInternalFaces();
//    for (int facei = 0; facei < nInternalFaces; ++facei) {
//        auto g = get_face(facei);
//        func(g(fs)...);
//    }

//    for (int patchi = 0; patchi < mesh.boundary().size(); ++patchi) {
//        const Foam::fvPatch& patch = mesh.boundary()[patchi];
//        auto patchfields = boost::mp11::tuple_transform(get_patch(patchi), fields);
//        int nFaces = patch.size();
//        for (int facei = 0; facei < nFaces; ++facei) {
//            auto values = boost::mp11::tuple_transform(get_index(facei), patchfields);
//            boost::mp11::tuple_apply(func, values);
//        }
//    }
//};

template <typename Callable, typename... Fields>
[[gnu::noinline]]
void for_each_cell(Callable&& func, Fields& ...fs) {
    auto fields = std::forward_as_tuple(fs...);
    const Foam::fvMesh& mesh = std::get<0>(fields).mesh();

    auto nCells = mesh.nCells();
    for (int celli = 0; celli < nCells; ++celli) {
        func(fs[celli]...);
    }

    for (int patchi = 0; patchi < mesh.boundary().size(); ++patchi) {
        const Foam::fvPatch& patch = mesh.boundary()[patchi];
        auto patchfields = boost::mp11::tuple_transform(get_patch(patchi), fields);
        int nFaces = patch.size();
        for (int facei = 0; facei < nFaces; ++facei) {
            auto values = boost::mp11::tuple_transform(get_index(facei), patchfields);
            boost::mp11::tuple_apply(func, values);
        }
    }
};

template <typename Callable, typename... Fields>
[[gnu::noinline]]
void for_each_cell_neighbours(Callable&& func, Fields& ...fs) {
    auto fields = std::forward_as_tuple(fs...);
    const Foam::fvMesh& mesh = std::get<0>(fields).mesh();
    const auto& weights = mesh.weights();
    const auto& owner = mesh.owner();
    const auto& neighbour = mesh.neighbour();

    auto nCells = mesh.nCells();
    for (int celli = 0; celli < nCells; ++celli) {
        const Foam::cell& c = mesh.cells()[celli];

        for (auto facei: c) {
            auto own = owner[facei];
            auto nei = neighbour[facei];
            auto other = (own == celli) ? nei : own;
            auto k =     (own == celli) ?   1 :  -1;

            func(celli, other, k, fs[celli]..., fs[other]...);
        }
    }

    for (int patchi = 0; patchi < mesh.boundary().size(); ++patchi) {
        const Foam::fvPatch& patch = mesh.boundary()[patchi];
        auto patchfields = boost::mp11::tuple_transform(get_patch(patchi), fields);
        int nFaces = patch.size();
        for (int facei = 0; facei < nFaces; ++facei) {
            auto values = boost::mp11::tuple_transform(get_index(facei), patchfields);
            boost::mp11::tuple_apply(func, values);
        }
    }
};

inline void test_1(const Foam::surfaceScalarField& f1, const Foam::surfaceVectorField& f2, Foam::surfaceVectorField& result)
{
    for_each_face_interp([](Foam::scalar f1, Foam::vector f2, Foam::vector& result){
        result = f1 * f2;
    }, f1, f2, result);
}
