/*
 *
 * Copyright: 2023 Ilya Popov <ilya.popov@isteq.nl>, ISTEQ BV
 *
 * SPDX-License-Identifier: GPL3.0-or-later
 *
 */

#pragma once

#include "fvMesh.H"
#include "MeshObject.H"

#include <algorithm>
#include <vector>

namespace Foam {
namespace fve {

struct index_range {
    struct iterator {
        int i;

        int operator*() const {return i;}
        iterator& operator++() {++i; return *this;}
        iterator operator++(int) {auto prev = *this; ++i; return prev;}
        bool operator==(const iterator& other) const {return i == other.i; }
        bool operator!=(const iterator& other) const {return i != other.i; }
    };

    int a;
    int b;

    iterator begin() const {return iterator{a}; }
    iterator end() const {return iterator{b}; }
    size_t size() const {return b - a; }
    bool empty() const {return a == b; }
    int operator[](size_t i) const {return a + i; }
    int front() const { return a; }
    int back() const { return b-1; }

    friend std::ostream& operator<<(std::ostream& s, const index_range& r)
    {
        return s << "{" << r.a << ", " << r.b << "}";
    }
};

struct microdomain {
    index_range cells;
    index_range internal_faces;
    index_range own_boundary_faces;
};

inline
bool is_contiguous(const std::vector<int>& v) {
    return std::is_sorted(v.begin(), v.end()) && v.back() - v.front() + 1 == v.size();
}

inline
index_range to_range(const std::vector<int>& v) {
    if (v.empty()) {
        return {-1, -1};
    }
    if (!is_contiguous(v)) {
        for (auto facei: v) {
            Foam::Info << facei << '\n';
        }
        Foam::FatalError << " start_face " << v.front() << " end_face " << v.back() << " n_faces " << v.size()
                         << " Faces for microdomains are not contiguous" << Foam::abort(Foam::FatalError);
    }
    return {v.front(), v.back()+1};
}

struct microdomains : public Foam::MeshObject<Foam::fvMesh, Foam::GeometricMeshObject, microdomains> {
    TypeName("microdomains");

    std::vector<int> cell_dist;
    std::vector<fve::microdomain> domains;

    explicit microdomains(const Foam::fvMesh& mesh);
    virtual ~microdomains();
};

} //namespace fve
} //namespace Foam
