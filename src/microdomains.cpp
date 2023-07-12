/*
 *
 * Copyright: 2023 Ilya Popov <ilya.popov@isteq.nl>, ISTEQ BV
 *
 * SPDX-License-Identifier: GPL3.0-or-later
 *
 */

#include "microdomains.hpp"

#include "volFields.H"
#include "defineDebugSwitch.H"

namespace Foam {
namespace fve {

defineTypeNameAndDebug(microdomains, 0);

} // namespace fve
} // namespace Foam

Foam::fve::microdomains::microdomains(const fvMesh &mesh)
    : Foam::MeshObject<Foam::fvMesh, Foam::GeometricMeshObject, microdomains>(mesh)
{
    Foam::volScalarField cellDist(Foam::IOobject("cellDist", mesh.time().constant(), mesh,
                                                 Foam::IOobject::MUST_READ, Foam::IOobject::NO_WRITE),
                                  mesh);
    Foam::volScalarField origCellID(Foam::IOobject("origCellID", mesh.time().timeName(), mesh,
                                                   Foam::IOobject::MUST_READ, Foam::IOobject::NO_WRITE),
                                    mesh);

    cell_dist.assign(mesh.nCells(), -1);

    Foam::Info << "Assigning cells to microdomains...\n";

    int current_domain = -1;
    for (Foam::label celli = 0; celli < mesh.nCells(); celli++) {
        Foam::label orig_celli = static_cast<Foam::label>(origCellID[celli]);
        Foam::label d = static_cast<Foam::label>(cellDist[orig_celli]);
        cell_dist[celli] = d;
        if (d != current_domain) {
            //                Foam::Info << "New microdomain " << d << " starting at cell " << celli << " origCellID " << orig_celli << "\n";
            if (d != domains.size()) {
                Foam::FatalError << "Microdomain " << d << " is out of order. Microdomains are not ordered properly" << Foam::abort(Foam::FatalError);
            }
            if (current_domain >= 0) {
                domains[current_domain].cells.b = celli;
            }
            domains.emplace_back(microdomain{{celli, mesh.nCells()}, {-1, -1}, {-1, -1}});
            current_domain = d;
        }
    }

    Foam::Info << "Assigning faces to microdomains...\n";

    std::vector<std::vector<Foam::label>> microdomain_internal_faces(domains.size());
    std::vector<std::vector<Foam::label>> microdomain_boundary_faces(domains.size());

    for (Foam::label facei = 0; facei < mesh.nInternalFaces(); facei++) {
        Foam::label own = mesh.owner()[facei];
        Foam::label md_own = cell_dist[own];
        Foam::label nei = mesh.neighbour()[facei];
        Foam::label md_nei = cell_dist[nei];
        if (md_own == md_nei) {
            microdomain_internal_faces[md_own].push_back(facei);
        }
        else if (md_own < md_nei) {
            //Foam::Info << facei << ' ' << md_own << ' ' << md_nei << '\n';
            microdomain_boundary_faces[md_own].push_back(facei);
        }
        else {
            Foam::FatalError << "Ordering is wrong." << Foam::abort(Foam::FatalError);
        }
    }

    for (size_t d = 0; d < domains.size(); ++d) {
        //            Foam::Info << "Processing faces for domain " << d << '\n';
        auto& md = domains[d];
        const auto& int_faces = microdomain_internal_faces[d];
        md.internal_faces = to_range(int_faces);
        md.own_boundary_faces = to_range(microdomain_boundary_faces[d]);
    }
}

Foam::fve::microdomains::~microdomains()
{

}
