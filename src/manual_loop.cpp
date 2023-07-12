/*
 *
 * Copyright: 2023 Ilya Popov <ilya.popov@isteq.nl>, ISTEQ BV
 *
 * SPDX-License-Identifier: GPL3.0-or-later
 *
 */

#include "manual_loop.hpp"

void Foam::compute_viscous_flux(surfaceVectorField &F_rhoU, const volTensorField &gradU, const volScalarField &mu)
{
    const fvMesh& mesh = gradU.mesh();
    const surfaceVectorField& Sf = mesh.Sf();
    const labelUList& owner = mesh.owner();
    const labelUList& neighbour = mesh.neighbour();
    const surfaceScalarField& weights = mesh.weights();

    for (label facei = 0; facei < mesh.nInternalFaces(); facei++) {
        label own = owner[facei];
        label nei = neighbour[facei];
        scalar w = weights[facei];

        scalar mu_f    = w*mu[own]    + (1.0-w)*mu[nei];
        tensor gradU_f = w*gradU[own] + (1.0-w)*gradU[nei];
        F_rhoU[facei] = (mu_f * dev(twoSymm(gradU_f))) & Sf[facei];
    }

    for (label patchi = 0; patchi < mesh.boundary().size(); patchi++) {
        const fvPatch& patch = mesh.boundary()[patchi];
        const fvPatchTensorField& p_gradU = gradU.boundaryField()[patchi];
        const fvsPatchVectorField& p_Sf = Sf.boundaryField()[patchi];
        const fvPatchScalarField& p_mu = mu.boundaryField()[patchi];
        fvsPatchVectorField& p_F_rhoU = F_rhoU.boundaryFieldRef()[patchi];

        if (!patch.coupled()) {
            label nFaces = patch.size();
            for (label facei = 0; facei < nFaces; facei++) {
                p_F_rhoU[facei] = (p_mu[facei] * dev(twoSymm(p_gradU[facei]))) & p_Sf[facei];
            }
        }
        else {
            auto p_mu_internal = p_mu.patchInternalField();
            auto p_mu_neighbour = p_mu.patchNeighbourField();
            auto p_gradU_internal = p_gradU.patchInternalField();
            auto p_gradU_neighbour = p_gradU.patchNeighbourField();
            const fvsPatchScalarField& p_w = weights.boundaryField()[patchi];

            label nFaces = patch.size();
            for (label facei = 0; facei < nFaces; facei++) {
                scalar w = p_w[facei];
                scalar mu_f    = w*p_mu_internal()[facei]    + (1.0-w)*p_mu_neighbour()[facei];
                tensor gradU_f = w*p_gradU_internal()[facei] + (1.0-w)*p_gradU_neighbour()[facei];
                p_F_rhoU[facei] = (mu_f * dev(twoSymm(gradU_f))) & p_Sf[facei];
            }
        }
    }
}


