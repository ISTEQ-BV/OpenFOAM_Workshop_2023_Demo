/*
 *
 * Copyright: 2023 Ilya Popov <ilya.popov@isteq.nl>, ISTEQ BV
 *
 * SPDX-License-Identifier: GPL3.0-or-later
 *
 */

#pragma once

#include "volFields.H"
#include "surfaceFields.H"
#include "fvMesh.H"

namespace Foam {

void compute_viscous_flux(surfaceVectorField & F_rhoU, const volTensorField& gradU, const volScalarField & mu);

} // namespace Foam
