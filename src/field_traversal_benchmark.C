/*
 *
 * Copyright: 2023 Ilya Popov <ilya.popov@isteq.nl>, ISTEQ BV
 *
 * SPDX-License-Identifier: GPL3.0-or-later
 *
 */

#include "expressions.hpp"
#include "map_expr.hpp"
#include "microdomains.hpp"
#include "traversal.hpp"
#include "grad_expr.hpp"
#include "process_microdomains.hpp"
#include "grad_expr_2.hpp"
#include "manual_loop.hpp"

#define ANKERL_NANOBENCH_IMPLEMENT
#include "nanobench.h"

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    volScalarField rho  ( IOobject ( "rho" , runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE  )
                       , mesh, dimensionedScalar("", dimDensity, 1.0) );
    volScalarField mu  ( IOobject ( "mu" , runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE  )
                      , mesh, dimensionedScalar("", dimensionSet(1, -1, -1, 0, 0), 1.0) );
    volVectorField U   (IOobject("U", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE)
                     , mesh, dimensionedVector("", dimensionSet(0, 1, -1, 0, 0), {1.0, 0.0, 0.0})
                     );
    volTensorField gradU (IOobject("gradU", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::AUTO_WRITE)
                         , mesh, dimensionedTensor("", dimensionSet(0, 0, -1, 0, 0), Foam::Zero)
                         );
    surfaceVectorField F_rhoU (IOobject("F_rhoU", runTime.timeName(), mesh, IOobject::NO_READ, IOobject::NO_WRITE),
                              (fvc::interpolate(rho*U) & mesh.Sf())*fvc::interpolate(U));

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    const fve::microdomains& mds = fve::microdomains::New(mesh);

    ankerl::nanobench::Bench b;
    b.title("Computing viscous flux")
        .unit("face")
        .batch(mesh.nFaces())
        .warmup(3)
        .minEpochIterations(5)
        .relative(true);
    b.performanceCounters(true);

    b.run("Standard OpenFOAM", [&] {

        F_rhoU = (fvc::interpolate(mu) * dev(twoSymm(fvc::interpolate(fvc::grad(U))))) & mesh.Sf();

    });

    b.run("manual loop", [&] {

        volTensorField gradU(fvc::grad(U));

        compute_viscous_flux(F_rhoU, gradU, mu);
    });


    b.run("for_each_face_interp", [&] {

        volTensorField gradU(fvc::grad(U));

        for_each_face_interp(
            [](const tensor& gradU_f, const scalar& mu_f, const vector& s_f, vector& F_rhoU_f) {
                F_rhoU_f = (mu_f * dev(twoSymm(gradU_f))) & s_f;
            },
            gradU, mu, mesh.Sf(), F_rhoU);

    });

    b.run("expression_templates", [&] {

        volTensorField gradU(fvc::grad(U));

        F_rhoU <<= (interpolate(fve::read(mu)) * dev(twoSymm(interpolate(fve::read(gradU))))) & fve::read(mesh.Sf());

    });

    b.run("map", [&] {

        volTensorField gradU(fvc::grad(U));

        F_rhoU <<= map(
            [](const tensor& gradU_f, const scalar& mu_f, const vector& s_f) {
                return (mu_f * dev(twoSymm(gradU_f))) & s_f;
            },
            interpolate(fve::read(gradU)), interpolate(fve::read(mu)), fve::read(mesh.Sf())
            );

    });

    b.run("grad_expr", [&] {

        F_rhoU <<= (interpolate(fve::read(mu)) * dev(twoSymm(interpolate(grad(fve::read(U)))))) & fve::read(mesh.Sf());

    });

    b.run("grad_expr_2", [&] {

        F_rhoU <<= (interpolate(fve::read(mu)) * dev(twoSymm(interpolate(grad(gradU, fve::read(U)))))) & fve::read(mesh.Sf());

    });


//    b.run("naive_grad_expr", [&] {

//        F_rhoU <<= (interpolate(fve::read(mu)) * dev(twoSymm(interpolate(gauss_grad(gradU, fve::read(U)))))) & fve::read(mesh.Sf());

//    });

    b.run("grad only", [&] {

        volTensorField gradU(fvc::grad(U));

    });

    std::ofstream csv("results.csv");
    b.render(ankerl::nanobench::templates::csv(), csv);


    Info << nl;
    runTime.printExecutionTime(Info);
    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
