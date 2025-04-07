#!/usr/bin/env python3

# Copyright: 2023 Ilya Popov <ilya.popov@isteq.nl>, ISTEQ BV
# SPDX-License-Identifier: GPL3.0-or-later

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from cycler import cycler

case_dir = Path("test_case")
data_dir = case_dir / "benchmarking"
plots_dir = Path("plots")
plots_dir.mkdir(parents=True, exist_ok=True)

methods = [
    ("Standard OpenFOAM", "Standard OpenFOAM"), 
    ("for_each_face", "for_each_face_interp"), 
    ("Expression templates", "expression_templates"),
#    "grad_expr",
    ("Expr.templ + grad", "grad_expr_2"),
]

series = [
#    ("desktop", data_dir, "results-[1-9]*.csv"),
#    ("desktop renumbered", data_dir/"desktop-block4k", "results-[1-9]*.csv"),
#    ("server renumbered", data_dir/"server-block4k", "results-[1-9]*.csv"),
    ("4k", data_dir/"server2-block4k", "results-[1-9]*.csv"),
    ("8k", data_dir/"server2-block8k", "results-[1-9]*.csv"),
    ("16k", data_dir/"server2-block16k", "results-[1-9]*.csv"),
    ("32k", data_dir/"server2-block32k", "results-[1-9]*.csv"),
#    ("laptop renumbered", data_dir/"laptop-block4k", "results-[1-9]*.csv"),
]

colors = cycler(color=["C0", "C1", "C3", "C4", "C5", "C6", "C7"])
styles = cycler(linestyle=["-", "-.", "--", ":", "-", "-.", "--", ":"])

def make_plots(methods, method_styles, series, series_styles, suffix):
    results = list()
    for series_label, directory, glob in series:

        data = pd.concat((pd.read_csv(filename, sep=";") for filename in directory.glob(glob)))
        data.sort_values("batch", inplace=True)

        results.append((series_label, data))

    # plt.clf()
    # for (method_label, method), color in zip(methods, colors):
    #     for (series_label, data), style in zip(results, styles):
    #         selected = data.loc[data.name == method]
    #         if not selected.empty:
    #             plt.plot(selected.batch, selected.elapsed, label=f"{method} ({series_label})", **color, **style)

    # plt.xlabel("Number of faces")
    # plt.ylabel("Time, s")
    # plt.xscale("log")
    # plt.yscale("log")
    # #plt.legend()
    # plt.grid()

    # plt.savefig(plots_dir / "runtime.pdf", metadata={"CreationDate": None})
    plt.clf()
    for (method_label, method), color in zip(methods, method_styles):
        need_label = True
        for (series_label, data), style in zip(results, series_styles):
            selected = data.loc[data.name == method]
            if not selected.empty:
                plt.plot(selected.batch, 1e-6*selected.batch/selected.elapsed, 
                    label= f"{method_label} {series_label}" if need_label else None, 
                    **color, **style)
                need_label = False

    plt.xlabel("Number of faces")
    plt.ylabel("Throughput, Mfaces/s")
    plt.xscale("log")
    plt.legend()
    plt.grid()
    plt.ylim(0, 120)

    plt.savefig(plots_dir / f"throughput{suffix}.pdf", metadata={"CreationDate": None})

make_plots(
    methods = [
            ("Standard OpenFOAM", "Standard OpenFOAM"), 
    ],
    method_styles=colors,
    series= [
            ("", data_dir/"server2-block4k", "results-[1-9]*.csv"),
    ],
    series_styles=styles,
    suffix="-OF"
)

make_plots(
    methods = [
            ("Standard OpenFOAM", "Standard OpenFOAM"), 
            ("for_each_face", "for_each_face_interp"), 
    ],
    method_styles=colors,
    series= [
            ("", data_dir/"server2-block4k", "results-[1-9]*.csv"),
    ],
    series_styles=styles,
    suffix="-for_each_face"
)

make_plots(
    methods = [
            ("Standard OpenFOAM", "Standard OpenFOAM"), 
            ("for_each_face", "for_each_face_interp"), 
            ("Expression templates", "expression_templates"),
    ],
    method_styles=colors,
    series= [
            ("", data_dir/"server2-block4k", "results-[1-9]*.csv"),
    ],
    series_styles=styles,
    suffix="-expr_templates"
)

make_plots(
    methods = [
            ("Standard OpenFOAM", "Standard OpenFOAM"), 
            ("for_each_face", "for_each_face_interp"), 
            ("Expression templates", "expression_templates"),
            ("With grad expr", "grad_expr_2"),
    ],
    method_styles=colors,
    series= [
            ("", data_dir/"server2-block4k", "results-[1-9]*.csv"),
    ],
    series_styles=styles,
    suffix="-grad"
)

make_plots(
    methods = [
#            ("Standard OpenFOAM", "Standard OpenFOAM"), 
#            ("for_each_face", "for_each_face_interp"), 
#            ("Expression templates", "expression_templates"),
            ("With grad expr", "grad_expr_2"),
    ],
    method_styles=colors,
    series= [
            ("4k", data_dir/"server2-block4k", "results-[1-9]*.csv"),
            ("8k", data_dir/"server2-block8k", "results-[1-9]*.csv"),
            ("16k", data_dir/"server2-block16k", "results-[1-9]*.csv"),
            ("32k", data_dir/"server2-block32k", "results-[1-9]*.csv"),
            ("64k", data_dir/"server2-block64k", "results-[1-9]*.csv"),
            ("128k", data_dir/"server2-block128k", "results-[1-9]*.csv"),
    ],
    series_styles=styles,
    suffix="-block_size"
)

make_plots(
    methods = [
            ("Standard OpenFOAM", "Standard OpenFOAM"), 
            ("for_each_face", "for_each_face_interp"), 
            ("Expression templates", "expression_templates"),
            ("With grad expr", "grad_expr_2"),
    ],
    method_styles=colors,
    series= [
            ("desktop", data_dir/"desktop-block4k", "results-[1-9]*.csv"),
            ("server 1", data_dir/"server-block4k", "results-[1-9]*.csv"),
            ("server 2", data_dir/"server2-block4k", "results-[1-9]*.csv"),
    ],
    series_styles=styles,
    suffix="-architecture"
)
