from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

from functions.pipeline_new_mask.report_plots import draw_single_seg_panel


def test_draw_single_seg_panel_renders_bars_and_labels():
    df = pd.DataFrame(
        [
            {"segment": "anterior", "target": "olfactory_cortex", "subject": "sub-a", "sift2_density": 1.0},
            {"segment": "anterior", "target": "olfactory_cortex", "subject": "sub-b", "sift2_density": 2.0},
            {"segment": "anterior", "target": "striatum", "subject": "sub-a", "sift2_density": 3.0},
            {"segment": "anterior", "target": "striatum", "subject": "sub-b", "sift2_density": 4.0},
            {"segment": "anterior", "target": "thalamus", "subject": "sub-a", "sift2_density": 5.0},
            {"segment": "anterior", "target": "thalamus", "subject": "sub-b", "sift2_density": 6.0},
            {"segment": "anterior", "target": "cerebellum", "subject": "sub-a", "sift2_density": 7.0},
            {"segment": "anterior", "target": "cerebellum", "subject": "sub-b", "sift2_density": 8.0},
        ]
    )

    targets_order = ["olfactory_cortex", "striatum", "thalamus", "cerebellum"]
    target_labels = {
        "olfactory_cortex": "Paleocortex",
        "striatum": "Striatum",
        "thalamus": "Thalamus",
        "cerebellum": "Cerebellum",
    }
    plot_subjects = df["subject"].unique()

    fig, ax = plt.subplots(figsize=(8, 5))
    draw_single_seg_panel(
        ax,
        df,
        "anterior",
        "#E24A33",
        "Anterior Graft → Destinations",
        targets_order=targets_order,
        target_labels=target_labels,
        plot_subjects=plot_subjects,
        marker_list=["o", "s"],
    )

    assert len(ax.patches) == 4
    assert [tick.get_text() for tick in ax.get_xticklabels()] == [
        "Paleocortex",
        "Striatum",
        "Thalamus",
        "Cerebellum",
    ]
    assert ax.get_title() == "Anterior Graft → Destinations"
    assert ax.get_legend() is not None

    plt.close(fig)