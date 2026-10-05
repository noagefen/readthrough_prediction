"""Reproducible model comparison for the readthrough prediction study.

The analysis uses the ten supplied 80/20 splits as repeated outer test splits and
performs all tuning inside each outer training set with three-fold cross-validation.
Preprocessing is fit independently in every inner and outer training split.
"""

from __future__ import annotations

import json
import os
import re
import warnings
from pathlib import Path

import numpy as np
import optuna
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import scipy.stats as stats
import scipy
import shap
import xgboost as xgb
from optuna.samplers import TPESampler
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import Lasso, Ridge
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import KFold
from sklearn.preprocessing import OneHotEncoder, StandardScaler


SEED = 20260918
INNER_FOLDS = 3
TRIALS = {"ridge": 20, "lasso": 20, "rf": 20, "xgb": 20}

NOVEL_FEATURES = [
    "CAI_30", "CAI_50", "CAI_70",
    "codon_m_6", "codon_m_5", "codon_m_4", "codon_m_3", "codon_m_2", "codon_m_1",
    "similar_to_stop_cnt", "similar_cnt_70", "similar_cnt_50", "similar_cnt_30",
    "stop_win_mfe", "stop_win_mfe_indx",
    "tAI_30", "tAI_50", "tAI_70",
]

METRICS = ["rmse", "mae", "r2", "spearman", "pearson"]

MODEL_LABELS = {
    "ridge": "Ridge",
    "lasso": "Lasso",
    "rf": "Random\nForest",
    "xgb": "XGBoost",
}

FEATURE_PALETTE = {
    "Similar to stop": "#c9ab88",
    "Codon identity": "#9aac8b",
    "Codon usage bias": "#d6b42f",
    "Exit tunnel AA properties": "#bd5f32",
    "MFE": "#d48ea8",
    "Nucleotide identity": "#808080",
    "Stop codon identity": "#8164b8",
    "AA identity": "#37aaa7",
    "First 3'UTR stop": "#8aca1b",
    "3'UTR length": "#bdbdbd",
}


def benjamini_hochberg(p_values):
    """Return Benjamini-Hochberg FDR-adjusted p-values in input order."""
    p_values = np.asarray(p_values, dtype=float)
    order = np.argsort(p_values)
    ranked = p_values[order]
    adjusted = ranked * len(ranked) / np.arange(1, len(ranked) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def locate_data() -> tuple[Path, Path, Path]:
    """Find WT25 and randomized-genome inputs locally or on Kaggle."""
    if "__file__" in globals():
        repository_root = Path(__file__).resolve().parents[1]
    else:
        repository_root = Path.cwd()
        if repository_root.name == "notebooks":
            repository_root = repository_root.parent
    candidates = [
        (Path("/kaggle/input/datasets/noagef/wt25-dataset"),
         Path("/kaggle/input/datasets/noagef/randgenome-x"), Path("/kaggle/working")),
        (Path("/kaggle/input/wt25-dataset"), Path("/kaggle/input/randgenome-x"), Path("/kaggle/working")),
        (repository_root / "WT25", repository_root / "data", repository_root),
        (repository_root.parent / "WT25", repository_root.parent / "data", repository_root.parent),
    ]
    for wt_dir, rand_dir, work_dir in candidates:
        if (wt_dir / "Xtbl_wt25.csv").exists() and (rand_dir / "Xtbl_rand_1.csv").exists():
            return wt_dir, rand_dir, work_dir
    raise FileNotFoundError("Could not locate WT25 and randomized-genome datasets")


WT_DIR, RAND_DIR, WORK_DIR = locate_data()
OUT_DIR = WORK_DIR / "nested_cv_results"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def make_preprocessor(X: pd.DataFrame) -> ColumnTransformer:
    categorical = [c for c in X.columns if X[c].dtype == "object" or isinstance(X[c].dtype, pd.CategoricalDtype)]
    numerical = [c for c in X.columns if c not in categorical]
    return ColumnTransformer(
        [
            ("cat", OneHotEncoder(handle_unknown="ignore", sparse_output=False), categorical),
            ("num", StandardScaler(), numerical),
        ],
        remainder="drop",
        sparse_threshold=0,
        verbose_feature_names_out=True,
    )


def clean_feature_name(name: str) -> str:
    return name.split("__", 1)[-1]


def raw_feature_name(encoded_name: str, categorical_columns: list[str]) -> str:
    clean = clean_feature_name(encoded_name)
    for col in sorted(categorical_columns, key=len, reverse=True):
        if clean == col or clean.startswith(col + "_"):
            return col
    return clean


def grouped_xgb_shap(model, transformed, preprocessor, raw_columns, categorical_columns):
    """Return XGBoost SHAP contributions summed back to raw input features."""
    contributions = model.get_booster().predict(
        xgb.DMatrix(transformed), pred_contribs=True
    )[:, :-1]
    encoded_names = preprocessor.get_feature_names_out()
    raw_to_index = {name: i for i, name in enumerate(raw_columns)}
    grouped = np.zeros((contributions.shape[0], len(raw_columns)), dtype=float)
    for encoded_index, encoded_name in enumerate(encoded_names):
        raw_name = raw_feature_name(encoded_name, categorical_columns)
        grouped[:, raw_to_index[raw_name]] += contributions[:, encoded_index]
    return grouped


def display_raw_feature(name: str) -> str:
    """Convert raw feature names to concise labels for interpretation figures."""
    name = re.sub(r"^tAI_(\d+)$", r"tAI, last \1 codons", name)
    name = re.sub(r"^CAI_(\d+)$", r"CAI, last \1 codons", name)
    name = re.sub(r"^codon_m_(\d+)$", r"codon identity at -\1", name)
    name = re.sub(r"^nt_m(\d+)$", r"nucleotide identity at -\1", name)
    name = re.sub(r"^nt_p(\d+)$", r"nucleotide identity at +\1", name)
    name = re.sub(r"^aa_m(\d+)$", r"amino acid identity at -\1", name)
    name = re.sub(r"([+-])0+(\d+)", r"\1\2", name)
    replacements = {
        "nis_stop": "first 3'UTR stop identity",
        "l_utr3": "3'UTR length",
        "stop_win_mfe": "stop-window MFE",
        "stop_win_mfe_indx": "stop-window MFE position",
        "MFE": "3'UTR local MFE",
        "dist_bp": "3'UTR structure distance",
        "similar_to_stop_cnt": "similar-to-stop codon count",
    }
    if name.startswith("similar_cnt_"):
        return f"similar-to-stop count, last {name.rsplit('_', 1)[-1]} codons"
    return replacements.get(name, name.replace("_", " "))


def feature_group(raw: str) -> str:
    """Assign a biological feature group for the Figure 4 color scheme."""
    if raw.startswith("similar_"):
        return "Similar to stop"
    if raw.startswith("codon_m_"):
        return "Codon identity"
    if raw.startswith(("CAI_", "tAI_")):
        return "Codon usage bias"
    if raw.startswith("tunnel_"):
        return "Exit tunnel AA properties"
    if raw in {"MFE", "dist_bp", "stop_win_mfe", "stop_win_mfe_indx"}:
        return "MFE"
    if raw.startswith("nt_"):
        return "Nucleotide identity"
    if raw == "stop_codon":
        return "Stop codon identity"
    if raw.startswith("aa_"):
        return "AA identity"
    if raw == "nis_stop":
        return "First 3'UTR stop"
    if raw == "l_utr3":
        return "3'UTR length"
    return "Nucleotide identity"


def display_encoded_feature(name: str) -> str:
    """Convert an encoded model feature name to a concise plot label."""
    name = re.sub(r"^tAI_(\d+)$", r"tAI, last \1 codons", name)
    name = re.sub(r"^CAI_(\d+)$", r"CAI, last \1 codons", name)
    name = re.sub(r"^codon_m_(\d+)_([A-Z]+)$", r"codon -\1: \2", name)
    name = re.sub(r"^nt_m(\d+)_([A-Z]+)$", r"nt -\1: \2", name)
    name = re.sub(r"^nt_p(\d+)_([A-Z]+)$", r"nt +\1: \2", name)
    name = re.sub(r"^aa_m(\d+)_([A-Z]+)$", r"AA -\1: \2", name)
    name = re.sub(r"^stop_codon_([A-Z]+)$", r"stop codon: \1", name)
    name = re.sub(r"^nis_stop_([A-Z]+)$", r"first 3'UTR stop: \1", name)
    replacements = {
        "nis_stop_nan": "first 3'UTR stop: missing",
        "l_utr3": "3'UTR length",
        "similar_to_stop_cnt": "similar-to-stop count",
        "stop_win_mfe": "stop-window MFE",
        "stop_win_mfe_indx": "stop-window MFE position",
    }
    if name in replacements:
        return replacements[name]
    name = re.sub(r"^similar_cnt_(\d+)$", r"similar-to-stop count, last \1 codons", name)
    return name.replace("_", " ")


def make_model(model_name: str, params: dict, seed: int):
    if model_name == "ridge":
        return Ridge(alpha=params["alpha"])
    if model_name == "lasso":
        return Lasso(alpha=params["alpha"], max_iter=20000, random_state=seed)
    if model_name == "rf":
        return RandomForestRegressor(**params, random_state=seed, n_jobs=-1)
    if model_name == "xgb":
        return xgb.XGBRegressor(
            **params,
            objective="reg:squarederror",
            eval_metric="rmse",
            random_state=seed,
            n_jobs=-1,
            tree_method="hist",
            verbosity=0,
        )
    raise ValueError(model_name)


def suggest_params(model_name: str, trial: optuna.Trial) -> dict:
    if model_name == "ridge":
        return {"alpha": trial.suggest_float("alpha", 0.01, 100.0, log=True)}
    if model_name == "lasso":
        return {"alpha": trial.suggest_float("alpha", 0.0001, 1.0, log=True)}
    if model_name == "rf":
        return {
            "n_estimators": trial.suggest_int("n_estimators", 100, 500, step=25),
            "max_depth": trial.suggest_int("max_depth", 3, 20),
            "min_samples_split": trial.suggest_int("min_samples_split", 2, 10),
            "min_samples_leaf": trial.suggest_int("min_samples_leaf", 1, 10),
            "max_features": trial.suggest_categorical("max_features", ["sqrt", "log2", None]),
        }
    if model_name == "xgb":
        return {
            "reg_lambda": trial.suggest_float("reg_lambda", 0.0, 10.0),
            "reg_alpha": trial.suggest_float("reg_alpha", 0.0, 10.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.5, 1.0),
            "gamma": trial.suggest_float("gamma", 0.0, 5.0),
            "subsample": trial.suggest_float("subsample", 0.5, 1.0),
            "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.3, log=True),
            "max_depth": trial.suggest_int("max_depth", 3, 10),
            "min_child_weight": trial.suggest_int("min_child_weight", 1, 10),
            "n_estimators": trial.suggest_int("n_estimators", 100, 1000, step=50),
        }
    raise ValueError(model_name)


def prepare_inner_splits(X: pd.DataFrame, y: pd.Series, seed: int):
    prepared = []
    cv = KFold(n_splits=INNER_FOLDS, shuffle=True, random_state=seed)
    for train_idx, valid_idx in cv.split(X):
        X_inner_train = X.iloc[train_idx]
        X_inner_valid = X.iloc[valid_idx]
        prep = make_preprocessor(X_inner_train)
        Xt = prep.fit_transform(X_inner_train)
        Xv = prep.transform(X_inner_valid)
        prepared.append((Xt, y.iloc[train_idx].to_numpy(), Xv, y.iloc[valid_idx].to_numpy()))
    return prepared


def tune_model(model_name: str, X: pd.DataFrame, y: pd.Series, seed: int):
    prepared = prepare_inner_splits(X.reset_index(drop=True), y.reset_index(drop=True), seed)

    def objective(trial: optuna.Trial) -> float:
        params = suggest_params(model_name, trial)
        fold_rmse = []
        for inner_id, (Xt, yt, Xv, yv) in enumerate(prepared):
            model = make_model(model_name, params, seed + inner_id)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                model.fit(Xt, yt)
            fold_rmse.append(float(np.sqrt(mean_squared_error(yv, model.predict(Xv)))))
        return float(np.mean(fold_rmse))

    study = optuna.create_study(direction="minimize", sampler=TPESampler(seed=seed))
    study.optimize(objective, n_trials=TRIALS[model_name], show_progress_bar=False)
    return study.best_params, float(study.best_value)


def calculate_metrics(y_true, y_pred) -> dict:
    return {
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "r2": float(r2_score(y_true, y_pred)),
        "spearman": float(stats.spearmanr(y_true, y_pred).statistic),
        "pearson": float(stats.pearsonr(y_true, y_pred).statistic),
    }


def load_inputs():
    features = pd.read_csv(WT_DIR / "Xtbl_wt25.csv")
    labels = pd.read_csv(WT_DIR / "Ytbl_wt25.csv")
    train_indices = pd.read_csv(WT_DIR / "train_inds_tbl.csv") - 1
    test_indices = pd.read_csv(WT_DIR / "test_inds_tbl.csv") - 1
    X_full = features.drop(columns=["transcript", "random_factor", "random_num"])
    X_baseline = X_full.drop(columns=NOVEL_FEATURES)
    y = labels["RE"].copy()
    assert len(features) == len(labels) == 829
    assert set(NOVEL_FEATURES).issubset(X_full.columns)
    return features, labels, X_full, X_baseline, y, train_indices, test_indices


def run_repeated_outer_analysis(features, X_sets, y, train_indices, test_indices):
    metric_rows, parameter_rows, prediction_rows, importance_rows = [], [], [], []
    for feature_set, X in X_sets.items():
        categorical = [c for c in X.columns if X[c].dtype == "object" or isinstance(X[c].dtype, pd.CategoricalDtype)]
        for split in range(10):
            train_idx = train_indices.iloc[:, split].to_numpy(dtype=int)
            test_idx = test_indices.iloc[:, split].to_numpy(dtype=int)
            X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
            y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
            prep = make_preprocessor(X_train)
            Xt = prep.fit_transform(X_train)
            Xv = prep.transform(X_test)
            encoded_names = list(prep.get_feature_names_out())
            for model_offset, model_name in enumerate(["ridge", "lasso", "rf", "xgb"]):
                seed = SEED + 10000 * (feature_set == "full") + 100 * split + model_offset
                best_params, inner_rmse = tune_model(model_name, X_train, y_train, seed)
                model = make_model(model_name, best_params, seed)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    model.fit(Xt, y_train)
                pred = model.predict(Xv)
                row = {"feature_set": feature_set, "outer_split": split + 1, "model": model_name,
                       "inner_cv_rmse": inner_rmse, **calculate_metrics(y_test, pred)}
                metric_rows.append(row)
                parameter_rows.append({"feature_set": feature_set, "outer_split": split + 1,
                                       "model": model_name, "inner_cv_rmse": inner_rmse,
                                       "best_params_json": json.dumps(best_params, sort_keys=True)})
                prediction_rows.extend(
                    {"feature_set": feature_set, "outer_split": split + 1, "model": model_name,
                     "row_index": int(idx), "transcript": features.iloc[idx]["transcript"],
                     "y_true": float(y.iloc[idx]), "y_pred": float(p)}
                    for idx, p in zip(test_idx, pred)
                )
                values = getattr(model, "feature_importances_", None)
                if values is None and hasattr(model, "coef_"):
                    values = np.abs(np.ravel(model.coef_))
                if values is not None:
                    for encoded, value in zip(encoded_names, values):
                        raw = raw_feature_name(encoded, categorical)
                        importance_rows.append({
                            "feature_set": feature_set, "outer_split": split + 1, "model": model_name,
                            "encoded_feature": clean_feature_name(encoded), "raw_feature": raw,
                            "is_novel": raw in NOVEL_FEATURES, "importance": float(value),
                        })
                print(f"completed {feature_set} split {split + 1}/10 {model_name}: RMSE={row['rmse']:.4f}")
            pd.DataFrame(metric_rows).to_csv(OUT_DIR / "fold_metrics.partial.csv", index=False)
            pd.DataFrame(parameter_rows).to_csv(OUT_DIR / "best_parameters.partial.csv", index=False)

    metrics = pd.DataFrame(metric_rows)
    parameters = pd.DataFrame(parameter_rows)
    predictions = pd.DataFrame(prediction_rows)
    importances = pd.DataFrame(importance_rows)
    metrics.to_csv(OUT_DIR / "fold_metrics.csv", index=False)
    parameters.to_csv(OUT_DIR / "best_parameters.csv", index=False)
    predictions.to_csv(OUT_DIR / "outer_test_predictions.csv", index=False)
    importances.to_csv(OUT_DIR / "feature_importance_by_fold.csv", index=False)
    return metrics, parameters, predictions, importances


def summarize_outer_results(metrics: pd.DataFrame, importances: pd.DataFrame):
    summary = metrics.groupby(["feature_set", "model"])[METRICS].agg(["mean", "std"])
    summary.columns = [f"{metric}_{stat}" for metric, stat in summary.columns]
    summary.reset_index().to_csv(OUT_DIR / "summary_metrics.csv", index=False)

    paired_rows = []
    for model in metrics["model"].unique():
        for metric in METRICS:
            pivot = metrics[metrics.model == model].pivot(index="outer_split", columns="feature_set", values=metric)
            delta = pivot["full"] - pivot["baseline"]
            paired_rows.append({
                "model": model, "metric": metric, "mean_full_minus_baseline": delta.mean(),
                "sd_difference": delta.std(ddof=1), "full_better_splits": int((delta < 0).sum())
                if metric in ["rmse", "mae"] else int((delta > 0).sum()),
                "wilcoxon_p_two_sided": float(stats.wilcoxon(delta).pvalue),
            })
    pd.DataFrame(paired_rows).to_csv(OUT_DIR / "paired_feature_set_comparison.csv", index=False)

    full = metrics[metrics.feature_set == "full"].copy()
    full["rmse_rank"] = full.groupby("outer_split")["rmse"].rank(method="min")
    full["spearman_rank"] = full.groupby("outer_split")["spearman"].rank(method="min", ascending=False)
    full.to_csv(OUT_DIR / "full_model_fold_rankings.csv", index=False)

    xgb_imp = importances[(importances.feature_set == "full") & (importances.model == "xgb")]
    grouped_fold = xgb_imp.groupby(["outer_split", "raw_feature", "is_novel"], as_index=False).importance.sum()
    grouped = grouped_fold.groupby(["raw_feature", "is_novel"]).importance.agg(["mean", "std"]).reset_index()
    grouped["rank"] = grouped["mean"].rank(method="min", ascending=False).astype(int)
    grouped.sort_values("rank").to_csv(OUT_DIR / "xgb_grouped_feature_importance.csv", index=False)

    individual = xgb_imp.groupby(["encoded_feature", "raw_feature", "is_novel"]).importance.agg(["mean", "std"]).reset_index()
    individual["rank"] = individual["mean"].rank(method="min", ascending=False).astype(int)
    individual.sort_values("rank").to_csv(OUT_DIR / "xgb_individual_feature_importance.csv", index=False)


def write_hyperparameter_supplement(parameters: pd.DataFrame):
    """Write the search spaces and split-specific selections for Supplementary Table S3."""
    search_spaces = pd.DataFrame([
        {"model": "Ridge", "parameter": "alpha", "search_space": "0.01 to 100", "sampling": "log-uniform"},
        {"model": "Lasso", "parameter": "alpha", "search_space": "0.0001 to 1", "sampling": "log-uniform"},
        {"model": "Random Forest", "parameter": "n_estimators", "search_space": "100 to 500 in steps of 25", "sampling": "integer"},
        {"model": "Random Forest", "parameter": "max_depth", "search_space": "3 to 20", "sampling": "integer"},
        {"model": "Random Forest", "parameter": "min_samples_split", "search_space": "2 to 10", "sampling": "integer"},
        {"model": "Random Forest", "parameter": "min_samples_leaf", "search_space": "1 to 10", "sampling": "integer"},
        {"model": "Random Forest", "parameter": "max_features", "search_space": "sqrt, log2, or all features", "sampling": "categorical"},
        {"model": "XGBoost", "parameter": "reg_lambda", "search_space": "0 to 10", "sampling": "uniform"},
        {"model": "XGBoost", "parameter": "reg_alpha", "search_space": "0 to 10", "sampling": "uniform"},
        {"model": "XGBoost", "parameter": "colsample_bytree", "search_space": "0.5 to 1", "sampling": "uniform"},
        {"model": "XGBoost", "parameter": "gamma", "search_space": "0 to 5", "sampling": "uniform"},
        {"model": "XGBoost", "parameter": "subsample", "search_space": "0.5 to 1", "sampling": "uniform"},
        {"model": "XGBoost", "parameter": "learning_rate", "search_space": "0.01 to 0.3", "sampling": "log-uniform"},
        {"model": "XGBoost", "parameter": "max_depth", "search_space": "3 to 10", "sampling": "integer"},
        {"model": "XGBoost", "parameter": "min_child_weight", "search_space": "1 to 10", "sampling": "integer"},
        {"model": "XGBoost", "parameter": "n_estimators", "search_space": "100 to 1000 in steps of 50", "sampling": "integer"},
    ])
    search_spaces.insert(2, "optuna_trials_per_outer_split", search_spaces.model.map(
        {"Ridge": TRIALS["ridge"], "Lasso": TRIALS["lasso"],
         "Random Forest": TRIALS["rf"], "XGBoost": TRIALS["xgb"]}
    ))
    search_spaces.insert(3, "inner_validation", f"{INNER_FOLDS}-fold cross-validation")
    search_spaces.to_csv(OUT_DIR / "supplementary_table_s3_search_spaces.csv", index=False)

    model_labels = {"ridge": "Ridge", "lasso": "Lasso", "rf": "Random Forest", "xgb": "XGBoost"}
    selected_rows = []
    for row in parameters.itertuples(index=False):
        for parameter, value in json.loads(row.best_params_json).items():
            if parameter == "max_features" and value is None:
                value = "all features"
            selected_rows.append({
                "feature_set": row.feature_set,
                "outer_split": int(row.outer_split),
                "model": model_labels[row.model],
                "inner_cv_rmse": float(row.inner_cv_rmse),
                "parameter": parameter,
                "selected_value": value,
            })
    pd.DataFrame(selected_rows).sort_values(
        ["feature_set", "outer_split", "model", "parameter"]
    ).to_csv(OUT_DIR / "supplementary_table_s3_selected_parameters.csv", index=False)


def generate_figure3(metrics: pd.DataFrame, predictions: pd.DataFrame):
    """Generate the revised three-panel model-performance figure."""
    model_order = ["ridge", "lasso", "rf", "xgb"]
    feature_sets = ["baseline", "full"]
    colors = {"baseline": "#c9a46d", "full": "#8fa877"}
    labels = {"baseline": "Without New Features", "full": "With New Features"}

    fig = plt.figure(figsize=(11.5, 10.5))
    grid = fig.add_gridspec(2, 6, height_ratios=[1, 1.15], hspace=0.72, wspace=0.78)
    axes = [fig.add_subplot(grid[0, :3]), fig.add_subplot(grid[0, 3:])]
    x = np.arange(len(model_order))
    width = 0.36
    for ax, metric, ylabel, title in [
        (axes[0], "r2", r"R$^2$", r"(a) R$^2$ score"),
        (axes[1], "spearman", "Spearman", "(b) Spearman correlation"),
    ]:
        for offset, feature_set in enumerate(feature_sets):
            grouped = (
                metrics[metrics.feature_set == feature_set]
                .groupby("model")[metric]
                .agg(["mean", "std"])
                .reindex(model_order)
            )
            ax.bar(
                x + (offset - 0.5) * width,
                grouped["mean"],
                width,
                yerr=grouped["std"],
                color=colors[feature_set],
                edgecolor="none",
                capsize=3,
                label=labels[feature_set],
            )
        ax.set_xticks(x, [MODEL_LABELS[m] for m in model_order])
        ax.set_ylabel(ylabel, fontsize=17, fontweight="bold")
        ax.set_title(title, loc="center", fontweight="bold", fontsize=20)
        ax.tick_params(axis="both", labelsize=15)
        for label in ax.get_xticklabels():
            label.set_fontweight("bold")
        ax.grid(axis="y", alpha=0.2)
        ax.set_axisbelow(True)

    handles, legend_labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles,
        legend_labels,
        loc="upper center",
        bbox_to_anchor=(0.5, 0.535),
        ncol=2,
        frameon=False,
        fontsize=15,
    )

    ax = fig.add_subplot(grid[1, 1:5])
    scatter = predictions[(predictions.feature_set == "full") & (predictions.model == "xgb")]
    ax.scatter(
        scatter.y_true,
        scatter.y_pred,
        s=28,
        alpha=0.5,
        color="#1f77b4",
        edgecolor="black",
        linewidth=0.45,
    )
    slope, intercept = np.polyfit(scatter.y_true, scatter.y_pred, 1)
    x_line = np.linspace(scatter.y_true.min(), scatter.y_true.max(), 200)
    ax.plot(x_line, slope * x_line + intercept, color="red", linewidth=2.2)
    ax.set_xlabel("True\nReadthrough", fontsize=17, fontweight="bold")
    ax.set_ylabel("Predicted Readthrough", fontsize=17, fontweight="bold")
    ax.set_title("(c) XGBoost True vs. Predicted", loc="center", fontweight="bold", fontsize=20)
    ax.tick_params(axis="both", labelsize=15)
    ax.grid(alpha=0.45)
    ax.set_axisbelow(True)
    fig.subplots_adjust(top=0.95, bottom=0.08, left=0.09, right=0.98)
    fig.savefig(OUT_DIR / "figure3_model_performance.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "figure3_model_performance.pdf", bbox_inches="tight")
    plt.close(fig)


def collect_full_xgb_oof_shap(
    X_full: pd.DataFrame,
    y: pd.Series,
    train_indices: pd.DataFrame,
    test_indices: pd.DataFrame,
    parameters: pd.DataFrame,
):
    """Refit the saved full XGBoost models and collect held-out SHAP values."""
    encoded_values: dict[str, dict[str, list[np.ndarray]]] = {}
    grouped_values: dict[str, list[np.ndarray]] = {name: [] for name in X_full.columns}
    raw_values: dict[str, list[np.ndarray]] = {name: [] for name in X_full.columns}
    for split in range(10):
        train_idx = train_indices.iloc[:, split].to_numpy(dtype=int)
        test_idx = test_indices.iloc[:, split].to_numpy(dtype=int)
        X_train, X_test = X_full.iloc[train_idx], X_full.iloc[test_idx]
        y_train = y.iloc[train_idx]
        row = parameters[
            (parameters.feature_set == "full")
            & (parameters.model == "xgb")
            & (parameters.outer_split == split + 1)
        ]
        if len(row) != 1:
            raise ValueError(f"Expected one saved full-XGBoost parameter row for split {split + 1}")
        params = json.loads(row.iloc[0].best_params_json)
        seed = SEED + 10000 + 100 * split + 3
        prep = make_preprocessor(X_train)
        transformed_train = prep.fit_transform(X_train)
        transformed_test = prep.transform(X_test)
        model = make_model("xgb", params, seed)
        model.fit(transformed_train, y_train)
        contributions = model.get_booster().predict(
            xgb.DMatrix(transformed_test), pred_contribs=True
        )[:, :-1]
        encoded_names = [clean_feature_name(name) for name in prep.get_feature_names_out()]
        for index, name in enumerate(encoded_names):
            store = encoded_values.setdefault(name, {"shap": [], "value": []})
            store["shap"].append(contributions[:, index])
            store["value"].append(np.asarray(transformed_test)[:, index])
        categorical = [
            c for c in X_train.columns
            if X_train[c].dtype == "object" or isinstance(X_train[c].dtype, pd.CategoricalDtype)
        ]
        grouped = grouped_xgb_shap(
            model, transformed_test, prep, list(X_full.columns), categorical
        )
        for index, name in enumerate(X_full.columns):
            grouped_values[name].append(grouped[:, index])
            if name in categorical:
                # A continuous color scale is not meaningful for unordered categories.
                raw_values[name].append(np.full(len(X_test), np.nan))
            else:
                raw_values[name].append(pd.to_numeric(X_test[name], errors="coerce").to_numpy())

    encoded = {
        name: {
            "shap": np.concatenate(values["shap"]),
            "value": np.concatenate(values["value"]),
        }
        for name, values in encoded_values.items()
    }
    grouped = {name: np.concatenate(values) for name, values in grouped_values.items()}
    raw = {name: np.concatenate(values) for name, values in raw_values.items()}
    encoded_summary = pd.DataFrame(
        {
            "encoded_feature": name,
            "mean_absolute_shap": float(np.mean(np.abs(values["shap"]))),
            "mean_shap": float(np.mean(values["shap"])),
        }
        for name, values in encoded.items()
    ).sort_values("mean_absolute_shap", ascending=False)
    grouped_summary = pd.DataFrame(
        {
            "raw_feature": name,
            "display_feature": display_raw_feature(name),
            "is_novel": name in NOVEL_FEATURES,
            "mean_absolute_shap": float(np.mean(np.abs(values))),
            "mean_shap": float(np.mean(values)),
        }
        for name, values in grouped.items()
    ).sort_values("mean_absolute_shap", ascending=False)
    encoded_summary.to_csv(OUT_DIR / "xgb_oof_shap_encoded_summary.csv", index=False)
    grouped_summary.to_csv(OUT_DIR / "xgb_oof_shap_grouped_summary.csv", index=False)
    return encoded, encoded_summary, grouped, raw, grouped_summary


def draw_feature_importance(ax, importances: pd.DataFrame):
    data = (
        importances[(importances.feature_set == "full") & (importances.model == "xgb")]
        .groupby(["encoded_feature", "raw_feature"], as_index=False)
        .agg(mean=("importance", "mean"), std=("importance", "std"))
        .nlargest(50, "mean")
        .sort_values("mean")
    )
    data["group"] = data.raw_feature.map(feature_group)
    data["label"] = data.encoded_feature.map(display_encoded_feature)
    ax.barh(
        data.label,
        data["mean"],
        xerr=data["std"],
        color=data.group.map(FEATURE_PALETTE),
        edgecolor="none",
        error_kw={"ecolor": "#555555", "elinewidth": 0.7, "capsize": 1.5},
    )
    ax.set_xlabel("Mean feature importance", fontsize=16, fontweight="bold")
    ax.set_title(
        "(a) Top 50 feature Importances for XGB",
        loc="left",
        fontweight="bold",
        fontsize=18,
    )
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    ax.grid(axis="x", alpha=0.2)
    ax.set_axisbelow(True)
    handles = [Patch(facecolor=color, label=label) for label, color in FEATURE_PALETTE.items()]
    legend = ax.legend(handles=handles, loc="lower right", frameon=False, fontsize=13)
    for label in legend.get_texts():
        label.set_fontweight("bold")


def draw_shap_beeswarm(ax, encoded, encoded_summary, top_n=20):
    # Reuse the original SHAP summary-plot function and its conventional beeswarm
    # appearance, while supplying the corrected held-out SHAP values and clean labels.
    top_names = encoded_summary.head(top_n).encoded_feature.tolist()
    shap_matrix = np.column_stack([encoded[name]["shap"] for name in top_names])
    feature_matrix = np.column_stack([encoded[name]["value"] for name in top_names])
    sample_size = min(200, len(shap_matrix))
    sample = np.random.default_rng(SEED).choice(len(shap_matrix), sample_size, replace=False)
    plt.sca(ax)
    shap.summary_plot(
        shap_matrix[sample],
        features=feature_matrix[sample],
        feature_names=[display_encoded_feature(name) for name in top_names],
        max_display=top_n,
        show=False,
        plot_size=None,
    )
    ax.set_xlabel(
        "SHAP value (contribution to predicted RE)",
        fontsize=16,
        fontweight="bold",
    )
    ax.set_title(
        "(b) SHAP Summary Plot for XGB",
        loc="left",
        fontweight="bold",
        fontsize=18,
    )
    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)
    colorbar_axes = [candidate for candidate in ax.figure.axes if candidate is not ax]
    if colorbar_axes:
        colorbar_ax = colorbar_axes[-1]
        colorbar_ax.set_ylabel("Feature value", fontsize=16, fontweight="bold")
        colorbar_ax.tick_params(labelsize=14)
        for label in colorbar_ax.get_yticklabels():
            label.set_fontweight("bold")


def generate_figure4(importances, encoded, encoded_summary, grouped, raw, grouped_summary):
    fig, ax = plt.subplots(figsize=(10.5, 12.5))
    draw_feature_importance(ax, importances)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "figure4a_feature_importance.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "figure4a_feature_importance.pdf", bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10.5, 8.5))
    draw_shap_beeswarm(ax, encoded, encoded_summary)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "figure4b_xgb_shap.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "figure4b_xgb_shap.pdf", bbox_inches="tight")
    plt.close(fig)

    fig = plt.figure(figsize=(11.5, 21))
    grid = fig.add_gridspec(2, 1, height_ratios=[1.65, 1], hspace=0.18)
    ax_a = fig.add_subplot(grid[0])
    draw_feature_importance(ax_a, importances)
    ax_b = fig.add_subplot(grid[1])
    draw_shap_beeswarm(ax_b, encoded, encoded_summary)
    fig.subplots_adjust(top=0.985, bottom=0.035, left=0.29, right=0.95)
    fig.savefig(OUT_DIR / "figure4_feature_importance_and_shap.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "figure4_feature_importance_and_shap.pdf", bbox_inches="tight")
    plt.close(fig)

    # Use the same conventional SHAP summary-plot format as Figure 4b, expanded
    # to every original model feature. One-hot contributions are combined under
    # their original categorical feature; those unordered categories are gray.
    all_names = grouped_summary.raw_feature.tolist()
    shap_matrix = np.column_stack([grouped[name] for name in all_names])
    feature_matrix = np.column_stack([raw[name] for name in all_names])
    sample_size = min(200, len(shap_matrix))
    sample = np.random.default_rng(SEED).choice(len(shap_matrix), sample_size, replace=False)
    fig, ax = plt.subplots(figsize=(12, 26))
    plt.sca(ax)
    shap.summary_plot(
        shap_matrix[sample],
        features=feature_matrix[sample],
        feature_names=[display_raw_feature(name) for name in all_names],
        max_display=len(all_names),
        show=False,
        plot_size=None,
    )
    ax.set_title("SHAP Summary Plot for XGB: all features", loc="left", fontweight="bold", fontsize=20)
    ax.set_xlabel("SHAP value (contribution to predicted RE)", fontsize=18, fontweight="bold")
    ax.tick_params(axis="x", labelsize=15)
    ax.tick_params(axis="y", labelsize=16)
    colorbar_axes = [candidate for candidate in fig.axes if candidate is not ax]
    if colorbar_axes:
        colorbar_ax = colorbar_axes[-1]
        colorbar_ax.set_ylabel("Feature value", fontsize=18, fontweight="bold")
        colorbar_ax.tick_params(labelsize=15)
        for label in colorbar_ax.get_yticklabels():
            label.set_fontweight("bold")
    fig.savefig(OUT_DIR / "supplementary_all_feature_shap.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "supplementary_all_feature_shap.pdf", bbox_inches="tight")
    plt.close(fig)


def generate_figure6(random_gene_predictions: pd.DataFrame, random_summary: pd.DataFrame):
    """Regenerate Figure 6 from matched out-of-fold native and randomized predictions."""
    native = (
        random_gene_predictions
        .drop_duplicates("row_index")
        .sort_values("row_index")
    )
    randomized_per_gene = (
        random_gene_predictions.groupby("row_index", as_index=False)
        .randomized_oof_prediction.mean()
        .sort_values("row_index")
    )
    if len(native) != len(randomized_per_gene):
        raise ValueError("Figure 6 requires one native and one randomized mean per gene")

    measured = native.measured_RE.to_numpy()
    native_oof = native.native_oof_prediction.to_numpy()
    randomized_mean = randomized_per_gene.randomized_oof_prediction.to_numpy()
    all_values = np.concatenate([measured, native_oof, randomized_mean])
    bins = np.linspace(all_values.min(), all_values.max(), 50)
    colors = {"measured": "#1f77b4", "native": "#ff7f0e", "randomized": "#2ca02c"}

    # Retain the submitted Figure 6 layout and styling; only the underlying
    # matched out-of-fold predictions and displayed values are updated.
    fig, axes = plt.subplots(2, 1, figsize=(10, 10), gridspec_kw={"height_ratios": [1.35, 1]})
    ax = axes[0]
    for values, color, label in [
        (measured, colors["measured"], "Measured RE"),
        (native_oof, colors["native"], "OOF predicted RE"),
        (randomized_mean, colors["randomized"], "Randomized-set mean RE"),
    ]:
        ax.hist(values, bins=bins, density=True, alpha=0.5, color=color, label=label)
    ax.axvline(measured.mean(), color=colors["measured"], linestyle="--", linewidth=2,
               label=f"Measured mean = {measured.mean():.2f}")
    ax.axvline(native_oof.mean(), color=colors["native"], linestyle="--", linewidth=2,
               label=f"OOF mean = {native_oof.mean():.2f}")
    ax.axvline(randomized_mean.mean(), color=colors["randomized"], linestyle="--", linewidth=2,
               label=f"Random mean = {randomized_mean.mean():.2f}")
    ax.set_xlabel("RE value", fontsize=15, fontweight="bold")
    ax.set_ylabel("Density", fontsize=15, fontweight="bold")
    ax.set_title(
        "(a) Distribution of RE values: measured vs OOF vs randomized (mean)",
        fontweight="bold",
        fontsize=17,
    )
    ax.tick_params(axis="both", labelsize=12)
    ax.legend(fontsize=11, loc="upper right")
    ax.grid(alpha=0.35)
    ax.set_axisbelow(True)

    ax = axes[1]
    genome = random_summary.randomized_genome.to_numpy(dtype=int)
    randomized_genome_mean = random_summary.randomized_oof_mean.to_numpy()
    native_mean = float(random_summary.native_oof_mean.iloc[0])
    ax.scatter(genome, randomized_genome_mean, color=colors["measured"], s=42, zorder=3)
    ax.axhline(native_mean, color=colors["measured"], linestyle="--", linewidth=2)
    ax.set_xlabel("Randomized transcript-set index", fontsize=15, fontweight="bold")
    ax.set_ylabel("Mean predicted RE (XGB)", fontsize=15, fontweight="bold")
    ax.set_title(
        "(b) Mean predicted RE per randomized transcript set vs OOF baseline",
        fontweight="bold",
        fontsize=17,
    )
    ax.tick_params(axis="both", labelsize=12)
    ax.set_xticks(np.arange(1, 21))
    ax.set_xlim(0.3, 20.7)
    ax.grid(alpha=0.35)
    ax.set_axisbelow(True)

    fig.subplots_adjust(top=0.96, bottom=0.07, left=0.12, right=0.98, hspace=0.48)
    fig.savefig(OUT_DIR / "figure6_randomized_genomes.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "figure6_randomized_genomes.pdf", bbox_inches="tight")
    plt.close(fig)


def generate_final_figures(
    X_full,
    y,
    train_indices,
    test_indices,
    metrics,
    parameters,
    predictions,
    importances,
):
    """Generate every main and supplementary figure affected by model corrections."""
    generate_figure3(metrics, predictions)
    encoded, encoded_summary, grouped, raw, grouped_summary = collect_full_xgb_oof_shap(
        X_full, y, train_indices, test_indices, parameters
    )
    generate_figure4(importances, encoded, encoded_summary, grouped, raw, grouped_summary)


def run_randomized_genome_analysis(features, X_full, y):
    """Nested 5-fold XGB OOF predictions for native and matched randomized genes."""
    transcript_to_rand = {}
    for genome in range(1, 21):
        table = pd.read_csv(RAND_DIR / f"Xtbl_rand_{genome}.csv")
        # RNALFold's parser leaves MFE=Inf and dist_bp=-1 when it emits no
        # candidate local structure. Zero is the physically meaningful MFE for
        # an unstructured sequence and matches the native feature's allowed max.
        invalid_mfe = np.isinf(table["MFE"].to_numpy())
        if invalid_mfe.any():
            if not (table.loc[invalid_mfe, "dist_bp"] == -1).all():
                raise ValueError(f"Unexpected infinite MFE values in randomized genome {genome}")
            table.loc[invalid_mfe, "MFE"] = 0.0
        table = table.set_index("transcript", drop=False)
        assert set(features.transcript) == set(table.index)
        transcript_to_rand[genome] = table

    native_oof = np.full(len(features), np.nan)
    random_oof = {genome: np.full(len(features), np.nan) for genome in range(1, 21)}
    raw_columns = list(X_full.columns)
    native_shap = np.full((len(features), len(raw_columns)), np.nan)
    random_shap = {
        genome: np.full((len(features), len(raw_columns)), np.nan)
        for genome in range(1, 21)
    }
    params_rows = []
    outer = KFold(n_splits=5, shuffle=True, random_state=SEED)
    for fold, (train_idx, test_idx) in enumerate(outer.split(X_full), start=1):
        X_train, X_test = X_full.iloc[train_idx], X_full.iloc[test_idx]
        y_train = y.iloc[train_idx]
        seed = SEED + 50000 + fold
        params, inner_rmse = tune_model("xgb", X_train, y_train, seed)
        prep = make_preprocessor(X_train)
        Xt = prep.fit_transform(X_train)
        model = make_model("xgb", params, seed)
        model.fit(Xt, y_train)
        categorical_columns = [
            c for c in X_train.columns
            if X_train[c].dtype == "object" or isinstance(X_train[c].dtype, pd.CategoricalDtype)
        ]
        X_test_transformed = prep.transform(X_test)
        native_oof[test_idx] = model.predict(X_test_transformed)
        native_shap[test_idx, :] = grouped_xgb_shap(
            model, X_test_transformed, prep, raw_columns, categorical_columns
        )
        held_out_transcripts = features.iloc[test_idx].transcript.tolist()
        for genome, table in transcript_to_rand.items():
            rand_rows = table.loc[held_out_transcripts]
            rand_X = rand_rows[X_full.columns]
            rand_transformed = prep.transform(rand_X)
            random_oof[genome][test_idx] = model.predict(rand_transformed)
            random_shap[genome][test_idx, :] = grouped_xgb_shap(
                model, rand_transformed, prep, raw_columns, categorical_columns
            )
        params_rows.append({"fold": fold, "inner_cv_rmse": inner_rmse,
                            "best_params_json": json.dumps(params, sort_keys=True)})
        print(f"completed randomized-genome nested fold {fold}/5")

    assert not np.isnan(native_oof).any()
    assert all(not np.isnan(values).any() for values in random_oof.values())
    assert not np.isnan(native_shap).any()
    assert all(not np.isnan(values).any() for values in random_shap.values())
    native_df = pd.DataFrame({"row_index": np.arange(len(features)), "transcript": features.transcript,
                              "measured_RE": y, "native_oof_prediction": native_oof})
    native_df.to_csv(OUT_DIR / "randomized_native_oof_predictions.csv", index=False)
    pd.DataFrame(params_rows).to_csv(OUT_DIR / "randomized_xgb_parameters.csv", index=False)

    long_rows = []
    summary_rows = []
    native_mean = float(native_oof.mean())
    for genome, pred in random_oof.items():
        delta_native = pred - native_oof
        delta_measured = pred - y.to_numpy()
        summary_rows.append({
            "randomized_genome": genome, "native_oof_mean": native_mean,
            "measured_mean": float(y.mean()), "randomized_oof_mean": float(pred.mean()),
            "mean_delta_vs_native_oof": float(delta_native.mean()),
            "median_delta_vs_native_oof": float(np.median(delta_native)),
            "genes_above_native_oof": int((delta_native > 0).sum()),
            "paired_wilcoxon_p_vs_native_oof_two_sided": float(
                stats.wilcoxon(delta_native, alternative="two-sided").pvalue
            ),
            "mean_delta_vs_measured": float(delta_measured.mean()),
            "median_delta_vs_measured": float(np.median(delta_measured)),
            "genes_above_measured": int((delta_measured > 0).sum()),
            "paired_wilcoxon_p_vs_measured_two_sided": float(
                stats.wilcoxon(delta_measured, alternative="two-sided").pvalue
            ),
        })
        long_rows.extend(
            {"randomized_genome": genome, "row_index": i, "transcript": features.iloc[i].transcript,
             "measured_RE": float(y.iloc[i]),
             "native_oof_prediction": float(native_oof[i]), "randomized_oof_prediction": float(pred[i]),
             "delta_vs_native_oof": float(delta_native[i]),
             "delta_vs_measured": float(delta_measured[i])}
            for i in range(len(features))
        )
    random_summary = pd.DataFrame(summary_rows)
    random_summary["paired_wilcoxon_fdr_vs_native_oof"] = benjamini_hochberg(
        random_summary["paired_wilcoxon_p_vs_native_oof_two_sided"]
    )
    random_summary["paired_wilcoxon_fdr_vs_measured"] = benjamini_hochberg(
        random_summary["paired_wilcoxon_p_vs_measured_two_sided"]
    )
    random_summary.to_csv(OUT_DIR / "randomized_genome_summary.csv", index=False)
    random_summary.to_csv(OUT_DIR / "supplementary_table_randomized_paired_tests.csv", index=False)
    random_gene_predictions = pd.DataFrame(long_rows)
    random_gene_predictions.to_csv(OUT_DIR / "randomized_gene_predictions.csv", index=False)
    generate_figure6(random_gene_predictions, random_summary)

    fig, axes = plt.subplots(2, 1, figsize=(12, 9), sharex=True)
    for ax, column, ylabel, title in [
        (axes[0], "delta_vs_native_oof", "Delta RE", "Randomized - native OOF"),
        (axes[1], "delta_vs_measured", "Delta RE", "Randomized - measured RE"),
    ]:
        values = [
            random_gene_predictions.loc[
                random_gene_predictions.randomized_genome == genome, column
            ].to_numpy()
            for genome in range(1, 21)
        ]
        violins = ax.violinplot(
            values, positions=np.arange(1, 21), widths=0.8,
            showmeans=False, showmedians=True, showextrema=False,
        )
        for body in violins["bodies"]:
            body.set_facecolor("#56B4E9")
            body.set_edgecolor("#176A8A")
            body.set_alpha(0.65)
        violins["cmedians"].set_color("black")
        quartiles = np.array([np.percentile(value, [25, 75]) for value in values])
        ax.vlines(np.arange(1, 21), quartiles[:, 0], quartiles[:, 1], color="black", linewidth=3)
        ax.axhline(0, color="black", linewidth=1, linestyle="--")
        ax.set_ylabel(ylabel)
        ax.set_title(title)
    axes[1].set_xlabel("Randomized genome")
    axes[1].set_xticks(range(1, 21))
    fig.tight_layout()
    fig.savefig(OUT_DIR / "randomized_delta_re_violin.png", dpi=300)
    fig.savefig(OUT_DIR / "randomized_delta_re_violin.pdf")
    plt.close(fig)

    shap_rows = []
    for genome in range(1, 21):
        delta_shap = random_shap[genome] - native_shap
        for feature_index, feature in enumerate(raw_columns):
            values = delta_shap[:, feature_index]
            shap_rows.append({
                "randomized_genome": genome,
                "feature": feature,
                "mean_absolute_delta_shap": float(np.mean(np.abs(values))),
                "mean_delta_shap": float(np.mean(values)),
                "median_delta_shap": float(np.median(values)),
            })
    shap_by_genome = pd.DataFrame(shap_rows)
    shap_by_genome.to_csv(OUT_DIR / "randomized_shap_by_genome_feature.csv", index=False)
    shap_summary = shap_by_genome.groupby("feature", as_index=False).agg(
        mean_absolute_delta_shap=("mean_absolute_delta_shap", "mean"),
        sd_absolute_delta_shap=("mean_absolute_delta_shap", "std"),
        mean_delta_shap=("mean_delta_shap", "mean"),
        sd_delta_shap=("mean_delta_shap", "std"),
    )
    shap_summary["rank"] = shap_summary["mean_absolute_delta_shap"].rank(
        method="min", ascending=False
    ).astype(int)
    shap_summary = shap_summary.sort_values("rank")
    shap_summary.to_csv(OUT_DIR / "randomized_shap_feature_summary.csv", index=False)

    top_signed = shap_summary.nlargest(20, "mean_delta_shap").sort_values("mean_delta_shap")
    labels = top_signed.feature.map(display_raw_feature)
    colors = np.where(top_signed.mean_delta_shap >= 0, "#D55E00", "#0072B2")
    fig, ax = plt.subplots(figsize=(10, 7.7))
    ax.barh(
        labels, top_signed.mean_delta_shap,
        xerr=top_signed.sd_delta_shap, color=colors, edgecolor="none",
        error_kw={"ecolor": "#444444", "elinewidth": 0.8, "capsize": 2},
    )
    ax.axvline(0, color="#444444", linewidth=1)
    ax.set_xlabel("Mean paired SHAP change (randomized - native OOF)")
    ax.set_title(
        "Net feature contributions to higher RE predictions after randomization",
        fontweight="bold",
    )
    ax.grid(axis="x", alpha=0.18)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "randomized_shap_net_contributions.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUT_DIR / "randomized_shap_net_contributions.pdf", bbox_inches="tight")
    plt.close(fig)

    overall = {
        "native_oof_mean": native_mean,
        "measured_mean": float(y.mean()),
        "mean_randomized_prediction": float(random_summary.randomized_oof_mean.mean()),
        "mean_randomized_minus_native": float(random_summary.mean_delta_vs_native_oof.mean()),
        "mean_randomized_minus_measured": float(random_summary.mean_delta_vs_measured.mean()),
        "randomized_genomes_above_native": int((random_summary.mean_delta_vs_native_oof > 0).sum()),
        "randomized_genomes_above_measured": int((random_summary.mean_delta_vs_measured > 0).sum()),
        "number_randomized_genomes": 20,
    }
    (OUT_DIR / "randomized_overall_summary.json").write_text(json.dumps(overall, indent=2))


def write_run_metadata():
    versions = {
        "seed": SEED, "inner_folds": INNER_FOLDS, "trials": TRIALS,
        "numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__,
        "matplotlib": matplotlib.__version__,
        "sklearn": __import__("sklearn").__version__, "xgboost": xgb.__version__, "optuna": optuna.__version__,
        "wt_dir": str(WT_DIR), "rand_dir": str(RAND_DIR),
    }
    (OUT_DIR / "run_metadata.json").write_text(json.dumps(versions, indent=2))


def main():
    optuna.logging.set_verbosity(optuna.logging.WARNING)
    features, labels, X_full, X_baseline, y, train_indices, test_indices = load_inputs()
    print(f"Using {X_full.shape[1]} full and {X_baseline.shape[1]} baseline raw features")
    if os.environ.get("FIGURES_ONLY") == "1":
        metrics = pd.read_csv(OUT_DIR / "fold_metrics.csv")
        parameters = pd.read_csv(OUT_DIR / "best_parameters.csv")
        predictions = pd.read_csv(OUT_DIR / "outer_test_predictions.csv")
        importances = pd.read_csv(OUT_DIR / "feature_importance_by_fold.csv")
        generate_final_figures(
            X_full, y, train_indices, test_indices,
            metrics, parameters, predictions, importances,
        )
        write_hyperparameter_supplement(parameters)
        random_gene_predictions = pd.read_csv(OUT_DIR / "randomized_gene_predictions.csv")
        random_summary = pd.read_csv(OUT_DIR / "randomized_genome_summary.csv")
        generate_figure6(random_gene_predictions, random_summary)
        random_summary.to_csv(OUT_DIR / "supplementary_table_randomized_paired_tests.csv", index=False)
        print(f"Final figures regenerated from saved corrected results in {OUT_DIR}")
        return
    if os.environ.get("RANDOM_ONLY") != "1":
        metrics, parameters, predictions, importances = run_repeated_outer_analysis(
            features, {"baseline": X_baseline, "full": X_full}, y, train_indices, test_indices
        )
        summarize_outer_results(metrics, importances)
        write_hyperparameter_supplement(parameters)
        generate_final_figures(
            X_full, y, train_indices, test_indices,
            metrics, parameters, predictions, importances,
        )
    run_randomized_genome_analysis(features, X_full, y)
    write_run_metadata()
    print(f"All outputs saved to {OUT_DIR}")


if __name__ == "__main__":
    main()
