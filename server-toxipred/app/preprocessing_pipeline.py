"""Reusable descriptor preprocessing pipeline components.

These transformers keep all fit-time statistics train-derived and are meant to be
used in sklearn pipelines to reduce leakage risk from manual preprocessing.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


class MedianImputer(BaseEstimator, TransformerMixin):
    """Impute missing values with train-set medians."""

    def fit(self, X, y=None):
        X_df = _ensure_dataframe(X)
        self.feature_names_in_ = X_df.columns.tolist()
        self.medians_ = X_df.median(numeric_only=True)
        return self

    def transform(self, X):
        X_df = _ensure_dataframe(X, columns=getattr(self, "feature_names_in_", None))
        return X_df.fillna(self.medians_)


@dataclass
class NearConstantFeatureRemover(BaseEstimator, TransformerMixin):
    """Drop features dominated by a single value above a threshold."""

    threshold: float = 0.80

    def fit(self, X, y=None):
        X_df = _ensure_dataframe(X)
        keep_cols = []
        for col in X_df.columns:
            value_fractions = X_df[col].value_counts(normalize=True, dropna=False)
            if value_fractions.empty:
                continue
            if value_fractions.iloc[0] < self.threshold:
                keep_cols.append(col)
        self.keep_cols_ = keep_cols
        self.feature_names_in_ = X_df.columns.tolist()
        return self

    def transform(self, X):
        X_df = _ensure_dataframe(X, columns=getattr(self, "feature_names_in_", None))
        return X_df.loc[:, self.keep_cols_]


@dataclass
class CorrelationFeatureRemover(BaseEstimator, TransformerMixin):
    """Drop highly correlated features, keeping first occurrences."""

    threshold: float = 0.70

    def fit(self, X, y=None):
        X_df = _ensure_dataframe(X)
        corr_matrix = X_df.corr().abs()
        upper = corr_matrix.where(np.triu(np.ones(corr_matrix.shape), k=1).astype(bool))
        drop_cols = [col for col in upper.columns if any(upper[col] > self.threshold)]
        self.keep_cols_ = [col for col in X_df.columns if col not in drop_cols]
        self.feature_names_in_ = X_df.columns.tolist()
        return self

    def transform(self, X):
        X_df = _ensure_dataframe(X, columns=getattr(self, "feature_names_in_", None))
        return X_df.loc[:, self.keep_cols_]


@dataclass
class IQRClipper(BaseEstimator, TransformerMixin):
    """Clip each column based on train-set IQR limits."""

    factor: float = 1.5

    def fit(self, X, y=None):
        X_df = _ensure_dataframe(X)
        self.feature_names_in_ = X_df.columns.tolist()
        self.limits_ = {}
        for col in X_df.columns:
            q1 = X_df[col].quantile(0.25)
            q3 = X_df[col].quantile(0.75)
            iqr = q3 - q1
            if iqr == 0 or np.isnan(iqr):
                continue
            lower = q1 - self.factor * iqr
            upper = q3 + self.factor * iqr
            self.limits_[col] = (lower, upper)
        return self

    def transform(self, X):
        X_df = _ensure_dataframe(
            X, columns=getattr(self, "feature_names_in_", None)
        ).copy()
        for col, (lower, upper) in self.limits_.items():
            if col in X_df.columns:
                X_df[col] = X_df[col].clip(lower, upper)
        return X_df


class DataFrameStandardScaler(BaseEstimator, TransformerMixin):
    """StandardScaler wrapper that returns a DataFrame with original labels."""

    def __init__(self):
        self.scaler = StandardScaler()

    def fit(self, X, y=None):
        X_df = _ensure_dataframe(X)
        self.feature_names_in_ = X_df.columns.tolist()
        self.scaler.fit(X_df)
        return self

    def transform(self, X):
        X_df = _ensure_dataframe(X, columns=getattr(self, "feature_names_in_", None))
        arr = self.scaler.transform(X_df)
        return pd.DataFrame(arr, columns=X_df.columns, index=X_df.index)


def create_descriptor_preprocessor(
    similarity_threshold: float = 0.80,
    correlation_threshold: float = 0.70,
    iqr_factor: float = 1.5,
) -> Pipeline:
    """Create the descriptor preprocessing pipeline with train-fitted state."""
    return Pipeline(
        steps=[
            ("imputer", MedianImputer()),
            (
                "near_constant",
                NearConstantFeatureRemover(threshold=similarity_threshold),
            ),
            ("correlation", CorrelationFeatureRemover(threshold=correlation_threshold)),
            ("iqr_clip", IQRClipper(factor=iqr_factor)),
            ("scaler", DataFrameStandardScaler()),
        ]
    )


def _ensure_dataframe(X, columns=None) -> pd.DataFrame:
    if isinstance(X, pd.DataFrame):
        return X
    return pd.DataFrame(X, columns=columns)
