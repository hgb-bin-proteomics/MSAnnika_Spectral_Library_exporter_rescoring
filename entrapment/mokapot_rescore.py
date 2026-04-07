#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import logging
import pandas as pd
import mokapot

from pyXLMS.parser.util import get_bool_from_value

from typing import Optional

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal

DEBUG = True
PQVALUES = [
    "PG.Cscore",
    "PG.Pvalue",
    "PG.PValue (Run-Wise)",
    "PG.Qvalue",
    "PG.QValue (Run-Wise)",
    "EG.GlobalPrecursorQvalue",
    "EG.MaxChannelQvalue",
    "EG.MinChannelQvalue",
    "EG.Qvalue",
    "EG.InSourceFragmentationParentQvalue",
    "EG.AvgProfileQvalue",
    "EG.MaxProfileQvalue",
    "EG.MinProfileQvalue",
    "EG.PercentileQvalue",
    "FG.Qvalue",
]
POSSIBLE_RESCORING_FEATURES = [
    "EG.Cscore",
    "FG.Charge",
    "PP.MatchedIonsA",
    "PP.TotalIonsA",
    "PP.MatchedIonsB",
    "PP.TotalIonsB",
    "PP.RelativeMatchScoreA",
    "PP.RelativeMatchScoreB",
    "PP.PartialCscoreA",
    "PP.PartialCscoreB",
    "PP.CompositeRelativeMatchScore",
    "PP.CompositePartialCscore",
    "PP.SequenceCoverageNTermAlpha",
    "PP.SequenceCoverageNTermBeta",
    "PP.SequenceCoverageNTermFull",
    "PP.SequenceCoverageCTermAlpha",
    "PP.SequenceCoverageCTermBeta",
    "PP.SequenceCoverageCTermFull",
    "PP.SequenceCoverageAlpha",
    "PP.SequenceCoverageBeta",
    "PP.SequenceCoverageFull",
    "PP.UniScoreAlpha",
    "PP.UniScoreBeta",
    "PP.UniScoreFull",
    "PP.PepLenAlpha",
    "PP.PepLenBeta",
    "PP.NumberCrosslinkFragmentsAlpha",
    "PP.NumberCrosslinkFragmentsBeta",
    "PP.NumberCrosslinkFragmentsFull",
    "PP.NormalizedCrosslinkFragmentsAlpha",
    "PP.NormalizedCrosslinkFragmentsBeta",
    "PP.NormalizedCrosslinkFragmentsFull",
]
POSSIBLE_RESCORING_FEATURES_ALPHA = [
    "EG.Cscore",
    "PP.MatchedIonsA",
    "PP.TotalIonsA",
    "PP.RelativeMatchScoreA",
    "PP.PartialCscoreA",
    "PP.SequenceCoverageNTermAlpha",
    "PP.SequenceCoverageCTermAlpha",
    "PP.SequenceCoverageAlpha",
    "PP.UniScoreAlpha",
    "PP.PepLenAlpha",
    "PP.NumberCrosslinkFragmentsAlpha",
    "PP.NormalizedCrosslinkFragmentsAlpha",
]
POSSIBLE_RESCORING_FEATURES_BETA = [
    "EG.Cscore",
    "PP.MatchedIonsB",
    "PP.TotalIonsB",
    "PP.RelativeMatchScoreB",
    "PP.PartialCscoreB",
    "PP.SequenceCoverageNTermBeta",
    "PP.SequenceCoverageCTermBeta",
    "PP.SequenceCoverageBeta",
    "PP.UniScoreBeta",
    "PP.PepLenBeta",
    "PP.NumberCrosslinkFragmentsBeta",
    "PP.NormalizedCrosslinkFragmentsBeta",
]

__version = "1.1.4"
logger = logging.getLogger(__name__)


def __do_mokapot_columns(
    orig_df: pd.DataFrame, use_p_and_q_values: bool, alpha: Optional[bool]
) -> pd.DataFrame:
    df = orig_df.copy()
    df["MP.Spectrum"] = df.apply(
        lambda row: (
            f"{row['R.FileName']}_{row['R.Condition']}_{row['PP.PseudoScanNumber']}"
        ),
        axis=1,
    )
    if alpha is None:
        df["MP.Target"] = df.apply(
            lambda row: (
                not (
                    get_bool_from_value(row["PP.IsDecoyA"])
                    or get_bool_from_value(row["PP.IsDecoyB"])
                )
            ),
            axis=1,
        )
        df["MP.Peptide"] = df.apply(
            lambda row: row["PP.PeptideA"] + row["PP.PeptideB"], axis=1
        )
        df["MP.Protein"] = df.apply(
            lambda row: ";".join(
                set(row["PP.ProteinA"].split(";")) | set(row["PP.ProteinB"].split(";"))
            ),
            axis=1,
        )
    elif alpha:
        df["MP.Target"] = df.apply(
            lambda row: not get_bool_from_value(row["PP.IsDecoyA"]), axis=1
        )
        df["MP.Peptide"] = df.apply(lambda row: row["PP.PeptideA"], axis=1)
        df["MP.Protein"] = df.apply(lambda row: row["PP.ProteinA"], axis=1)
    else:
        df["MP.Target"] = df.apply(
            lambda row: not get_bool_from_value(row["PP.IsDecoyB"]), axis=1
        )
        df["MP.Peptide"] = df.apply(lambda row: row["PP.PeptideB"], axis=1)
        df["MP.Protein"] = df.apply(lambda row: row["PP.ProteinB"], axis=1)
    possible_rescoring_features = POSSIBLE_RESCORING_FEATURES.copy()
    if use_p_and_q_values:
        possible_rescoring_features += PQVALUES
    dropped_features = [c for c in possible_rescoring_features if df[c].isna().any()]  # pyright: ignore[reportGeneralTypeIssues]
    if len(dropped_features) > 0:
        df.drop(columns=dropped_features, inplace=True)
    logger.info(
        f"Removed the following features because of missing values: {', '.join(dropped_features)}"
    )
    return df


def __rescore_csms(
    orig_df: pd.DataFrame, use_p_and_q_values: bool = False
) -> pd.DataFrame:
    df = __do_mokapot_columns(orig_df, use_p_and_q_values, None)
    possible_rescoring_features = POSSIBLE_RESCORING_FEATURES.copy()
    if use_p_and_q_values:
        possible_rescoring_features += PQVALUES
    rescoring_features = [
        feature for feature in possible_rescoring_features if feature in df
    ]
    if DEBUG:
        df.to_csv("pre.csv", index=False)
    psms = mokapot.dataset.LinearPsmDataset(
        psms=df,
        target_column="MP.Target",
        spectrum_columns="MP.Spectrum",  # pyright: ignore[reportArgumentType]
        peptide_column="MP.Peptide",
        protein_column="MP.Protein",
        feature_columns=rescoring_features,
        copy_data=True,
        rng=1337,
    )
    logger.info(psms)
    models, scores = mokapot.brew([psms], rng=1337)
    df["Mokapot Score"] = scores[0]
    if DEBUG:
        df.to_csv("post.csv", index=False)
    conf = mokapot.confidence.assign_confidence([psms], scores)
    psms_conf: pd.DataFrame = conf[0].psms
    logger.info(
        f"Target hits with q-value < 0.01: {psms_conf[psms_conf['mokapot_qvalue'] < 0.01].shape[0]}"
    )
    return df


def __rescore_psms_separately(orig_df: pd.DataFrame) -> pd.DataFrame:
    df = orig_df.copy()
    df_alpha = __do_mokapot_columns(df, use_p_and_q_values=False, alpha=True)
    df_beta = __do_mokapot_columns(df, use_p_and_q_values=False, alpha=False)
    rescoring_features_alpha = [
        feature for feature in POSSIBLE_RESCORING_FEATURES_ALPHA if feature in df
    ]
    rescoring_features_beta = [
        feature for feature in POSSIBLE_RESCORING_FEATURES_BETA if feature in df
    ]
    if DEBUG:
        df.to_csv("pre.csv", index=False)
        df_alpha.to_csv("pre_alpha.csv", index=False)
        df_beta.to_csv("pre_beta.csv", index=False)
    psms_alpha = mokapot.dataset.LinearPsmDataset(
        psms=df_alpha,
        target_column="MP.Target",
        spectrum_columns="MP.Spectrum",  # pyright: ignore[reportArgumentType]
        peptide_column="MP.Peptide",
        protein_column="MP.Protein",
        feature_columns=rescoring_features_alpha,
        copy_data=True,
        rng=1337,
    )
    logger.info(psms_alpha)
    psms_beta = mokapot.dataset.LinearPsmDataset(
        psms=df_beta,
        target_column="MP.Target",
        spectrum_columns="MP.Spectrum",  # pyright: ignore[reportArgumentType]
        peptide_column="MP.Peptide",
        protein_column="MP.Protein",
        feature_columns=rescoring_features_beta,
        copy_data=True,
        rng=1337,
    )
    logger.info(psms_beta)
    models_alpha, scores_alpha = mokapot.brew([psms_alpha], rng=1337)
    models_beta, scores_beta = mokapot.brew([psms_beta], rng=1337)
    df["Mokapot Score Alpha"] = scores_alpha[0]
    df["Mokapot Score Beta"] = scores_beta[0]
    df["Mokapot Score"] = df.apply(
        lambda row: min(row["Mokapot Score Alpha"], row["Mokapot Score Beta"]), axis=1
    )
    if DEBUG:
        df.to_csv("post.csv", index=False)
    return df


def __rescore_psms_merged(orig_df: pd.DataFrame) -> pd.DataFrame:
    df = orig_df.copy()
    df_alpha = __do_mokapot_columns(df, use_p_and_q_values=False, alpha=True)
    df_beta = __do_mokapot_columns(df, use_p_and_q_values=False, alpha=False)
    df_alpha["MERGE.CSMID"] = range(df_alpha.shape[0])
    df_beta["MERGE.CSMID"] = range(df_beta.shape[0])
    df_alpha.rename(
        columns={
            "PP.MatchedIonsA": "PP.MatchedIons",
            "PP.TotalIonsA": "PP.TotalIons",
            "PP.RelativeMatchScoreA": "PP.RelativeMatchScore",
            "PP.PartialCscoreA": "PP.PartialCscore",
            "PP.SequenceCoverageNTermAlpha": "PP.SequenceCoverageNTerm",
            "PP.SequenceCoverageCTermAlpha": "PP.SequenceCoverageCTerm",
            "PP.SequenceCoverageAlpha": "PP.SequenceCoverage",
            "PP.UniScoreAlpha": "PP.UniScore",
            "PP.PepLenAlpha": "PP.PepLen",
            "PP.NumberCrosslinkFragmentsAlpha": "PP.NumberCrosslinkFragments",
            "PP.NormalizedCrosslinkFragmentsAlpha": "PP.NormalizedCrosslinkFragments",
        },
        inplace=True,
    )
    df_beta.rename(
        columns={
            "PP.MatchedIonsB": "PP.MatchedIons",
            "PP.TotalIonsB": "PP.TotalIons",
            "PP.RelativeMatchScoreB": "PP.RelativeMatchScore",
            "PP.PartialCscoreB": "PP.PartialCscore",
            "PP.SequenceCoverageNTermBeta": "PP.SequenceCoverageNTerm",
            "PP.SequenceCoverageCTermBeta": "PP.SequenceCoverageCTerm",
            "PP.SequenceCoverageBeta": "PP.SequenceCoverage",
            "PP.UniScoreBeta": "PP.UniScore",
            "PP.PepLenBeta": "PP.PepLen",
            "PP.NumberCrosslinkFragmentsBeta": "PP.NumberCrosslinkFragments",
            "PP.NormalizedCrosslinkFragmentsBeta": "PP.NormalizedCrosslinkFragments",
        },
        inplace=True,
    )
    df_alpha.drop(
        columns=list(
            set(df_alpha.columns.tolist())
            - {
                "EG.Cscore",
                "PP.MatchedIons",
                "PP.TotalIons",
                "PP.RelativeMatchScore",
                "PP.PartialCscore",
                "PP.SequenceCoverageNTerm",
                "PP.SequenceCoverageCTerm",
                "PP.SequenceCoverage",
                "PP.UniScore",
                "PP.PepLen",
                "PP.NumberCrosslinkFragments",
                "PP.NormalizedCrosslinkFragments",
                "MP.Target",
                "MP.Spectrum",
                "MP.Peptide",
                "MP.Protein",
                "MERGE.CSMID",
            }
        ),
        inplace=True,
    )
    df_beta.drop(
        columns=list(
            set(df_beta.columns.tolist())
            - {
                "EG.Cscore",
                "PP.MatchedIons",
                "PP.TotalIons",
                "PP.RelativeMatchScore",
                "PP.PartialCscore",
                "PP.SequenceCoverageNTerm",
                "PP.SequenceCoverageCTerm",
                "PP.SequenceCoverage",
                "PP.UniScore",
                "PP.PepLen",
                "PP.NumberCrosslinkFragments",
                "PP.NormalizedCrosslinkFragments",
                "MP.Target",
                "MP.Spectrum",
                "MP.Peptide",
                "MP.Protein",
                "MERGE.CSMID",
            }
        ),
        inplace=True,
    )
    psms_df = pd.concat([df_alpha, df_beta], ignore_index=True)
    psms_df["MERGE.PSMID"] = range(psms_df.shape[0])
    psms_df["MP.Spectrum"] = psms_df.apply(
        lambda row: f"{row['MP.Spectrum']}_{row['MERGE.PSMID']}", axis=1
    )
    dropped_features = [
        c
        for c in [
            "EG.Cscore",
            "PP.MatchedIons",
            "PP.TotalIons",
            "PP.RelativeMatchScore",
            "PP.PartialCscore",
            "PP.SequenceCoverageNTerm",
            "PP.SequenceCoverageCTerm",
            "PP.SequenceCoverage",
            "PP.UniScore",
            "PP.PepLen",
            "PP.NumberCrosslinkFragments",
            "PP.NormalizedCrosslinkFragments",
        ]
        if psms_df[c].isna().any()  # pyright: ignore[reportGeneralTypeIssues]
    ]
    if len(dropped_features) > 0:
        psms_df.drop(columns=dropped_features, inplace=True)
    logger.info(
        f"Removed the following features because of missing values: {', '.join(dropped_features)}"
    )
    rescoring_features = [
        feature
        for feature in [
            "EG.Cscore",
            "PP.MatchedIons",
            "PP.TotalIons",
            "PP.RelativeMatchScore",
            "PP.PartialCscore",
            "PP.SequenceCoverageNTerm",
            "PP.SequenceCoverageCTerm",
            "PP.SequenceCoverage",
            "PP.UniScore",
            "PP.PepLen",
            "PP.NumberCrosslinkFragments",
            "PP.NormalizedCrosslinkFragments",
        ]
        if feature in psms_df
    ]
    if DEBUG:
        psms_df.to_csv("pre.csv", index=False)
    psms = mokapot.dataset.LinearPsmDataset(
        psms=psms_df,
        target_column="MP.Target",
        spectrum_columns="MP.Spectrum",  # pyright: ignore[reportArgumentType]
        peptide_column="MP.Peptide",
        protein_column="MP.Protein",
        feature_columns=rescoring_features,
        copy_data=True,
        rng=1337,
    )
    logger.info(psms)
    models, scores = mokapot.brew([psms], rng=1337)
    scores_a = list()
    scores_b = list()
    scores_csm = list()
    psms_df["Mokapot Score"] = scores[0]
    if DEBUG:
        psms_df.to_csv("post.csv", index=False)
    for i in range(df.shape[0]):
        alpha: pd.Series = psms_df.iloc[i]
        beta: pd.Series = psms_df.iloc[df.shape[0] + i]
        assert int(alpha["MERGE.CSMID"]) == int(beta["MERGE.CSMID"])
        scores_a.append(alpha["Mokapot Score"])
        scores_b.append(beta["Mokapot Score"])
        scores_csm.append(min([alpha["Mokapot Score"], beta["Mokapot Score"]]))
    df["Mokapot Score Alpha"] = scores_a
    df["Mokapot Score Beta"] = scores_b
    df["Mokapot Score"] = scores_csm
    return df


def rescore(
    data: str | pd.DataFrame,
    level: Literal["PSM", "CSM"] = "PSM",
    use_p_and_q_values: bool = False,
    **kwargs,
) -> pd.DataFrame:
    r"""Rescore an annotated Spectronaut result file with Mokapot.

    Reads an annotated Spectronaut result file with grouped crosslink-spectrum-matches and rescores
    them using Mokapot.

    Parameters
    ----------
    data : str, or pandas.DataFrame
        The name/path of the Spectronaut result file or a pandas DataFrame.
    level : "PSM", or "CSM", default = "PSM"
        Whether the whole CSMs (crosslink-spectrum-matches) should be rescored or
        the individual PSMs (peptide-spectrum-matches).
    use_p_and_q_values : bool, default = False
        Whether to use p- and q-values from Spectronaut for rescoring. Only applies
        to level "CSM".

    Returns
    -------
    pandas.DataFrame
        Returns a copy of the original result file including a new column named ``Mokapot Score``.

    Examples
    --------
    >>> from mokapot_rescore import rescore
    >>> df = rescore(
    ...     "data/THIDDIAXL003_DIAmethodEval_SN20c4_Report_FM_crosslinking_plusDecoy_req_DIA12_CV48.csv_annotated.csv_grouped_by_residue_pair.csv",
    ...     level="PSM",
    ... )
    >>> df["Mokapot Score"]
    0       -0.679182
    1       -1.037381
    2       -0.197376
    3       -0.796810
    4        0.318801
               ...
    14475   -2.460086
    14476   -0.099634
    14477   -0.300225
    14478   -0.293951
    14479   -0.325833
    Name: Mokapot Score, Length: 14480, dtype: float64
    """
    if isinstance(data, str):
        df = pd.read_csv(data, low_memory=False, **kwargs)
    else:
        df = data
    if level == "CSM":
        return __rescore_csms(df, use_p_and_q_values)
    return __rescore_psms_merged(df)
