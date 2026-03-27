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

POSSIBLE_RESCORING_FEATURES = [
    "PG.Cscore",
    "EG.Cscore",
    "FG.Charge",
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

__version = "1.0.1"
logger = logging.getLogger(__name__)


def __do_mokapot_columns(
    orig_df: pd.DataFrame, alpha: Optional[bool] = None
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
    df.dropna(axis=1, subset=POSSIBLE_RESCORING_FEATURES, inplace=True)
    return df


def __rescore_csms(orig_df: pd.DataFrame) -> pd.DataFrame:
    df = __do_mokapot_columns(orig_df)
    rescoring_features = [
        feature for feature in POSSIBLE_RESCORING_FEATURES if feature in df
    ]
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
    models, scores = mokapot.brew([psms])
    df["Mokapot Score"] = scores[0]
    conf = mokapot.confidence.assign_confidence([psms], scores)
    psms_conf: pd.DataFrame = conf[0].psms
    logger.info(
        f"Target hits with q-value < 0.01: {psms_conf[psms_conf['mokapot_qvalue'] < 0.01].shape[0]}"
    )
    return df


def __resore_psms(orig_df: pd.DataFrame) -> pd.DataFrame:
    df = orig_df.copy()
    df_alpha = __do_mokapot_columns(df, alpha=True)
    df_beta = __do_mokapot_columns(df, alpha=False)
    rescoring_features_alpha = [
        feature for feature in POSSIBLE_RESCORING_FEATURES_ALPHA if feature in df
    ]
    rescoring_features_beta = [
        feature for feature in POSSIBLE_RESCORING_FEATURES_BETA if feature in df
    ]
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
    models_alpha, scores_alpha = mokapot.brew([psms_alpha])
    models_beta, scores_beta = mokapot.brew([psms_beta])
    df["Mokapot Score Alpha"] = scores_alpha[0]
    df["Mokapot Score Beta"] = scores_beta[0]
    df["Mokapot Score"] = df.apply(
        lambda row: min(row["Mokapot Score Alpha"], row["Mokapot Score Beta"]), axis=1
    )
    return df


def rescore(
    data: str | pd.DataFrame,
    level: Literal["PSM", "CSM"] = "PSM",
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
    0       -0.570398
    1       -1.001960
    2       -0.193422
    3       -0.849390
    4        0.201500
               ...
    14475   -2.197316
    14476   -0.022281
    14477   -0.333449
    14478   -0.159475
    14479   -0.494221
    Name: Mokapot Score, Length: 14480, dtype: float64
    """
    if isinstance(data, str):
        df = pd.read_csv(data, low_memory=False, **kwargs)
    else:
        df = data
    if level == "CSM":
        return __rescore_csms(df)
    return __resore_psms(df)
