#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

import pandas as pd
from tqdm import tqdm

from pyXLMS.data import create_csm
from pyXLMS.data import create_parser_result
from pyXLMS.parser.util import format_sequence
from pyXLMS.parser.util import get_bool_from_value
from pyXLMS.parser.util import __serialize_pandas_series

from typing import Dict
from typing import Any

try:
    from typing import Literal
except ImportError:
    from typing_extensions import Literal


def read(
    data: str | pd.DataFrame,
    score: Literal[
        "EG.Cscore",
        "PP.CompositeRelativeMatchScore",
        "PP.CompositePartialCscore",
        "PP.UniScoreFull",
        "Mokapot Score",
    ],
    **kwargs,
) -> Dict[str, Any]:
    if isinstance(data, str):
        df = pd.read_csv(data, low_memory=False, **kwargs)
    else:
        df = data
    csms = list()
    for i, row in tqdm(df.iterrows(), total=df.shape[0], desc="Reading CSMs..."):
        if score == "PP.CompositeRelativeMatchScore":
            score_a = float(row["PP.RelativeMatchScoreA"])
            score_b = float(row["PP.RelativeMatchScoreB"])
            score_csm = float(row["PP.CompositeRelativeMatchScore"])
        elif score == "PP.CompositePartialCscore":
            score_a = float(row["PP.PartialCscoreA"])
            score_b = float(row["PP.PartialCscoreB"])
            score_csm = float(row["PP.CompositePartialCscore"])
        elif score == "PP.UniScoreFull":
            score_a = float(row["PP.UniScoreAlpha"])
            score_b = float(row["PP.UniScoreBeta"])
            score_csm = float(row["PP.UniScoreFull"])
        elif score == "Mokapot Score":
            score_a = float(row["Mokapot Score Alpha"])
            score_b = float(row["Mokapot Score Beta"])
            score_csm = float(row["Mokapot Score"])
        else:
            score_a = float(row["EG.Cscore"])
            score_b = float(row["EG.Cscore"])
            score_csm = float(row["EG.Cscore"])
        csm = create_csm(
            peptide_a=format_sequence(str(row["PP.PeptideA"])),
            modifications_a=None,
            xl_position_peptide_a=int(row["PP.CrosslinkPositionPeptideA"]),
            proteins_a=str(row["PP.ProteinA"]).split(";"),
            xl_position_proteins_a=[
                int(pos) for pos in str(row["PP.CrosslinkPositionProteinA"]).split(";")
            ],
            pep_position_proteins_a=[
                int(pos) for pos in str(row["PP.PeptidePositionProteinA"]).split(";")
            ],
            score_a=score_a,
            decoy_a=get_bool_from_value(row["PP.IsDecoyA"]),
            peptide_b=format_sequence(str(row["PP.PeptideB"])),
            modifications_b=None,
            xl_position_peptide_b=int(row["PP.CrosslinkPositionPeptideB"]),
            proteins_b=str(row["PP.ProteinB"]).split(";"),
            xl_position_proteins_b=[
                int(pos) for pos in str(row["PP.CrosslinkPositionProteinB"]).split(";")
            ],
            pep_position_proteins_b=[
                int(pos) for pos in str(row["PP.PeptidePositionProteinB"]).split(";")
            ],
            score_b=score_b,
            decoy_b=get_bool_from_value(row["PP.IsDecoyB"]),
            score=score_csm,
            spectrum_file=f"{str(row['R.Condition']).strip()}:{str(row['R.FileName']).strip()}",
            scan_nr=int(row["PP.PseudoScanNumber"]),
            charge=int(row["FG.Charge"]),
            rt=None,
            im_cv=None,
            additional_information={"source": __serialize_pandas_series(row)},
        )
        csms.append(csm)
    if len(csms) == 0:
        raise RuntimeError(
            "No crosslink-spectrum-matches were parsed! If this is unexpected, please file a bug report!"
        )
    return create_parser_result(
        search_engine="Spectronaut",
        csms=csms,
        crosslinks=None,
    )
