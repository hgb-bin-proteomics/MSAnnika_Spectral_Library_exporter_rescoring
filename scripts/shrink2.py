#!/usr/bin/env python3
#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "polars",
# ]
# ///

# Essential Data Extractor
# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

# REQUIREMENTS
# pip install polars
import polars as pl

# version tracking
__version = "1.0.2"
__date = "2026-03-13"

INPUT = "THIDDIAXL003_DIAmethodEval_SN20c4_Report_FM_crosslinking_plusDecoy_allCol.tsv"
OUTPUT = "THIDDIAXL003_DIAmethodEval_SN20c4_Report_FM_crosslinking_plusDecoy_req_DIA12_CV48.csv"
CONDITIONS = ["DIA12_CV48"]
COLS = [
    "R.FileName",
    "R.Condition",
    "PG.ProteinNames",
    "PG.Cscore",
    "EG.Library",
    "EG.PrecursorId",
    "EG.Cscore",
    "FG.Charge",
    "FG.Comment",
    "F.CalibratedMz",
]
PQVALUES = [
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


def shrink(csv: str, sep: str, conditions: list[str], out: str) -> int:
    q = (
        pl.scan_csv(csv, separator=sep)
        .select(COLS+PQVALUES)
        .filter(pl.col("R.Condition").is_in(conditions))
    )
    df = q.collect()
    df.write_csv(out)
    return 0


def main():
    return shrink(INPUT, sep="\t", conditions=CONDITIONS, out=OUTPUT)


if __name__ == "__main__":
    exit(main())
