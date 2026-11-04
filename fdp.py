#!/usr/bin/env python3

# 2026 (c) Micha Johannes Birklbauer
# https://github.com/michabirklbauer/
# micha.birklbauer@gmail.com

from __future__ import annotations

LIB_SIZE_SPYOGENES = 4503 + 13509  # target precursors + decoy precursors
LIB_SIZE_CELEGANS = 355 + 1065  # target precursors + decoy precursors
FDP_R = LIB_SIZE_CELEGANS / LIB_SIZE_SPYOGENES


def lower_bound(
    nr_entrapment_hits: int, nr_target_hits: int, as_percentage: bool = True
) -> float:
    r"""Calculates the FDP lower bound.

    Parameters
    ----------
    nr_entrapment_hits : int
        The number of entrapment hits.
    nr_target_hits : int
        The number of target hits.
    as_percentage : bool, default = True
        If 'True' (default) returns the FDP lower bound as a percentage, otherwise
        the FDP lower bound is returned as a fraction.

    Returns
    -------
    float
        The FDP lower bound.
    """
    fdp_lower_bound = nr_entrapment_hits / (nr_target_hits + nr_entrapment_hits)
    if as_percentage:
        return fdp_lower_bound * 100.0
    return fdp_lower_bound


def upper_bound(
    nr_entrapment_hits: int,
    nr_target_hits: int,
    ratio: float = FDP_R,
    as_percentage: bool = True,
) -> float:
    r"""Calculates the FDP upper bound.

    Parameters
    ----------
    nr_entrapment_hits : int
        The number of entrapment hits.
    nr_target_hits : int
        The number of target hits.
    ratio : float, default = FDP_R (the ratio between C. elegans and S. pyogenes)
        The ratio between the entrapment database size and the target database size.
    as_percentage : bool, default = True
        If 'True' (default) returns the FDP upper bound as a percentage, otherwise
        the FDP upper bound is returned as a fraction.

    Returns
    -------
    float
        The FDP upper bound.
    """
    fdp_upper_bound = (nr_entrapment_hits * (1 + 1 / ratio)) / (
        nr_target_hits + nr_entrapment_hits
    )
    if as_percentage:
        return fdp_upper_bound * 100.0
    return fdp_upper_bound
