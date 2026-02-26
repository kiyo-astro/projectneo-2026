#--------------------------------------------------------------------------------------------------#
# api_data.py | ASTR 302 WINTER 2026 Project NEO                                                   #
# Developed by Kiyoaki Okudaira * University of Washington                                         #
#--------------------------------------------------------------------------------------------------#
# Description                                                                                      #
#--------------------------------------------------------------------------------------------------#
# Retrieve asteroids data set from NASA SENTRY / ESA AEGIS API and parse data set                  #
#--------------------------------------------------------------------------------------------------#
# History                                                                                          #
#--------------------------------------------------------------------------------------------------#
# coding 2026.02.03: 1st coding                                                                    #
#--------------------------------------------------------------------------------------------------#

#--------------------------------------------------------------------------------------------------#
# Libraries                                                                                        #
#--------------------------------------------------------------------------------------------------#
import requests
from astropy.table import Table
import astropy.units as u
import re

from typing import Dict, List, Tuple

#--------------------------------------------------------------------------------------------------#
# Global constants                                                                                 #
#--------------------------------------------------------------------------------------------------#
NASA_SENTRY_API_ROOT = "https://ssd-api.jpl.nasa.gov/sentry.api"
ESA_AEGIS_API_ROOT = "https://neo.ssa.esa.int/PSDB-portlet/download"

#--------------------------------------------------------------------------------------------------#
# Main                                                                                             #
#--------------------------------------------------------------------------------------------------#
def retrieve_NASA_sentry_data(
        ip_min : bool = None
        ):
    """
    Retrieve asteroids data set from NASA SENTRY API

    Parameters
    ----------
    ip_min: `float` or `int` or `bool`
        Minimum impact probability. Default is None (not filter any data by impact probability)

    Returns
    -------
    query_result: `dict`
        Asteroids data set

    Notes
    -----
        (c) 2026 Kiyoaki Okudaira - University of Washington
    """
    query_params = {"all": "1"}
    if ip_min is not None:
        query_params["ip-min"] = f"{ip_min:g}"

    r = requests.get(NASA_SENTRY_API_ROOT, params=query_params)
    r.raise_for_status()
    query_result = r.json()

    if "error" in query_result:
        raise RuntimeError(f"Sentry API error: {query_result['error']}")

    return query_result

def retrieve_NASA_sentry_summary(
        ip_min : bool = None
        ):
    """
    Retrieve asteroids summary data set from NASA SENTRY API

    Parameters
    ----------
    ip_min: `float` or `int` or `bool`
        Minimum impact probability. Default is None (not filter any data by impact probability)

    Returns
    -------
    query_result: `dict`
        Asteroids data set

    Notes
    -----
        (c) 2026 Kiyoaki Okudaira - University of Washington
    """
    query_params = {}
    if ip_min is not None:
        query_params["ip-min"] = f"{ip_min:g}"

    r = requests.get(NASA_SENTRY_API_ROOT, params=query_params)
    r.raise_for_status()
    query_result = r.json()

    if "error" in query_result:
        raise RuntimeError(f"Sentry API error: {query_result['error']}")

    return query_result

def retrieve_ESA_AEGIS_risk_list():
    """
    Retrieve asteroids risk list from ESA AEGIS API

    Returns
    -------
    query_result: `str`
        Asteroids risk list

    Notes
    -----
        (c) 2026 Kiyoaki Okudaira - University of Washington
    """

    r = requests.get(ESA_AEGIS_API_ROOT + "?file=esa_risk_list")
    r.raise_for_status()
    query_result = r.text

    if "error" in query_result:
        raise RuntimeError(f"Sentry API error: {query_result['error']}")

    return query_result

def dict_list2table(
        data : list
        ):
    """
    Parse list with dict-format dataset and return data with astropy Table format

    Parameters
    ----------
    data: `list`
        List with dict-format dataset

    Returns
    -------
    astropy_table: `astropy.table.table.Table`
        Data with astropy Table format

    Notes
    -----
        (c) 2026 Kiyoaki Okudaira - University of Washington
    """
    all_keys = {k for row in data for k in row.keys()}
    rows = [{k: row.get(k, None) for k in all_keys} for row in data]
    astropy_table = Table(rows=rows)
    return astropy_table

def _to_int(x):
    x = x.strip()
    if x == "" or x.lower() in {"n/a", "na"}:
        return None
    try:
        return int(x)
    except Exception:
        return None

def _to_float(x):
    x = x.strip()
    if x == "" or x.lower() in {"n/a", "na"}:
        return None
    try:
        return float(x)
    except Exception:
        return None

def parse_ESA_AEGIS_risk_list(
        query_result : str
        ):
    """
    Parse asteroids risk list from ESA AEGIS API and return data with astropy Table format

    Parameters
    ----------
    data: `str`
        Asteroids risk list from ESA AEGIS API

    Returns
    -------
    astropy_table: `astropy.table.table.Table`
        Data with astropy Table format

    Notes
    -----
        (c) 2026 Kiyoaki Okudaira - University of Washington
    """
    lines = [ln.rstrip("\n") for ln in query_result.splitlines()]

    # Last Update rows
    last_update = None
    for ln in lines:
        if ln.startswith("Last Update:"):
            last_update = ln.replace("Last Update:", "").strip()
            break

    # Header row
    start_idx = None
    for i, ln in enumerate(lines):
        if ln.strip().startswith("AAAAAAAA"):
            start_idx = i + 1
            break
    if start_idx is None:
        raise RuntimeError("Could not find the header template line (AAAAAAAA...). Format may have changed.")

    # Data rows
    data_lines = [ln for ln in lines[start_idx:] if ln.strip() != ""]

    # Parse
    rows = []
    for ln in data_lines:
        parts = [p.strip() for p in ln.split("|")]

        if len(parts) < 4:
            continue

        # Object
        obj_tokens = parts[0].split()
        num_des = obj_tokens[0] if obj_tokens else ""
        name = " ".join(obj_tokens[1:]) if len(obj_tokens) > 1 else ""

        # Diameter
        diameter_m = _to_float(parts[1]) if len(parts[1]) >= 1 else None
        diameter_est = ("*" in parts[2])

        # VI Max
        vi_datetime_utc_tokens = parts[3].split()
        vi_datetime_utc = " ".join(vi_datetime_utc_tokens[0:2]) if len(vi_datetime_utc_tokens) >= 2 else ""
        ip_max = _to_float(parts[4]) if len(parts[4]) >= 1 else None
        ps_max = _to_float(parts[5]) if len(parts[5]) >= 1 else None
        ts = _to_int(parts[6]) if len(parts[6]) >= 1 else None
        vel_kms = _to_float(parts[7]) if len(parts[7]) >= 1 else None

        # VIs
        years = parts[8] if len(parts[8]) >= 1 else ""
        ip_cum = _to_float(parts[9]) if len(parts[9]) >= 1 else None
        ps_cum = _to_float(parts[10]) if len(parts[10]) >= 1 else None

        rows.append(
            (num_des, name, diameter_m, diameter_est, vi_datetime_utc, ip_max, ps_max, ts, vel_kms, years, ip_cum, ps_cum)
        )

    astropy_table = Table(
        rows=rows,
        names=[
            "des", "name",
            "diameter_m", "diameter_estimated",
            "vi_datetime_utc",
            "ip_max", "ps_max", "ts", "vel_km_s",
            "years_range", "ip_cum", "ps_cum",
        ],
        meta={"last_update": last_update},
    )
    
    return astropy_table

def retrieve_ESA_AEGIS_data(
        des : str
        ):
    """
    Retrieve asteroid data set from ESA AEGIS API

    Parameters
    ----------
    des: `str`
        Asteroid designator

    Returns
    -------
    query_result: `str`
        Asteroid data set

    Notes
    -----
        (c) 2026 Kiyoaki Okudaira - University of Washington
    """

    r = requests.get(f"https://neo.ssa.esa.int/PSDB-portlet/download?file={des}.risk")
    r.raise_for_status()
    query_result = r.text

    if "error" in query_result:
        raise RuntimeError(f"Sentry API error: {query_result['error']}")

    return query_result

def parse_ESA_AEGIS_data(
        query_result : str,
        risk_list_row : Table = None
        ):

    lines = query_result.splitlines()

    # metadata
    sep_idx = None
    for i, ln in enumerate(lines):
        if re.match(r"^\s*-{10,}\s*$", ln):
            sep_idx = i
            break
    if sep_idx is None:
        raise RuntimeError("Undefined file format.")

    data_lines = []
    for ln in lines[sep_idx + 1 :]:
        if not ln.strip(): 
            break
        if "<" in ln and ">" in ln:
            break
        data_lines.append(ln.strip())

    names = (
        "des", "name", "diameter_m", "diameter_estimated",
        "vi_datetime_utc", "ip_max", "ps_max",
        "vel_km_s", "years_range", "ip_cum", "ps_cum",
        "date", "mjd", "sigma", "sigimp",
        "dist_re", "width_re", "stretch_re",
        "p_re", "exp_energy_mt", "ps", "ts"
    )

    if not data_lines:
        return Table(rows=[], names=names)

    rows = []
    for ln in data_lines:
        parts = ln.split()

        if len(parts) < 12 or parts[5] != "+/-":
            raise RuntimeError(f"Unexpected table row format:\n{ln}")

        rows.append((
            risk_list_row["des"],
            risk_list_row["name"],
            risk_list_row["diameter_m"],
            risk_list_row["diameter_estimated"],
            risk_list_row["vi_datetime_utc"],
            risk_list_row["ip_max"],
            risk_list_row["ps_max"],
            risk_list_row["vel_km_s"],
            risk_list_row["years_range"],
            risk_list_row["ip_cum"],
            risk_list_row["ps_cum"],
            parts[0],                 # date
            _to_float(parts[1]),          # mjd
            _to_float(parts[2]),          # sigma
            _to_float(parts[3]),          # sigimp
            _to_float(parts[4]),          # dist
            _to_float(parts[6]),          # width
            _to_float(parts[7]),          # stretch
            _to_float(parts[8]),          # p_RE
            _to_float(parts[9]),          # exp. en.
            _to_float(parts[10]),         # PS
            None if _to_float(parts[11]) is None else int(_to_float(parts[11])),    # TS
        ))

    return Table(rows=rows, names=names, dtype=[object]*len(names))

def retrieve_ESA_AEGIS_NEA_cat():
    """
    Retrieve catalogue of NEAs including their Keplerian orbital elements 

    Returns
    -------
    query_result: `str`
        Asteroids risk list

    Notes
    -----
        (c) 2026 Kiyoaki Okudaira - University of Washington
    """

    r = requests.get(ESA_AEGIS_API_ROOT + "?file=neo_kc.cat")
    r.raise_for_status()
    query_result = r.text

    if "error" in query_result:
        raise RuntimeError(f"Sentry API error: {query_result['error']}")

    return query_result

def parse_ESA_AEGIS_NEA_cat(
    text: str,
    strict: bool = True,
    add_units: bool = True,
) -> Tuple[Table, Dict[str, str]]:
    """
    Parse asteroids catalogue of NEAs and return data with astropy Table format

    Parameters
    ----------
    text : str
        Raw text content downloaded from ESA AEGIS.
    strict : bool, default True
        If True, raise ValueError on malformed data lines.
        If False, skip malformed lines.
    add_units : bool, default True
        If True, attach common units to columns (a in AU, angles in deg, epoch in d).

    Returns
    -------
    table : astropy.table.Table
        Parsed table.
    header_kv : dict
        Key/value pairs parsed from header lines like "key = value ! comment".
    """
    lines = text.splitlines()

    # ---- split header / body ----
    header_kv: Dict[str, str] = {}
    body_start = None
    for idx, raw in enumerate(lines):
        line = raw.strip()
        if not line:
            continue
        if line == "END_OF_HEADER":
            body_start = idx + 1
            break
        if line.startswith("!"):
            continue
        # Parse "key = value" (optionally followed by "! comment")
        if "=" in line:
            left, right = line.split("=", 1)
            key = left.strip()
            val = right.split("!", 1)[0].strip()
            if key:
                header_kv[key] = val

    if body_start is None:
        raise ValueError("END_OF_HEADER not found in text.")

    # ---- prepare columns ----
    names = [
        "name",
        "epoch_mjd",
        "a_au",
        "e",
        "i_deg",
        "raan_deg",   # long. node
        "aop_deg",    # arg. peric.
        "M_deg",      # mean anomaly
        "H",
        "G",
        "nongrav",
    ]

    rows: List[tuple] = []
    bad_lines: List[str] = []

    for raw in lines[body_start:]:
        line = raw.strip()
        if (not line) or line.startswith("!"):
            continue

        parts = line.split()
        if len(parts) < 11:
            msg = f"Malformed data line (expected >=11 fields, got {len(parts)}): {line}"
            if strict:
                raise ValueError(msg)
            bad_lines.append(line)
            continue

        # take first 11 fields; if extra tokens exist, ignore them
        parts = parts[:11]

        try:
            row = (
                parts[0],           # name (string)
                float(parts[1]),    # epoch mjd
                float(parts[2]),    # a (AU)
                float(parts[3]),    # e
                float(parts[4]),    # i (deg)
                float(parts[5]),    # long. node (deg)
                float(parts[6]),    # arg. peric. (deg)
                float(parts[7]),    # mean anomaly (deg)
                float(parts[8]),    # H
                float(parts[9]),    # G
                float(parts[10]),   # non-grav param (often 0)
            )
            rows.append(row)
        except Exception as e:
            msg = f"Failed to parse line: {line} ({e})"
            if strict:
                raise ValueError(msg)
            bad_lines.append(line)

    tbl = Table(rows=rows, names=names)

    if add_units:
        # These are "display/metadata" units; values remain plain floats in the table.
        tbl["epoch_mjd"].unit = u.d
        tbl["a_au"].unit = u.AU
        tbl["i_deg"].unit = u.deg
        tbl["raan_deg"].unit = u.deg
        tbl["aop_deg"].unit = u.deg
        tbl["M_deg"].unit = u.deg
        # H, G are dimensionless; nongrav could be something else depending on format,
        # but typically 0 in your sample so leave as dimensionless.

    # If you want to keep track of skipped lines in non-strict mode:
    if (not strict) and bad_lines:
        tbl.meta["skipped_lines"] = bad_lines

    # Also store header info in meta
    tbl.meta.update(header_kv)

    return tbl, header_kv

#--------------------------------------------------------------------------------------------------#
# Test                                                                                             #
#--------------------------------------------------------------------------------------------------#