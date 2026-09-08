"""Summarize program CSV output without third-party dependencies.

Usage: python tools/summarize_results.py build/multi.csv --start 2 --end 6
Position errors are relative to the RINEX approximate coordinate, not ground truth.
Velocity RMS assumes the supplied station is static.
"""
import argparse
import csv
import json
import math
from collections import Counter


def summarize(path, start=0, end=24):
    with open(path, encoding="utf-8-sig", newline="") as source:
        rows = [r for r in csv.DictReader(source)
                if start <= float(r["tow"]) % 86400 / 3600 <= end]
    metrics = {}
    for column in ("de", "dn", "du", "ve", "vn", "vu"):
        values = [float(r[column]) for r in rows if math.isfinite(float(r[column]))]
        metrics[column] = {
            "count": len(values),
            "rms": math.sqrt(sum(v*v for v in values) / len(values)) if values else None,
            "mean": sum(values) / len(values) if values else None,
            "max_abs": max(map(abs, values)) if values else None,
        }
    return {"epochs": len(rows), "statuses": dict(Counter(r["status"] for r in rows)),
            "metrics": metrics}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv")
    parser.add_argument("--start", type=float, default=0)
    parser.add_argument("--end", type=float, default=24)
    args = parser.parse_args()
    if not 0 <= args.start <= args.end <= 24:
        parser.error("expected 0 <= start <= end <= 24")
    print(json.dumps(summarize(args.csv, args.start, args.end), indent=2, allow_nan=False))
