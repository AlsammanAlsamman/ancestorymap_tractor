#!/usr/bin/env python3
import argparse
import csv
import math
import sys


def to_float(value):
    if value is None:
        return None
    text = str(value).strip()
    if text == "" or text.upper() == "NA":
        return None
    try:
        return float(text)
    except ValueError:
        return None


def read_pvar_alleles(pvar_path):
    allele_map = {}
    with open(pvar_path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"ID", "REF", "ALT"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"PVAR missing required columns {sorted(missing)}: {pvar_path}")
        for row in reader:
            snp_id = str(row["ID"])
            allele_map[snp_id] = (str(row["REF"]), str(row["ALT"]))
    return allele_map


def resolve_ea_freq(ea, ref, alt, alt_freq, ref_freq):
    if ea == alt:
        return alt_freq, ref
    if ea == ref:
        return ref_freq, alt
    return None, None


def fmt(value):
    if value is None:
        return ""
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return ""
    return str(value)


def main():
    parser = argparse.ArgumentParser(description="Build standardized effect-allele GWAS tables.")
    parser.add_argument("--input-gwas", required=True, help="Per-ancestry GWAS table with case/control frequencies.")
    parser.add_argument("--input-pvar", required=True, help="Per-ancestry merged PVAR with REF/ALT.")
    parser.add_argument("--output-full", required=True, help="Output full table.")
    parser.add_argument("--output-significant", required=True, help="Output subset table with P < threshold.")
    parser.add_argument("--p-threshold", type=float, default=5e-3, help="Significance threshold for subset output.")
    args = parser.parse_args()

    allele_map = read_pvar_alleles(args.input_pvar)

    out_fields = [
        "RSID",
        "CHR",
        "POS",
        "REF",
        "ALT",
        "EA",
        "NEA",
        "P",
        "BETA",
        "OR",
        "SE",
        "OR_SE",
        "EA_FREQ_CASE",
        "EA_FREQ_CONTROL",
        "EA_FREQ_TOTAL",
    ]

    skipped_missing_alleles = 0
    skipped_unmatched_effect = 0
    written = 0
    written_sig = 0

    with open(args.input_gwas, "r", encoding="utf-8") as src, \
        open(args.output_full, "w", encoding="utf-8", newline="") as dst_full, \
        open(args.output_significant, "w", encoding="utf-8", newline="") as dst_sig:

        reader = csv.DictReader(src, delimiter="\t")
        required = {
            "ID",
            "CHR",
            "POS",
            "A1",
            "P",
            "ALT_FREQ_CASE",
            "REF_FREQ_CASE",
            "OBS_CT_CASE",
            "ALT_FREQ_CTRL",
            "REF_FREQ_CTRL",
            "OBS_CT_CTRL",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Input GWAS missing required columns {sorted(missing)}: {args.input_gwas}")

        full_writer = csv.DictWriter(dst_full, fieldnames=out_fields, delimiter="\t")
        sig_writer = csv.DictWriter(dst_sig, fieldnames=out_fields, delimiter="\t")
        full_writer.writeheader()
        sig_writer.writeheader()

        for row in reader:
            snp_id = str(row["ID"])
            if snp_id not in allele_map:
                skipped_missing_alleles += 1
                continue

            ref, alt = allele_map[snp_id]
            ea = str(row["A1"])

            alt_freq_case = to_float(row.get("ALT_FREQ_CASE"))
            ref_freq_case = to_float(row.get("REF_FREQ_CASE"))
            alt_freq_ctrl = to_float(row.get("ALT_FREQ_CTRL"))
            ref_freq_ctrl = to_float(row.get("REF_FREQ_CTRL"))
            obs_case = to_float(row.get("OBS_CT_CASE"))
            obs_ctrl = to_float(row.get("OBS_CT_CTRL"))

            ea_freq_case, nea = resolve_ea_freq(ea, ref, alt, alt_freq_case, ref_freq_case)
            ea_freq_ctrl, _ = resolve_ea_freq(ea, ref, alt, alt_freq_ctrl, ref_freq_ctrl)
            if nea is None:
                skipped_unmatched_effect += 1
                continue

            ea_freq_total = None
            if ea_freq_case is not None and ea_freq_ctrl is not None and obs_case is not None and obs_ctrl is not None:
                denom = obs_case + obs_ctrl
                if denom > 0:
                    ea_freq_total = (ea_freq_case * obs_case + ea_freq_ctrl * obs_ctrl) / denom
            elif ea_freq_case is not None and ea_freq_ctrl is None:
                ea_freq_total = ea_freq_case
            elif ea_freq_ctrl is not None and ea_freq_case is None:
                ea_freq_total = ea_freq_ctrl

            p_value = to_float(row.get("P"))
            beta = row.get("BETA", "")
            odds_ratio = row.get("OR", "")
            se = row.get("SE", "")
            or_se = se if str(odds_ratio).strip() not in {"", "NA"} else ""

            out_row = {
                "RSID": snp_id,
                "CHR": row.get("CHR", ""),
                "POS": row.get("POS", ""),
                "REF": ref,
                "ALT": alt,
                "EA": ea,
                "NEA": nea,
                "P": fmt(p_value),
                "BETA": beta,
                "OR": odds_ratio,
                "SE": se,
                "OR_SE": or_se,
                "EA_FREQ_CASE": fmt(ea_freq_case),
                "EA_FREQ_CONTROL": fmt(ea_freq_ctrl),
                "EA_FREQ_TOTAL": fmt(ea_freq_total),
            }

            full_writer.writerow(out_row)
            written += 1
            if p_value is not None and p_value < args.p_threshold:
                sig_writer.writerow(out_row)
                written_sig += 1

    print(
        f"Wrote {written} rows to {args.output_full}; "
        f"Wrote {written_sig} rows to {args.output_significant}; "
        f"Skipped missing_alleles={skipped_missing_alleles}, unmatched_effect={skipped_unmatched_effect}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
