import os
import re
import sys
import traceback

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_default_resource, get_results_dir, get_step5_successful_chromosomes


def _cfg(path, default=None):
    try:
        return get_analysis_value(path)
    except KeyError:
        return default


def _safe_name(value):
    value = str(value).strip()
    value = re.sub(r"[^A-Za-z0-9._-]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    if not value:
        raise ValueError("Collapsed group name resolved to an empty token")
    return value


def _parse_group(group_cfg, group_key):
    if not isinstance(group_cfg, dict):
        raise ValueError(f"collapsed_tractor.{group_key} must be a mapping")

    raw_name = str(group_cfg.get("name", "")).strip()
    ancestries = group_cfg.get("ancestries", [])

    if not raw_name:
        raise ValueError(f"collapsed_tractor.{group_key}.name must be non-empty")
    if not isinstance(ancestries, list) or not ancestries:
        raise ValueError(f"collapsed_tractor.{group_key}.ancestries must be a non-empty list")

    ancestry_labels = [str(item).strip() for item in ancestries if str(item).strip()]
    if len(ancestry_labels) != len(set(ancestry_labels)):
        raise ValueError(f"collapsed_tractor.{group_key}.ancestries contains duplicates: {ancestry_labels}")

    return {
        "name_raw": raw_name,
        "name_safe": _safe_name(raw_name),
        "ancestries": ancestry_labels,
    }


def _resolve_prefix(pattern, chr_name):
    return str(pattern).replace("{chr}", str(chr_name))


RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")

STEP7_DONE = os.path.join(RESULTS_DIR, "tractor", "extract_tracts.done")
CHROMOSOMES = [str(chrom) for chrom in get_step5_successful_chromosomes()]

SOURCE_DIR = str(_cfg("collapsed_tractor.input_dir", get_analysis_value("tractor.output_dir")))
OUTPUT_DIR = str(_cfg("collapsed_tractor.output_dir", os.path.join(RESULTS_DIR, "tractor_collapsed")))
OUTPUT_PREFIX_PATTERN = str(_cfg("collapsed_tractor.output_prefix_pattern", "chr{chr}.phased"))
CHUNK_ROWS = int(_cfg("collapsed_tractor.chunk_rows", 50000))

if CHUNK_ROWS <= 0:
    raise ValueError("collapsed_tractor.chunk_rows must be a positive integer")

GROUP1 = _parse_group(get_analysis_value("collapsed_tractor.group1"), "group1")
GROUP2 = _parse_group(get_analysis_value("collapsed_tractor.group2"), "group2")

if GROUP1["name_safe"] == GROUP2["name_safe"]:
    raise ValueError(
        "collapsed_tractor.group1.name and collapsed_tractor.group2.name must be different after filename normalization"
    )

label_to_code = get_analysis_value("step5.sample_map_label_codes")
if not isinstance(label_to_code, dict) or not label_to_code:
    raise ValueError("step5.sample_map_label_codes must be a non-empty mapping")
label_to_code = {str(label): str(code) for label, code in label_to_code.items()}

all_requested = GROUP1["ancestries"] + GROUP2["ancestries"]
unknown = sorted({label for label in all_requested if label not in label_to_code})
if unknown:
    raise ValueError(
        "collapsed_tractor groups reference unknown ancestry labels: "
        + ", ".join(unknown)
        + "; available labels: "
        + ", ".join(sorted(label_to_code.keys()))
    )

overlap = sorted(set(GROUP1["ancestries"]).intersection(set(GROUP2["ancestries"])))
if overlap:
    raise ValueError(
        "collapsed_tractor.group1 and collapsed_tractor.group2 must be disjoint; overlap: " + ", ".join(overlap)
    )

GROUP1_CODES = [label_to_code[label] for label in GROUP1["ancestries"]]
GROUP2_CODES = [label_to_code[label] for label in GROUP2["ancestries"]]

GROUP1_DOSAGE_PATTERN = os.path.join(OUTPUT_DIR, OUTPUT_PREFIX_PATTERN + f".{GROUP1['name_safe']}.dosage.txt")
GROUP2_DOSAGE_PATTERN = os.path.join(OUTPUT_DIR, OUTPUT_PREFIX_PATTERN + f".{GROUP2['name_safe']}.dosage.txt")
GROUP1_HAPCOUNT_PATTERN = os.path.join(OUTPUT_DIR, OUTPUT_PREFIX_PATTERN + f".{GROUP1['name_safe']}.hapcount.txt")
GROUP2_HAPCOUNT_PATTERN = os.path.join(OUTPUT_DIR, OUTPUT_PREFIX_PATTERN + f".{GROUP2['name_safe']}.hapcount.txt")


def source_dosage_path(chr_name, ancestry_code):
    return os.path.join(SOURCE_DIR, f"chr{chr_name}.phased.anc{ancestry_code}.dosage.txt")


def source_hapcount_path(chr_name, ancestry_code):
    return os.path.join(SOURCE_DIR, f"chr{chr_name}.phased.anc{ancestry_code}.hapcount.txt")


def _collapse_ancestry_tables_to_path(paths, table_kind, chr_name, output_path, chunk_rows):
    import pandas as pd

    if not paths:
        raise ValueError(f"No input files supplied for collapsed {table_kind}, chr{chr_name}")

    def _iter_chunks(path):
        # Inputs are tab-delimited; stream chunks to avoid whole-file memory spikes.
        return pd.read_csv(path, sep="\t", low_memory=False, chunksize=chunk_rows)

    def _coerce_numeric(frame, value_cols, path, chunk_idx):
        converted = frame[value_cols].apply(pd.to_numeric, errors="coerce")
        bad = converted.isna() & frame[value_cols].notna()
        if bad.any().any():
            bad_rows = bad.any(axis=1)
            first_row = int(bad_rows.idxmax()) if bad_rows.any() else -1
            first_col = next((col for col in value_cols if bool(bad.loc[first_row, col])), "<unknown>")
            raise ValueError(
                f"Non-numeric value while collapsing {table_kind} for chr{chr_name}: "
                f"input={path}, chunk={chunk_idx}, row={first_row}, col={first_col}"
            )
        return converted

    readers = [_iter_chunks(path) for path in paths]
    colnames = None
    key_cols = None
    value_cols = None
    chunk_idx = 0
    rows_written = 0

    while True:
        chunks = []
        exhausted = []
        for reader in readers:
            try:
                chunks.append(next(reader))
                exhausted.append(False)
            except StopIteration:
                chunks.append(None)
                exhausted.append(True)

        if all(exhausted):
            break
        if any(exhausted):
            raise ValueError(
                f"Row-count mismatch while collapsing {table_kind} for chr{chr_name}; one file ended early"
            )

        base_chunk = chunks[0]
        if colnames is None:
            colnames = list(base_chunk.columns)
            key_cols = [col for col in ["CHROM", "POS", "ID", "REF", "ALT"] if col in colnames]
            value_cols = [col for col in colnames if col not in key_cols]
            if not value_cols:
                raise ValueError(f"No numeric value columns found while collapsing {table_kind} for chr{chr_name}")

        if list(base_chunk.columns) != colnames:
            raise ValueError(
                f"Column mismatch while collapsing {table_kind} for chr{chr_name} in input[0]"
            )

        summed_values = _coerce_numeric(base_chunk, value_cols, paths[0], chunk_idx)

        for idx, (path, frame) in enumerate(zip(paths[1:], chunks[1:]), start=1):
            if list(frame.columns) != colnames:
                raise ValueError(
                    f"Column mismatch while collapsing {table_kind} for chr{chr_name} between input[0] and input[{idx}]"
                )

            if key_cols and not frame[key_cols].equals(base_chunk[key_cols]):
                raise ValueError(
                    f"Variant row mismatch while collapsing {table_kind} for chr{chr_name} between input[0] and input[{idx}] in chunk {chunk_idx}"
                )

            frame_values = _coerce_numeric(frame, value_cols, path, chunk_idx)
            summed_values = summed_values.add(frame_values, fill_value=0)

        out_chunk = base_chunk.copy()
        out_chunk[value_cols] = summed_values
        out_chunk = out_chunk[colnames]

        out_chunk.to_csv(
            output_path,
            sep="\t",
            index=False,
            mode="w" if chunk_idx == 0 else "a",
            header=(chunk_idx == 0),
        )
        rows_written += len(out_chunk)
        chunk_idx += 1

    if chunk_idx == 0:
        raise ValueError(f"No rows found while collapsing {table_kind} for chr{chr_name}")

    return rows_written


rule collapse_ancestry_dosage_counts_chr:
    input:
        step7_done=STEP7_DONE,
        g1_dosage=lambda wildcards: [source_dosage_path(wildcards.chr, code) for code in GROUP1_CODES],
        g2_dosage=lambda wildcards: [source_dosage_path(wildcards.chr, code) for code in GROUP2_CODES],
        g1_hapcount=lambda wildcards: [source_hapcount_path(wildcards.chr, code) for code in GROUP1_CODES],
        g2_hapcount=lambda wildcards: [source_hapcount_path(wildcards.chr, code) for code in GROUP2_CODES],
    output:
        g1_dosage=os.path.join(OUTPUT_DIR, OUTPUT_PREFIX_PATTERN + f".{GROUP1['name_safe']}.dosage.txt"),
        g2_dosage=os.path.join(OUTPUT_DIR, OUTPUT_PREFIX_PATTERN + f".{GROUP2['name_safe']}.dosage.txt"),
        g1_hapcount=os.path.join(OUTPUT_DIR, OUTPUT_PREFIX_PATTERN + f".{GROUP1['name_safe']}.hapcount.txt"),
        g2_hapcount=os.path.join(OUTPUT_DIR, OUTPUT_PREFIX_PATTERN + f".{GROUP2['name_safe']}.hapcount.txt"),
    log:
        os.path.join(LOG_DIR, "collapse_ancestry_dosage_counts_chr{chr}.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    run:
        os.makedirs(OUTPUT_DIR, exist_ok=True)
        os.makedirs(LOG_DIR, exist_ok=True)

        with open(log[0], "w", encoding="utf-8") as handle:
            handle.write(
                "Starting collapsed ancestry build for chr{}\n"
                "group1={} (labels={}, codes={})\n"
                "group2={} (labels={}, codes={})\n".format(
                    wildcards.chr,
                    GROUP1["name_raw"],
                    ",".join(GROUP1["ancestries"]),
                    ",".join(GROUP1_CODES),
                    GROUP2["name_raw"],
                    ",".join(GROUP2["ancestries"]),
                    ",".join(GROUP2_CODES),
                )
            )

            try:
                handle.write(f"chunk_rows={CHUNK_ROWS}\n")

                g1_dose_rows = _collapse_ancestry_tables_to_path(
                    list(input.g1_dosage), "dosage/group1", wildcards.chr, output.g1_dosage, CHUNK_ROWS
                )
                g2_dose_rows = _collapse_ancestry_tables_to_path(
                    list(input.g2_dosage), "dosage/group2", wildcards.chr, output.g2_dosage, CHUNK_ROWS
                )
                g1_hap_rows = _collapse_ancestry_tables_to_path(
                    list(input.g1_hapcount), "hapcount/group1", wildcards.chr, output.g1_hapcount, CHUNK_ROWS
                )
                g2_hap_rows = _collapse_ancestry_tables_to_path(
                    list(input.g2_hapcount), "hapcount/group2", wildcards.chr, output.g2_hapcount, CHUNK_ROWS
                )

                handle.write(
                    "Rows written: g1_dosage={}, g2_dosage={}, g1_hapcount={}, g2_hapcount={}\n".format(
                        g1_dose_rows, g2_dose_rows, g1_hap_rows, g2_hap_rows
                    )
                )
                handle.write("Completed collapsed ancestry build successfully\n")
            except Exception as exc:
                handle.write(f"ERROR: {exc}\n")
                handle.write(traceback.format_exc())
                raise


rule collapse_ancestry_dosage_counts_done:
    input:
        step7_done=STEP7_DONE,
        g1_dosage=expand(GROUP1_DOSAGE_PATTERN, chr=CHROMOSOMES),
        g2_dosage=expand(GROUP2_DOSAGE_PATTERN, chr=CHROMOSOMES),
        g1_hapcount=expand(GROUP1_HAPCOUNT_PATTERN, chr=CHROMOSOMES),
        g2_hapcount=expand(GROUP2_HAPCOUNT_PATTERN, chr=CHROMOSOMES),
    output:
        done=os.path.join(OUTPUT_DIR, "collapse_ancestry_dosage_counts.done"),
    log:
        os.path.join(LOG_DIR, "collapse_ancestry_dosage_counts_done.log"),
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2),
    params:
        g1_name=GROUP1["name_safe"],
        g2_name=GROUP2["name_safe"],
    shell:
        """
        mkdir -p {OUTPUT_DIR}
        touch {output.done}
        echo "Collapsed ancestry dosage/hapcount build complete: {params.g1_name} vs {params.g2_name}" > {log}
        """
