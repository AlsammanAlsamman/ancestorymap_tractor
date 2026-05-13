import os
import sys

sys.path.append("utils")
from bioconfigme import get_analysis_value, get_default_resource, get_results_dir

PANEL_FILE = get_analysis_value("first_step.panel_file")
PANEL_OUT_DIR = get_analysis_value("first_step.output_dir")
RESULTS_DIR = get_results_dir()
LOG_DIR = os.path.join(RESULTS_DIR, "log")
DONE_DIR = os.path.join(RESULTS_DIR, "reference_panel_prep")

POPULATIONS = [str(p) for p in get_analysis_value("first_step.populations")]
ANC_LABELS = {str(k): str(v) for k, v in get_analysis_value("first_step.ancestry_labels").items()}
COMBINED_SUFFIX = "_".join(POPULATIONS)

if not POPULATIONS:
    raise ValueError("first_step.populations is empty; configure at least one reference population")

MISSING_LABELS = [pop for pop in POPULATIONS if pop not in ANC_LABELS]
if MISSING_LABELS:
    raise ValueError(
        f"Missing first_step.ancestry_labels for populations: {', '.join(MISSING_LABELS)}"
    )

POP_OUTPUTS = [os.path.join(PANEL_OUT_DIR, f"{pop}.samples.txt") for pop in POPULATIONS]
POP_ARGS = " ".join(
    f"--population {pop}:{ANC_LABELS[pop]}:{os.path.join(PANEL_OUT_DIR, pop + '.samples.txt')}"
    for pop in POPULATIONS
)

rule build_reference_sample_lists:
    input:
        panel=PANEL_FILE
    output:
        per_population=POP_OUTPUTS,
        combined=os.path.join(PANEL_OUT_DIR, f"{COMBINED_SUFFIX}.samples.txt"),
        sample_map=os.path.join(PANEL_OUT_DIR, f"{COMBINED_SUFFIX}.sample_map.txt"),
        done=os.path.join(DONE_DIR, "build_reference_sample_lists.done")
    log:
        os.path.join(LOG_DIR, "build_reference_sample_lists.log")
    resources:
        mem_mb=get_default_resource("mem_mb", 32000),
        time=get_default_resource("time", "00:30:00"),
        cores=get_default_resource("cores", 2)
    params:
        helper="scripts/example_helper.py",
        population_args=POP_ARGS,
    shell:
        """
        mkdir -p {PANEL_OUT_DIR} {LOG_DIR} {DONE_DIR}
        python {params.helper} \
            --panel {input.panel} \
            {params.population_args} \
            --combined {output.combined} \
            --sample-map {output.sample_map} \
            > {log} 2>&1
        touch {output.done}
        """
