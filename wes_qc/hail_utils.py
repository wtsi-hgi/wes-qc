"""
A service module to work with Hail and Spark - spin up, shut down, work with temporary folders
"""

import sys
import pyspark
import hail as hl  # type: ignore
import shutil
import os


def clear_temp_folder(tmp_dir: str) -> None:
    if not tmp_dir.startswith("file://"):
        return
    tmp_dir = tmp_dir.replace("file://", "")
    print(f"=== Cleaning up temporary folder {tmp_dir}")
    shutil.rmtree(tmp_dir)
    os.makedirs(tmp_dir, exist_ok=True)


def init_hl(tmp_dir: str) -> pyspark.SparkContext:
    clear_temp_folder(tmp_dir)
    print("=== Checking for SparkContext ===")
    sc = pyspark.SparkContext.getOrCreate()

    # sc = pyspark.SparkContext()
    print("=== Initializing Hail ===")
    hl.init(sc=sc, tmp_dir=tmp_dir, idempotent=True)
    hl.default_reference("GRCh38")
    return sc


def stop_hl(sc: pyspark.SparkContext) -> None:
    hl.stop()
    sc.stop()


def currentFuncName(n: int = 0) -> str:
    return sys._getframe(n + 1).f_code.co_name


def trace_analysis_step(mt: hl.MatrixTable, description: str = "") -> hl.MatrixTable:
    """
    The function records the name of the caller function into the global annotation of the matrixtable
    """
    caller_name = currentFuncName(1)
    step_description = hl.struct(step_name=caller_name, step_description=description)
    print(f"=== Annotating step {caller_name}")
    if "analysis_steps" in mt.globals:
        steps = mt.globals.analysis_steps
        steps = steps.append(step_description)
    else:
        steps = hl.array([step_description])
    print_analysis_steps(steps)
    mt = mt.annotate_globals(analysis_steps=steps)
    return mt


def print_analysis_steps(analysis_steps: hl.ArrayExpression) -> None:
    for n, step in enumerate(hl.eval(analysis_steps)):
        print(f"Step {n:4}\t{step.step_name:20}\t{step.step_description}")
