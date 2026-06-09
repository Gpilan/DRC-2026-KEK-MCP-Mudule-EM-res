#!/usr/bin/env bash

set -u

MAX_EVENT="${1:--1}"
FAILED_RUNS=()

if [[ ! -x "./calib_DRC" ]]; then
  echo "[ERROR] ./calib_DRC 실행 파일이 없습니다."
  echo "        먼저 analysis 디렉토리에서 컴파일하세요. (예: ./compile.sh calib_DRC.cc)"
  exit 1
fi

mkdir -p "./Calib"

run_one() {
  local run="$1"
  local aux_file="./AUX_ref/AUX_ref_Run_${run}.root"

  echo "============================================================"
  echo "[START] Run ${run} / calib_DRC (maxEvent=${MAX_EVENT})"

  if [[ ! -f "${aux_file}" ]]; then
    echo "[FAIL ] Run ${run} - AUX ref 파일 없음: ${aux_file}"
    FAILED_RUNS+=("${run}")
    return
  fi

  ./calib_DRC "${run}" "${MAX_EVENT}"
  local rc=$?
  if [[ ${rc} -ne 0 ]]; then
    echo "[FAIL ] Run ${run} (exit=${rc})"
    FAILED_RUNS+=("${run}")
  else
    echo "[DONE ] Run ${run}"
  fi
}

run_range() {
  local start="$1"
  local end="$2"
  local r
  for ((r=start; r<=end; r++)); do
    run_one "${r}"
  done
}

run_range 14127 14135
run_range 14222 14224
run_range 14226 14240

echo "============================================================"
if [[ ${#FAILED_RUNS[@]} -eq 0 ]]; then
  echo "[SUMMARY] calib_DRC: all success"
  exit 0
fi

echo "[SUMMARY] calib_DRC failed runs: ${FAILED_RUNS[*]}"
exit 1
