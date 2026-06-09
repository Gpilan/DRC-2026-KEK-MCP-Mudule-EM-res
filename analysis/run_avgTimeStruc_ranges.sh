#!/usr/bin/env bash

set -u

MAX_EVENT="${1:--1}"
FAILED_RUNS=()

if [[ ! -x "./avgTimeStruc" ]]; then
  echo "[ERROR] ./avgTimeStruc 실행 파일이 없습니다."
  echo "        먼저 analysis 디렉토리에서 컴파일하세요. (예: ./compile.sh avgTimeStruc.cc)"
  exit 1
fi

run_one() {
  local run="$1"
  echo "============================================================"
  echo "[START] Run ${run} (maxEvent=${MAX_EVENT})"
  ./avgTimeStruc "${run}" "${MAX_EVENT}"
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
  echo "[SUMMARY] 모든 런 처리 완료"
  exit 0
fi

echo "[SUMMARY] 실패한 런: ${FAILED_RUNS[*]}"
exit 1
