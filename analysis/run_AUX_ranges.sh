#!/usr/bin/env bash

set -u

MAX_EVENT="${1:--1}"
FAILED_AUX=()

if [[ ! -x "./draw_AUX" ]]; then
  echo "[ERROR] ./draw_AUX 실행 파일이 없습니다."
  echo "        먼저 analysis 디렉토리에서 컴파일하세요. (예: ./compile.sh draw_AUX.cc)"
  exit 1
fi

run_one() {
  local run="$1"

  echo "============================================================"

  echo "[RUN  ] ${run} / draw_AUX (maxEvent=${MAX_EVENT})"
  ./draw_AUX "${run}" "${MAX_EVENT}"
  local rc_aux=$?
  if [[ ${rc_aux} -ne 0 ]]; then
    echo "[FAIL ] draw_AUX Run ${run} (exit=${rc_aux})"
    FAILED_AUX+=("${run}")
  else
    echo "[DONE ] draw_AUX Run ${run}"
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

if [[ ${#FAILED_AUX[@]} -eq 0 ]]; then
  echo "[SUMMARY] draw_AUX: all success"
else
  echo "[SUMMARY] draw_AUX failed runs: ${FAILED_AUX[*]}"
fi

if [[ ${#FAILED_AUX[@]} -eq 0 ]]; then
  exit 0
fi

exit 1
