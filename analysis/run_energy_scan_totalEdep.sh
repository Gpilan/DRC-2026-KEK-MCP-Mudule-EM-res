#!/usr/bin/env bash

set -u

if [[ ! -x "./energy_scan_totalEdep" ]]; then
  echo "[ERROR] ./energy_scan_totalEdep 실행 파일이 없습니다."
  echo "        먼저 analysis 디렉토리에서 컴파일하세요. (예: ./compile.sh energy_scan_totalEdep.cc)"
  exit 1
fi

echo "============================================================"
echo "[START] energy_scan_totalEdep"
./energy_scan_totalEdep
rc=$?

if [[ ${rc} -ne 0 ]]; then
  echo "[FAIL ] energy_scan_totalEdep (exit=${rc})"
  exit "${rc}"
fi

echo "[DONE ] energy_scan_totalEdep"
exit 0
