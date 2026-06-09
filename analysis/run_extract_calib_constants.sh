#!/usr/bin/env bash

set -u

if [[ ! -x "./extract_calib_constants" ]]; then
  echo "[ERROR] ./extract_calib_constants 실행 파일이 없습니다."
  echo "        먼저 analysis 디렉토리에서 컴파일하세요. (예: ./compile.sh extract_calib_constants.cc)"
  exit 1
fi

echo "============================================================"
echo "[START] extract_calib_constants"
./extract_calib_constants
rc=$?

if [[ ${rc} -ne 0 ]]; then
  echo "[FAIL ] extract_calib_constants (exit=${rc})"
  exit "${rc}"
fi

echo "[DONE ] extract_calib_constants"
exit 0
