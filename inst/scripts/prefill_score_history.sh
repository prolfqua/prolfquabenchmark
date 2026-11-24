#!/usr/bin/env bash
# Prefills inst/extdata/score_history.tsv with older prolfqua versions (done once, 2026-11-23):
#   1.6.2  2026-05-28  after the first AI-assisted wave
#   1.9.0  2026-11-02  before the rOpenSci self-review
# Each version is installed into its own library together with the prolfquapp, prolfquatakeoutmissing and
# prolfquabenchmark of its time (CRAN and Bioconductor packages come from the normal library). The prolfqua
# benchmark vignettes of that time are rendered, the core prolfqua steps are timed with today's
# time_prolfqua_steps.R, and the results are appended with that version and its commit date.
#
# Record of a run: it called inst/scripts/update_score_history.R, timing_input.R and time_prolfqua_steps.R as of
# commit ae5316f (the last two now live in workflow/scripts, the first was replaced by the pipeline's score_history
# rule); inst/scripts/migrate_score_history.R then relabelled its rows to the pipeline's names. The run also reran
# prolfqua 1.4.0 (2025-11-09, with render_harvest.R); its rows were dropped on 2026-11-24 because its models and
# vignettes share too little with today's.
#
# Run from the prolfquabenchmark root, next to ../prolfqua, ../prolfquapp and ../prolfquatakeoutmissing:
#   inst/scripts/prefill_score_history.sh render <work_dir>    # install and render; the machine may be busy
#   inst/scripts/prefill_score_history.sh collect <work_dir>   # time the steps and append; needs a quiet machine
set -euo pipefail

STEP=$1
mkdir -p "$2"
WORK=$(cd "$2" && pwd)
ROOT=$(pwd)
SCRIPTS=$ROOT/inst/scripts

# name prolfqua prolfquapp prolfquatakeoutmissing(- if none) prolfquabenchmark
POINTS=(
  "1.6.2 f4dfb793 2972ef65 - 09c3cd2"
  "1.9.0 f4dfb7f3 18914707 7811e5a 9f3665e"
)

render_point() {
  local name=$1 prolfqua=$2 app=$3 missing=$4 bench=$5
  local dir=$WORK/$name lib=$WORK/$name/lib
  rm -rf "$dir"
  mkdir -p "$lib"
  extract() { mkdir -p "$dir/$1" && git -C "$ROOT/../$1" archive "$2" | tar -x -C "$dir/$1"; }
  extract prolfqua "$prolfqua"
  extract prolfquapp "$app"
  [ "$missing" != "-" ] && extract prolfquatakeoutmissing "$missing"
  extract prolfquabenchmark "$bench"
  local packages="prolfqua prolfquapp prolfquatakeoutmissing prolfquabenchmark"
  for package in $packages; do
    if [ -d "$dir/$package" ]; then
      R_LIBS=$lib R CMD INSTALL -l "$lib" "$dir/$package" > "$dir/install_$package.log" 2>&1
    fi
  done
  local results=$dir/prolfquabenchmark/inst/Benchresults
  mkdir -p "$results"
  touch "$results/.score_run"
  cd "$dir/prolfquabenchmark/vignettes"
  case $name in
    1.6.2)
      R_LIBS=$lib Rscript -e "rmarkdown::render('Benchmark_prolfqua.Rmd', quiet = TRUE)"
      ;;
    1.9.0)
      R_LIBS=$lib Rscript -e "rmarkdown::render('Benchmark_prolfqua.Rmd', quiet = TRUE)"
      # this version kept the hurdle results in tempdir(); keep them like today's vignette does
      sed -i '' 's|file.path(tempdir(), paste0("hurdle_IonStar_MQ_"|file.path("../inst/Benchresults", paste0("hurdle_IonStar_MQ_"|' \
        Benchmark_hurdle.Rmd
      R_LIBS=$lib Rscript -e "rmarkdown::render('Benchmark_hurdle.Rmd', quiet = TRUE)"
      ;;
  esac
  cd "$ROOT"
}

collect_point() {
  local name=$1 prolfqua=$2
  local dir=$WORK/$name lib=$WORK/$name/lib
  local results=$dir/prolfquabenchmark/inst/Benchresults
  R_LIBS=$lib Rscript "$SCRIPTS/time_prolfqua_steps.R" "$WORK/timing_input.rds" "$results/timing"
  Rscript "$SCRIPTS/update_score_history.R" --results-dir="$results" --version="$name" \
    --date="$(git -C "$ROOT/../prolfqua" log -1 --format=%ad --date=short "$prolfqua")"
}

case $STEP in
  render)
    for point in "${POINTS[@]}"; do
      echo "== render $point"
      render_point $point
    done
    ;;
  collect)
    Rscript "$SCRIPTS/timing_input.R" "$WORK/timing_input.rds"
    for point in "${POINTS[@]}"; do
      echo "== collect $point"
      collect_point $point
    done
    ;;
esac
