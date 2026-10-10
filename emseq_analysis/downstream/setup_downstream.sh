#!/usr/bin/env bash
# One-time setup for the downstream steps (05 regions, 06 tissue-of-origin).
# Needs internet, so run it on the LOGIN node:  bash downstream/setup_downstream.sh
#   1. creates the conda env $DOWNSTREAM_ENV (python + pandas/numpy/scipy/matplotlib,
#      samtools, htslib, bedtools, a C++ compiler)
#   2. clones and compiles wgbstools + UXM (Loyfer et al. 2023) into OUTDIR/tools
#   3. downloads hg38 CpG-island and RefSeq gene annotation from UCSC into OUTDIR/annotation
source "${PIPE_DIR:-$(dirname "$0")/..}/scripts/common.sh"

# ---- 1. conda env ----------------------------------------------------------------
[[ -n "${CONDA_SETUP:-}" ]] && eval "$CONDA_SETUP"
set +u; source "$(conda info --base)/etc/profile.d/conda.sh"; set -u
if conda env list | awk '{print $1}' | grep -qx "$DOWNSTREAM_ENV"; then
    log "conda env $DOWNSTREAM_ENV exists"
else
    log "creating conda env $DOWNSTREAM_ENV (a few minutes)"
    conda create -y -n "$DOWNSTREAM_ENV" -c conda-forge -c bioconda \
        "python=3.11" pandas numpy scipy matplotlib samtools htslib bedtools cxx-compiler
fi
activate_env "$DOWNSTREAM_ENV"

# ---- 2. wgbstools + UXM ------------------------------------------------------------
mkdir -p "$TOOLS_DIR"
if [[ ! -d "$TOOLS_DIR/wgbs_tools" ]]; then
    git clone --depth 1 https://github.com/nloyfer/wgbs_tools.git "$TOOLS_DIR/wgbs_tools"
fi
if [[ ! -x "$TOOLS_DIR/wgbs_tools/src/pipeline_wgbs/patter" ]]; then
    log "compiling wgbstools"
    # GCC >= 13 no longer pulls in <cstdint> implicitly, which wgbstools relies on;
    # wrap g++ so every compile gets it.
    mkdir -p "$TOOLS_DIR/gxx_wrap"
    printf '#!/bin/sh\nexec "%s" -include cstdint "$@"\n' "$(command -v g++)" > "$TOOLS_DIR/gxx_wrap/g++"
    chmod +x "$TOOLS_DIR/gxx_wrap/g++"
    (cd "$TOOLS_DIR/wgbs_tools" && PATH="$TOOLS_DIR/gxx_wrap:$PATH" python setup.py)
fi
[[ -x "$TOOLS_DIR/wgbs_tools/src/pipeline_wgbs/patter" ]] \
    || die "wgbstools did not compile (patter missing); see messages above"
if [[ ! -d "$TOOLS_DIR/UXM_deconv" ]]; then
    git clone --depth 1 https://github.com/nloyfer/UXM_deconv.git "$TOOLS_DIR/UXM_deconv"
fi
[[ -s "$TOOLS_DIR/UXM_deconv/supplemental/$UXM_ATLAS" ]] || die "atlas $UXM_ATLAS not found in UXM_deconv/supplemental"

# ---- 3. annotation -------------------------------------------------------------------
mkdir -p "$ANNOT_DIR"
UCSC=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database
for f in cpgIslandExt.txt.gz ncbiRefSeqCurated.txt.gz; do
    [[ -s "$ANNOT_DIR/$f" ]] || curl -fsSL "$UCSC/$f" -o "$ANNOT_DIR/$f" || die "download of $f failed"
done

log "downstream setup done. Next: bash run_downstream.sh"
