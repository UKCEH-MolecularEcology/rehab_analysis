# Source this before any `pixi` command in this repo (installing envs or
# running `pixi r -e <env> ...` from a rule's shell block). Keeps pixi's
# binary, package cache, and temp files project-local under tools/ and
# tmp/ -- never $HOME, which has a tight disk quota on this host (see
# README's "Tools" section for the earlier quota incident this avoids).
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PIXI_HOME="${REPO_DIR}/tools/pixi_home"
export PIXI_CACHE_DIR="${REPO_DIR}/tools/pixi_home/cache"
export RATTLER_CACHE_DIR="${REPO_DIR}/tools/pixi_home/cache/rattler"
export TMPDIR="${REPO_DIR}/tmp"
export PATH="${REPO_DIR}/tools/pixi_home/bin:$PATH"
