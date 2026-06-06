#!/usr/bin/env bash
# =============================================================================
# split_pr.sh
#
# Splits the large working→sandialabs/oso-airfoils PR into ~17 smaller PRs,
# each targeting a logical subset of the changed files.
#
# Usage:
#   ./split_pr.sh              # Run for real
#   ./split_pr.sh --dry-run    # Preview file counts without making any changes
#
# Prerequisites:
#   - Must be run from the repo root with the 'working' branch checked out
#   - git
#   - gh (GitHub CLI, authenticated; needs write access to codykarcher/oso-airfoils
#          and permission to open PRs against sandialabs/oso-airfoils)
#
# What this script does:
#   1. Adds 'sandialabs' as a remote (if not present) and fetches it.
#   2. For each PR group, creates a new branch off sandialabs/main.
#   3. Checks out only the relevant paths from the 'working' branch.
#   4. Commits, pushes to origin (codykarcher/oso-airfoils), and opens a PR.
#
# All branches target sandialabs/main independently — they do NOT chain.
# This means each PR can be reviewed and merged in any order without conflicts,
# since the changed file sets are disjoint.
#
# After all split PRs are merged, close the original large PR #10 without merging.
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Configuration — edit these if your setup differs
# ---------------------------------------------------------------------------
UPSTREAM_REMOTE="sandialabs"
UPSTREAM_URL="https://github.com/sandialabs/oso-airfoils.git"
ORIGIN="origin"
WORKING="working"
BASE="${UPSTREAM_REMOTE}/main"
TARGET_REPO="sandialabs/oso-airfoils"
FORK_USER="codykarcher"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
DRY_RUN=false
[[ "${1:-}" == "--dry-run" ]] && DRY_RUN=true

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
log()  { echo ""; echo ">>> $*"; }
info() { echo "    $*"; }
die()  { echo ""; echo "ERROR: $*" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Preflight checks
# ---------------------------------------------------------------------------
log "Preflight checks"
command -v git >/dev/null 2>&1 || die "git not found in PATH"
command -v gh  >/dev/null 2>&1 || die "gh not found in PATH"

CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)
[[ "${CURRENT_BRANCH}" == "${WORKING}" ]] || \
    die "Must be on the '${WORKING}' branch (currently on '${CURRENT_BRANCH}')."

if $DRY_RUN; then
    info "DRY RUN MODE — no branches will be created, no commits made, no PRs opened."
fi

# ---------------------------------------------------------------------------
# Set up upstream remote and fetch
# (Always fetches even in dry-run, so file counts are accurate)
# ---------------------------------------------------------------------------
if ! git remote | grep -q "^${UPSTREAM_REMOTE}$"; then
    log "Adding remote '${UPSTREAM_REMOTE}' → ${UPSTREAM_URL}"
    git remote add "${UPSTREAM_REMOTE}" "${UPSTREAM_URL}"
fi

log "Fetching ${UPSTREAM_REMOTE}..."
git fetch "${UPSTREAM_REMOTE}"
info "Base commit: $(git rev-parse --short ${BASE})"

# ---------------------------------------------------------------------------
# make_pr BRANCH TITLE BODY PATH [PATH ...]
#
# Creates one PR branch containing the given paths from the working branch.
# Handles both additions/modifications (git checkout) and deletions (git rm).
# ---------------------------------------------------------------------------
make_pr() {
    local branch="$1"
    local title="$2"
    local body="$3"
    shift 3
    local paths=("$@")

    log "PR: ${branch}"
    info "Title: ${title}"
    info "Paths: ${paths[*]}"

    # In dry-run mode: count affected files and return early
    if $DRY_RUN; then
        local total=0
        for p in "${paths[@]}"; do
            local n
            n=$(git diff --name-only "${BASE}...${WORKING}" -- "${p}" 2>/dev/null | wc -l | tr -d ' ')
            total=$((total + n))
        done
        info "[dry-run] ~${total} files would be in this PR."
        return 0
    fi

    # Create the branch off the upstream base
    git checkout -b "${branch}" "${BASE}"

    # Stage additions and modifications from the working branch
    for p in "${paths[@]}"; do
        git checkout "${WORKING}" -- "${p}" 2>/dev/null || true
    done

    # Stage deletions: files that exist in BASE but were removed in WORKING
    for p in "${paths[@]}"; do
        while IFS= read -r deleted_file; do
            [[ -z "${deleted_file}" ]] && continue
            info "Staging deletion: ${deleted_file}"
            git rm -f --ignore-unmatch "${deleted_file}" >/dev/null 2>&1 || true
        done < <(git diff --name-only --diff-filter=D "${BASE}..${WORKING}" -- "${p}" 2>/dev/null)
    done

    git add .

    # Skip this PR if there's nothing new vs the upstream base
    if git diff --cached --quiet; then
        info "No net changes vs ${BASE} — skipping this branch."
        git checkout "${WORKING}"
        git branch -d "${branch}"
        return 0
    fi

    local count
    count=$(git diff --cached --name-only | wc -l | tr -d ' ')
    info "Committing ${count} changed files..."
    git commit -m "${title}"

    info "Pushing '${branch}' to ${ORIGIN}..."
    git push "${ORIGIN}" "${branch}"

    info "Opening PR against ${TARGET_REPO}..."
    gh pr create \
        --repo "${TARGET_REPO}" \
        --base main \
        --head "${FORK_USER}:${branch}" \
        --title "${title}" \
        --body "${body}"

    # Return to working branch before creating the next PR branch
    git checkout "${WORKING}"
    info "Done — PR opened for '${branch}'."
}


# ===========================================================================
# PR DEFINITIONS
#
# Approximate file counts (from git diff --name-only sandialabs/main...working):
#
#   postprocessing/cases/  : ~36,000 files  (GA population run data)
#     cases_71_to_80/      :   1,147   → PR 04
#     cases_101_to_110/
#       case_105 + 106     :   1,955   → PR 05
#       case_107           :   4,798   → PR 06a (~3,180) + 06b (~1,618)
#       case_108           :   1,486   → PR 07
#       case_109           :   4,626   → PR 08a (~2,001) + 08b (~2,625)
#       case_110           :   8,786   → PR 09 (~3,097) + 10 (~3,528) + 11 (~1,523)
#     cases_111_to_120/
#       case_111           :   ~3,922  → PR 12
#       case_112           :   ~4,661  → PR 13
#       case_113           :   ~1,977  → PR 14
#       case_114           :   ~5,254  → PR 15a (~2,465) + 15b (~2,789)
#       case_115           :   ~1,006  → PR 16
#   released_designs/      :   1,342   → PR 17
#   runfiles/              :     466   → PR 03
#   historical_airfoils/   :      17   → PR 01
#   oso-logo/              :       3   → PR 01
#   README.md              :       1   → PR 01
#
# All PRs target sandialabs/main independently.  Max per PR: ~3,528 files.
# ===========================================================================


# ---------------------------------------------------------------------------
# PR 01 — Documentation, logos, historical airfoils         (~21 files)
# ---------------------------------------------------------------------------
make_pr \
    "pr/01-docs-logos-historical" \
    "Add docs, OSO logo, historical MHKF1 airfoil data, and split_pr script" \
    "Updates to top-level README, OSO logo files (ipe/png/svg), historical MHKF1 airfoil digitization data with processing scripts, and the split_pr.sh utility script." \
    "README.md" \
    "oso-logo/" \
    "historical_airfoils/" \
    "split_pr.sh"


# ---------------------------------------------------------------------------
# PR 02 — Postprocessing scripts, config, and caselog       (~16 files)
# ---------------------------------------------------------------------------
make_pr \
    "pr/02-postprocessing-scripts" \
    "Update postprocessing scripts, configuration, and caselog" \
    "Updates to postprocessing Python scripts, JSON pareto config, README, cached aerodynamic data (cp_clean_*.json), caselog, and empty cases_100 placeholder." \
    "postprocessing/README.md" \
    "postprocessing/compare_airfoils_og.py" \
    "postprocessing/compare_airfoils.py" \
    "postprocessing/find_shape_functions.py" \
    "postprocessing/generate_gif.py" \
    "postprocessing/geometry_functions.py" \
    "postprocessing/kulfan.py" \
    "postprocessing/neuralfoil_wrapper_noprint.py" \
    "postprocessing/pareto_plot_config.json" \
    "postprocessing/plot_paretos.py" \
    "postprocessing/postprocess.py" \
    "postprocessing/rainbow_plot_with_comparison.py" \
    "postprocessing/wt_objective_nsga2.py" \
    "postprocessing/xfoil_wrapper_noprint.py" \
    "postprocessing/cached_data/" \
    "postprocessing/pareto_plot.png" \
    "postprocessing/cases/caselog.txt" \
    "postprocessing/cases/cases_100/"


# ---------------------------------------------------------------------------
# PR 03 — Runfiles: run scripts and configurations          (~466 files)
# ---------------------------------------------------------------------------
make_pr \
    "pr/03-runfiles" \
    "Update runfiles and optimization case run scripts" \
    "Updates to runfiles: per-case optimization run scripts (c*.py), JSON run configurations, shell runners, updated GA scripts, revised objective function (wt_objective_nsga2.py replacing wt_objective.py), and updated xfoil/neuralfoil wrappers." \
    "runfiles/"


# ---------------------------------------------------------------------------
# PR 04 — Case data: cases 71–80                          (~1,147 files)
# ---------------------------------------------------------------------------
make_pr \
    "pr/04-cases-71-80" \
    "Add case data: cases 71–80 (case_73 runs)" \
    "GA population snapshot files for optimization case 73 (cases 71–80 group)." \
    "postprocessing/cases/cases_71_to_80/"


# ---------------------------------------------------------------------------
# PR 05 — Case data: cases 105 & 106                      (~1,955 files)
# ---------------------------------------------------------------------------
make_pr \
    "pr/05-cases-105-106" \
    "Add case data: cases 105 and 106" \
    "GA population snapshot files for optimization cases 105 and 106." \
    "postprocessing/cases/cases_101_to_110/case_105/" \
    "postprocessing/cases/cases_101_to_110/case_106/"


# ---------------------------------------------------------------------------
# PR 06a/06b — Case data: case 107  (split 2 ways, ~3,180 / ~1,618 files)
#
# Group A: all t21 runs (4 run dirs, including n564 early run and n752 full runs)
# Group B: remaining thickness sweeps t24–t36
# ---------------------------------------------------------------------------
_C107="postprocessing/cases/cases_101_to_110/case_107"

make_pr \
    "pr/06a-case-107-t21" \
    "Add case data: case 107 (t21 runs)" \
    "GA population data for case 107 thickness t=21% runs (all run directories, ~3,180 files)." \
    "${_C107}/c107_t21_l15_k16_g2000_n564__2025_12_15_12-40/" \
    "${_C107}/c107_t21_l15_k16_g2000_n752__2025_12_16_10-15/" \
    "${_C107}/c107_t21_l15_k16_g2000_n752__2025_12_16_10-24/" \
    "${_C107}/c107_t21_l15_k16_g2000_n752__2025_12_16_17-11/"

make_pr \
    "pr/06b-case-107-t24-t36" \
    "Add case data: case 107 (t24–t36 runs)" \
    "GA population data for case 107 thickness sweep runs t=24% through t=36% (~1,618 files)." \
    "${_C107}/c107_t24_l14_k16_g2000_n752__2025_12_20_01-39/" \
    "${_C107}/c107_t27_l13_k16_g2000_n752__2025_12_21_21-30/" \
    "${_C107}/c107_t30_l12_k16_g2000_n752__2025_12_19_00-28/" \
    "${_C107}/c107_t33_l12_k16_g2000_n752__2025_12_22_20-35/" \
    "${_C107}/c107_t36_l12_k16_g2000_n752__2025_12_20_16-46/"


# ---------------------------------------------------------------------------
# PR 07 — Case data: case 108                             (~1,486 files)
# ---------------------------------------------------------------------------
make_pr \
    "pr/07-case-108" \
    "Add case data: case 108" \
    "GA population snapshot files for optimization case 108." \
    "postprocessing/cases/cases_101_to_110/case_108/"


# ---------------------------------------------------------------------------
# PR 08a/08b — Case data: case 109  (split 2 ways, ~2,001 / ~2,625 files)
#
# Group A: t18 run (largest single run dir at 2,001 files)
# Group B: t21, t24, t27 runs
# ---------------------------------------------------------------------------
_C109="postprocessing/cases/cases_101_to_110/case_109"

make_pr \
    "pr/08a-case-109-t18" \
    "Add case data: case 109 (t18 run)" \
    "GA population data for case 109 thickness t=18% run (~2,001 files)." \
    "${_C109}/c109_t18_l15_k16_g2000_n752__2025_12_31_16-57/"

make_pr \
    "pr/08b-case-109-t21-t27" \
    "Add case data: case 109 (t21–t27 runs)" \
    "GA population data for case 109 thickness sweep runs t=21%, t=24%, and t=27% (~2,625 files)." \
    "${_C109}/c109_t21_l15_k16_g2000_n752__2026_01_07_12-47/" \
    "${_C109}/c109_t24_l14_k16_g2000_n752__2026_01_10_14-28/" \
    "${_C109}/c109_t27_l13_k16_g2000_n752__2026_01_14_19-33/"


# ---------------------------------------------------------------------------
# PR 09–11 — Case data: case 110  (split 3 ways, ~3,097 / ~3,528 / ~1,523)
#
# case_110 has 8,786 files across 5 thickness-sweep run directories.
# Split by thickness value to keep each PR under ~3,500 files.
# Run directories confirmed from: git diff --name-only BASE...working -- case_110/
# ---------------------------------------------------------------------------
_C110="postprocessing/cases/cases_101_to_110/case_110"

make_pr \
    "pr/09-case-110a" \
    "Add case data: case 110 runs t21 and t24" \
    "GA population data for case 110 thickness sweep runs t=21% and t=24% (~3,097 files)." \
    "${_C110}/c110_t21_l15_k16_g2000_n752_m14_p14__2026_01_17_08-45/" \
    "${_C110}/c110_t24_l14_k16_g2000_n752_m14_p14__2026_01_22_11-18/"

make_pr \
    "pr/10-case-110b" \
    "Add case data: case 110 runs t27 and t30" \
    "GA population data for case 110 thickness sweep runs t=27% and t=30% (~3,528 files)." \
    "${_C110}/c110_t27_l13_k16_g2000_n752_m14_p14__2026_02_21_12-43/" \
    "${_C110}/c110_t30_l12_k16_g2000_n752_m14_p14__2026_02_27_13-34/"

make_pr \
    "pr/11-case-110c" \
    "Add case data: case 110 run t33" \
    "GA population data for case 110 thickness sweep run t=33% (~1,523 files)." \
    "${_C110}/c110_t33_l12_k16_g2000_n752_m14_p14__2026_03_05_10-37/"


# ---------------------------------------------------------------------------
# PR 12 — Case data: case 111                             (~3,922 files)
# ---------------------------------------------------------------------------
make_pr \
    "pr/12-case-111" \
    "Add case data: case 111" \
    "GA population snapshot files for optimization case 111." \
    "postprocessing/cases/cases_111_to_120/case_111/"


# ---------------------------------------------------------------------------
# PR 13 — Case data: case 112                             (~4,661 files)
# ---------------------------------------------------------------------------
make_pr \
    "pr/13-case-112" \
    "Add case data: case 112" \
    "GA population snapshot files for optimization case 112." \
    "postprocessing/cases/cases_111_to_120/case_112/"


# ---------------------------------------------------------------------------
# PR 14 — Case data: case 113                             (~1,977 files)
# ---------------------------------------------------------------------------
make_pr \
    "pr/14-case-113" \
    "Add case data: case 113" \
    "GA population snapshot files for optimization case 113." \
    "postprocessing/cases/cases_111_to_120/case_113/"


# ---------------------------------------------------------------------------
# PR 15a/15b — Case data: case 114  (split 2 ways, ~2,465 / ~2,789 files)
#
# Group A: t18, t21, t24, t27 thickness runs
# Group B: t30, t33, t36 runs (including _bad/_meh/_no variants) + t18x exploratory runs
# ---------------------------------------------------------------------------
_C114="postprocessing/cases/cases_111_to_120/case_114"

make_pr \
    "pr/15a-case-114-t18-t27" \
    "Add case data: case 114 (t18–t27 runs)" \
    "GA population data for case 114 thickness runs t=18% through t=27% (~2,465 files)." \
    "${_C114}/c114_t18_k16_n752_l13_e15__2026_05_13_18-25-5282/" \
    "${_C114}/c114_t21_k16_n752_l13_e15__2026_05_14_02-12-5450/" \
    "${_C114}/c114_t24_k16_n752_l13_e15__2026_05_14_09-51-5547/" \
    "${_C114}/c114_t27_k16_n752_l13_e15__2026_05_14_17-37-2982/" \
    "${_C114}/c114_t27_k16_n752_l13_e15__2026_05_14_17-37-2982_collapsed/"

make_pr \
    "pr/15b-case-114-t30-t36-t18x" \
    "Add case data: case 114 (t30–t36 and t18x exploratory runs)" \
    "GA population data for case 114 t=30% through t=36% runs (including intermediate variants) and t18x exploratory runs (~2,789 files)." \
    "${_C114}/c114_t30_k16_n752_l13_e15__2026_05_15_01-36-3995_bad/" \
    "${_C114}/c114_t30_k16_n752_l13_e15__2026_05_15_11-40-3530_no/" \
    "${_C114}/c114_t30_k16_n752_l13_e15__2026_05_15_12-49-2085_meh/" \
    "${_C114}/c114_t30_k16_n752_l13_e15__2026_05_16_13-50-1665/" \
    "${_C114}/c114_t33_k16_n752_l13_e15__2026_05_15_09-20-2590_bad/" \
    "${_C114}/c114_t33_k16_n752_l13_e15__2026_05_15_22-44-0428/" \
    "${_C114}/c114_t36_k16_n752_l13_e15__2026_05_16_06-50-1929/" \
    "${_C114}/c114x_t18_k16_n752_l13_e15__2026_05_13_16-24-2135/" \
    "${_C114}/c114x_t18_k16_n752_l13_e15__2026_05_13_18-21-0586/"


# ---------------------------------------------------------------------------
# PR 16 — Case data: case 115                             (~1,006 files)
# ---------------------------------------------------------------------------
make_pr \
    "pr/16-case-115" \
    "Add case data: case 115" \
    "GA population snapshot files for optimization case 115." \
    "postprocessing/cases/cases_111_to_120/case_115/"


# ---------------------------------------------------------------------------
# PR 17 — Released design packages                        (~1,342 files)
# ---------------------------------------------------------------------------
make_pr \
    "pr/17-released-designs" \
    "Add released design packages (OSO_2025_WT2, OSO_2026_WT2S, OSO_2026_WT3)" \
    "Released design packages: OSO_2025_WT2, OSO_2026_WT2S (WT2 family, based on case 109), and OSO_2026_WT3." \
    "released_designs/"


# ---------------------------------------------------------------------------
# Done
# ---------------------------------------------------------------------------
log "All PRs processed."
echo ""
echo "    View open PRs : https://github.com/${TARGET_REPO}/pulls"
echo ""
echo "    Next steps:"
echo "      1. Wait for each split PR to be reviewed and approved."
echo "      2. Merge them into sandialabs/main one at a time (any order is fine"
echo "         since the file sets are disjoint and each branch targets the same base)."
echo "      3. Once all split PRs are merged, CLOSE (do not merge) the original"
echo "         large PR #10: https://github.com/${TARGET_REPO}/pull/10"
echo ""
