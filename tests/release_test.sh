#!/usr/bin/env bash
#
# Release gate for CanSNPer2 - the scripted half of MANUAL_TEST_COMMANDS.txt.
# Run it from anywhere, it finds the project root itself.
#
#   bash tests/release_test.sh          # full gate, all real alignments (~30 min)
#   bash tests/release_test.sh fast     # no alignments: unit tests, database, exports,
#                                       # failure paths                     (~2 min)
#
# It prints one line per check and a summary, and exits 1 if anything failed. A check that
# cannot run here (no network, no progressiveMauve) is reported as WARN, not FAIL, so the
# gate stays honest on machines that are not the reference one.
#
# Data it needs from the project (issue #185, the data lives next to the repository):
#   ../CanSNPer2-data/source/francisella_tularensis_{tree,snp,references}.txt
#   references/*.fna in this repository, and ../../test_data/*.fasta for the known answers
# Override the root with ROOT=/path/to/project when the layout is different.

set -u

TIER=${1:-full}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/.." && pwd)
ROOT=${ROOT:-$(cd "$REPO/.." && pwd)}
PKG=$REPO
SRC=$ROOT/CanSNPer2-data/source
PY=${PYTHON:-python}
OUT=${OUT:-$ROOT/test-run/release_gate}

export PYTHONPATH=$PKG
rm -rf "$OUT"
mkdir -p "$OUT/logs"

LASTLOG=""
PASSED=0
FAILED=0
WARNED=0
STEP=0

cmd() { # cmd <name> <expect-ok|expect-fail> <grep-for> <command...>
    local name=$1 want=$2 needle=$3; shift 3
    STEP=$((STEP+1))
    local log="$OUT/logs/$(printf '%02d' "$STEP")-$name.log"
    LASTLOG=$log
    local rc=0
    "$@" > "$log" 2>&1 || rc=$?
    LASTRC=$rc
    if [ "$want" = any ]; then
        printf '  run   %-52s exit %s, %s\n' "$name" "$rc" "$(basename "$log")"; return
    fi
    if [ "$want" = ok ] && [ $rc -ne 0 ]; then
        printf '  FAIL  %-52s exit %s, see %s\n' "$name" "$rc" "$(basename "$log")"; FAILED=$((FAILED+1)); return
    fi
    if [ "$want" = fail ] && [ $rc -eq 0 ]; then
        printf '  FAIL  %-52s should have failed, exit 0\n' "$name"; FAILED=$((FAILED+1)); return
    fi
    if [ -n "$needle" ] && ! grep -qi -- "$needle" "$log"; then
        printf '  FAIL  %-52s exit %s but no "%s" in %s\n' "$name" "$rc" "$needle" "$(basename "$log")"
        FAILED=$((FAILED+1)); return
    fi
    printf '  ok    %-52s exit %s, %s\n' "$name" "$rc" "$(basename "$log")"; PASSED=$((PASSED+1))
}

expect() { # expect <name> <want> <have>
    if [ "$2" = "$3" ]; then
        printf '  ok    %-52s %s\n' "$1" "$3"; PASSED=$((PASSED+1))
    else
        printf '  FAIL  %-52s wanted [%s], got [%s]\n' "$1" "$2" "$3"; FAILED=$((FAILED+1))
    fi
}

warn() { printf '  warn  %-52s %s\n' "$1" "$2"; WARNED=$((WARNED+1)); }
say()  { printf '\n%s\n' "$1"; }

say "== the package itself"
cmd compileall ok "" "$PY" -m compileall -q "$PKG/CanSNPer2"
cmd unit-tests ok "" "$PY" -m pytest -q "$REPO/tests"
cmd selftest ok "required python packages installed" "$PY" -m CanSNPer2.selftest

say "== the documentation names every command line option"
DOCMISS=0
for tool in CanSNPerTree SNPDatabase DownloadGenomes; do
    for opt in $("$PY" "$PKG/CanSNPer2/$tool.py" --help 2>&1 | grep -o -- '--[a-z_]*' | sort -u); do
        if grep -q -- "$opt" "$REPO/README.md"; then :; else
            printf '  FAIL  %-52s %s is not documented in README.md\n' "$tool $opt" "$opt"
            DOCMISS=$((DOCMISS+1)); FAILED=$((FAILED+1))
        fi
    done
done
if [ "$DOCMISS" = 0 ]; then
    printf '  ok    %-52s all options of the three commands appear in README.md\n' "README matches --help"
    PASSED=$((PASSED+1))
fi

say "== database build and exports"
DB=$OUT/gate.db
cmd build-db ok "" "$PY" "$PKG/CanSNPer2/SNPDatabase.py" --database "$DB" --create \
    --tree "$SRC/francisella_tularensis_tree.txt" --annotation "$SRC/francisella_tularensis_snp.txt" \
    --references "$SRC/francisella_tularensis_references.txt" --source_type CanSNPer \
    --logs "$OUT/logs/db"
COUNTS=$("$PY" -c "import sqlite3;c=sqlite3.connect('$DB');print(' '.join(str(c.execute('SELECT COUNT(*) FROM '+t).fetchone()[0]) for t in ('nodes','tree','rank','snp_annotation','snp_references')))")
expect "row counts nodes/tree/rank/annotation/references" "329 329 2 328 5" "$COUNTS"
cmd export-tab ok "" "$PY" "$PKG/CanSNPer2/SNPDatabase.py" --database "$DB" --export \
    --export_format tab -o "$OUT/export_tab" --logs "$OUT/logs/export"
expect "tab export annotation.dmp" "329" "$(wc -l < "$OUT/export_tab/annotation.dmp")"
expect "tab export snptree.dmp" "330" "$(wc -l < "$OUT/export_tab/snptree.dmp")"
cmd export-newick ok "" "$PY" "$PKG/CanSNPer2/SNPDatabase.py" --database "$DB" --export \
    --export_format newick -o "$OUT/export_nwk" --logs "$OUT/logs/export"
NW=$(ls "$OUT/export_nwk"/*.nwk 2>/dev/null | head -1)
expect "newick starts with the root ($(basename "${NW:-none}"))" "((((M.1," "$(cut -c1-8 "${NW:-/dev/null}")"

say "== the id 0 regression (issues #202, #222)"
cat > "$OUT/id0_check.py" <<'PYEOF'
import logging, os, shutil, sqlite3, sys
logging.basicConfig(level=logging.WARNING)
from CanSNPer2.modules.DatabaseConnection import CreateDatabase, DatabaseFunctions, CanSNPdbFunctions

work = sys.argv[1]
if os.path.isdir(work):
    shutil.rmtree(work)
os.makedirs(work)
db = os.path.join(work, "idzero.db")
CreateDatabase().create_database(database=db)
conn = sqlite3.connect(db)
conn.execute("INSERT INTO rank (rank_i,rank) VALUES (1,'species')")
conn.commit(); conn.close()

fns = DatabaseFunctions(db); fns.add_node("ZERO_NODE", id=0); fns.add_rank("subtaxon", id=0); fns.commit()
dbf = CanSNPdbFunctions(db); dbf.add_node("FIVE_NODE", id=5); dbf.add_rank("subtaxon2", id=7); dbf.commit()

rows = sqlite3.connect(db).execute("SELECT id,name FROM nodes ORDER BY id").fetchall()
ranks = sqlite3.connect(db).execute("SELECT rank_i,rank FROM rank ORDER BY rank_i").fetchall()
print("nodes:", rows)
print("ranks:", ranks)
want_nodes = [(0, "ZERO_NODE"), (5, "FIVE_NODE")]
want_ranks = [(0, "subtaxon"), (7, "subtaxon2")]
missing = [w for w in want_nodes + want_ranks if w not in rows + ranks]
if missing:
    print("MISSING:", missing)
    sys.exit(1)
print("supplied ids kept as given, including id 0")
PYEOF
cmd id-zero ok "supplied ids kept" "$PY" "$OUT/id0_check.py" "$OUT/id0"

mkdir -p "$OUT/refs"
for r in FSC200 FTNF002-00 LVS OSU18 SCHUS4.2; do cp "$PKG/references/$r.fna" "$OUT/refs/"; done
mkdir -p "$OUT/inputs" && cp "$PKG/references/FSC200.fna" "$OUT/inputs/query.fna"

say "== the failure paths must be loud (issues #176, #177, #179, #181, #190, #191)"
mkdir -p "$OUT/empty_refs"
cmd no-references fail "No reference genomes" "$PY" -m CanSNPer2.modules.runMauve \
    "$PKG/references/FSC200.fna" --refdir "$OUT/empty_refs" --tmpdir "$OUT/t_norefs" \
    --logs "$OUT/logs/norefs" --outdir "$OUT/r_norefs"
cmd missing-query-file fail "was not found" "$PY" -m CanSNPer2.modules.runMauve \
    "$OUT/no_such_file.fna" --refdir "$OUT/empty_refs" --tmpdir "$OUT/t_nofile" \
    --logs "$OUT/logs/nofile" --outdir "$OUT/r_nofile"
cmd skip-mauve-without-alignments fail "No alignments were produced" "$PY" -m CanSNPer2.modules.runMauve \
    "$OUT/inputs/query.fna" --refdir "$OUT/refs" --tmpdir "$OUT/t_empty" \
    --logs "$OUT/logs/empty" --outdir "$OUT/r_empty" --skip_mauve

say "== download failure is reported and fails the command (issue #218)"
DLLOG=$OUT/logs/download
DLRC=$("$PY" "$PKG/CanSNPer2/DownloadGenomes.py" -db "$DB" -o "$OUT/dl_refs" --logs "$DLLOG" >/dev/null 2>&1; echo $?)
DERR=$(grep -rh "Could not download" "$DLLOG" 2>/dev/null | wc -l | tr -d ' ')
if [ "$DLRC" = 1 ] && [ "$DERR" -gt 0 ]; then
    printf '  ok    %-52s exit 1, %s failures reported\n' "download with no network" "$DERR"; PASSED=$((PASSED+1))
elif [ "$DLRC" = 0 ]; then
    warn "download with no network" "this machine reached NCBI, the failure path was not exercised"
else
    printf '  FAIL  %-52s exit %s, %s failures reported\n' "download with no network" "$DLRC" "$DERR"; FAILED=$((FAILED+1))
fi

if [ "$TIER" != full ]; then
    say "== summary ($TIER tier)"
    printf '   %s passed, %s failed, %s warnings - logs in %s\n' "$PASSED" "$FAILED" "$WARNED" "$OUT/logs"
    [ $FAILED -eq 0 ] || exit 1
    exit 0
fi

say "== real alignments: the known answers"
cmd align-query ok "5 alignments written" "$PY" -m CanSNPer2.modules.runMauve "$OUT/inputs/query.fna" \
    --refdir "$OUT/refs" --tmpdir "$OUT/tmp_align" --logs "$OUT/logs/align" --outdir "$OUT/align" --keep_temp
expect "five alignments" "5" "$(ls "$OUT/align" | wc -l | tr -d ' ')"
expect "mauve left no sslist in the input folders" "0" "$(find "$OUT/refs" "$OUT/inputs" -name '*.sslist' | wc -l | tr -d ' ')"
cmd type-query ok "query: B.22" "$PY" "$PKG/CanSNPer2/CanSNPerTree.py" "$OUT/inputs/query.fna" -db "$DB" \
    --refdir "$OUT/refs" -o "$OUT/res_query" --tmpdir "$OUT/tmp_type" --logdir "$OUT/logs/type" --save_tree
expect "query_snps.txt lines" "19" "$(wc -l < "$OUT/res_query/query_snps.txt")"
expect "query_not_called.txt lines" "313" "$(wc -l < "$OUT/res_query/query_not_called.txt")"
if [ -s "$OUT/res_query/query_tree.pdf" ]; then
    printf '  ok    %-52s %s bytes\n' "query_tree.pdf written" "$(wc -c < "$OUT/res_query/query_tree.pdf" | tr -d ' ')"
    PASSED=$((PASSED+1))
else
    printf '  FAIL  %-52s missing or empty\n' "query_tree.pdf written"; FAILED=$((FAILED+1))
fi

mkdir -p "$OUT/tmp_skip/CanSNPer2" && cp "$OUT"/align/*.xmfa "$OUT/tmp_skip/CanSNPer2/"
cmd skip-mauve-reuse ok "query: B.22" "$PY" "$PKG/CanSNPer2/CanSNPerTree.py" "$OUT/inputs/query.fna" -db "$DB" \
    --refdir "$OUT/refs" -o "$OUT/res_skip" --tmpdir "$OUT/tmp_skip" --logdir "$OUT/logs/skip" --skip_mauve
expect "same snps through --skip_mauve" "19" "$(wc -l < "$OUT/res_skip/query_snps.txt")"
expect "same not-called through --skip_mauve" "313" "$(wc -l < "$OUT/res_skip/query_not_called.txt")"

if ls "$ROOT"/test_data/*.fasta >/dev/null 2>&1; then
    ls "$ROOT"/test_data/*.fasta > "$OUT/real_list.txt"
    cmd real-data-batch ok "B.75" "$PY" "$PKG/CanSNPer2/CanSNPerTree.py" "$OUT/real_list.txt" -db "$DB" \
        --refdir "$OUT/refs" -o "$OUT/res_real" --tmpdir "$OUT/tmp_real" --logdir "$OUT/logs/real" --save_tree
    for pair in "A271_1_FDC408_complete.fasta B.75 20 312" "FSC458.fasta B.77 17 315"; do
        set -- $pair
        LABEL=$(grep -h "^$(basename $1 .fasta):" "$LASTLOG" | tail -1 | sed 's/.*: //')
        expect "known answer $1" "$2" "$LABEL"
        expect "  snps for $1" "$3" "$(wc -l < "$OUT/res_real/$(basename $1 .fasta)_snps.txt" 2>/dev/null || echo 0)"
        expect "  not called for $1" "$4" "$(wc -l < "$OUT/res_real/$(basename $1 .fasta)_not_called.txt" 2>/dev/null || echo 0)"
    done
else
    warn "known answer batch" "no test_data/*.fasta found in $ROOT"
fi

say "== leave-one-out: every reference typed without its own genome"
for q in FSC200 FTNF002-00 LVS OSU18 SCHUS4.2; do
    mkdir -p "$OUT/loo_refs_$q"
    for r in FSC200 FTNF002-00 LVS OSU18 SCHUS4.2; do
        [ "$r" = "$q" ] || cp "$PKG/references/$r.fna" "$OUT/loo_refs_$q/"
    done
    cmd "loo-$q" any "" "$PY" "$PKG/CanSNPer2/CanSNPerTree.py" "$PKG/references/$q.fna" -db "$DB" \
        --refdir "$OUT/loo_refs_$q" -o "$OUT/res_loo_$q" --tmpdir "$OUT/tmp_loo_$q" \
        --logdir "$OUT/logs/loo_$q"
    SUB=$(grep -h "$q:" "$LASTLOG" | tail -1 | sed "s/.*: //")
    SNPS=$(wc -l < "$OUT/res_loo_$q/${q}_snps.txt" 2>/dev/null | tr -d ' ')
    echo "leave-one-out $q -> ${SUB:-no result} (snps ${SNPS:-0})" >> "$OUT/leave_one_out_baseline.txt"
    if [ -n "$SUB" ] && [ "$SUB" != NA ] && [ "${SNPS:-0}" -gt 0 ]; then
        printf '  ok    %-52s %s, %s snps\n' "leave-one-out $q" "$SUB" "$SNPS"; PASSED=$((PASSED+1))
    elif [ "$SUB" = NA ]; then
        if grep -q "No sample was typed" "$LASTLOG"; then
            warn "leave-one-out $q" "NA, said so loudly (see $OUT/leave_one_out_baseline.txt)"
        else
            printf '  FAIL  %-52s NA without any error message\n' "leave-one-out $q"; FAILED=$((FAILED+1))
        fi
    else
        printf '  FAIL  %-52s no subtaxon, no result files\n' "leave-one-out $q"; FAILED=$((FAILED+1))
    fi
done
echo "the leave-one-out pairs above are the regression baseline, they belong in"
echo "MANUAL_TEST_COMMANDS.txt step 7 - they were recorded in $OUT/leave_one_out_baseline.txt"

say "== summary ($TIER tier)"
printf '   %s passed, %s failed, %s warnings - logs in %s\n' "$PASSED" "$FAILED" "$WARNED" "$OUT/logs"
[ $FAILED -eq 0 ] || exit 1
