# =====================================================================
#  40_staging.smk — D-05, the manual steps made declarative
#
#  Three rules that between them replace two undocumented manual `cp`
#  campaigns and add one verification gate that never existed.
#
#  ------------------------------------------------------------------
#  1. `stage_coolbox_bigwigs` — SYMLINKS, not copies
#  ------------------------------------------------------------------
#  `/mnt/scratch/nikolai/Antonia/bigwigs_for_coolbox` (165 files, 47 GB)
#  is byte-identical to the `merged/` (20 bp) tracks it was copied from
#  [the ground-truth inventory §7, verified: both 209,327,288 B for
#  H3K27me3_E6_WT_Xi.bw]. No script in the original repository creates
#  it. We recreate it as symlinks, which is semantically identical for
#  every consumer (`01_08_coolbox_diff_hic_only.py` only opens the file)
#  and reclaims ~47 GB.
#
#  The inventory is reproduced exactly, INCLUDING its one deliberate
#  omission: 33 AcMe3 + 27 CTCF + 33 H3K27ac + 33 H3K27me3 + 39 RNA-Seq
#  = 165, and **no Rad21**. Rad21 exists as a merged20 track (6 of them,
#  B1621 only) and is deliberately not staged, because coolbox never
#  draws it. The count is asserted below, so a sample-sheet edit that
#  silently changes the browser panel set is visible at DAG build
#  instead of at figure time.
#
#  ------------------------------------------------------------------
#  2. `stage_coolbox_boundaries` — symlinks + a step from NO original
#  ------------------------------------------------------------------
#  `boundaries_for_coolbox` (108 files) holds `*_valleys.bed` / `.bw`
#  copies PLUS `*_valleys.bed.bgz` and `.bed.bgz.tbi`. **Nothing in the
#  original repository produces the .bgz/.tbi pair** — no script, no
#  notebook. It was done by hand with bgzip/tabix and is reconstructed
#  here from the file format alone (BGZF + a tabix index of a 0-based
#  BED). Recording that explicitly is the point of the rule: it is the
#  one place in the pipeline reproducing an artefact whose command line
#  was never written down.
#
#  Only the **Xa / Xi** stems are indexed (22 of them). The 11 `Gall`
#  stems get `.bed` + `.bw` and no index, exactly as on disk — coolbox
#  only ever range-queries the allele-split tracks.
#
#  Expected counts here: 33 bed + 33 bw + 22 bgz + 22 tbi = 110, against
#  the ground truth's 32 + 32 + 22 + 22 = 108. The +2 is one extra
#  track, `H3K27me3_F3_CTCF-NodTAG_Gall`, which is the known "18 degron
#  valley BEDs where ground truth has 17" deviation endorsed in
#  the design review §8. It is not a defect in this module.
#
#  ------------------------------------------------------------------
#  3. `check_fixtures` — the gate on the human-judgement inputs
#  ------------------------------------------------------------------
#  Nine irreducible human inputs live in `resources/fixtures/` (§7.1,
#  the archaeology notes §J.7 / A-8). If one is edited, every downstream number moves
#  and no other rule would notice. `sha256sum -c` against the manifest
#  runs before anything consumes them.
#
#  NO FIGURES HERE. This module emits symlinks, index files and a
#  verification stamp.
# =====================================================================

STAGING_EXCLUDE_MARKS = ("Rad21",)      # deliberate — see header
STAGING_EXPECTED_BIGWIGS = 165          # the ground-truth inventory §7
STAGING_INDEXED_ALLELES = ("Xa", "Xi")  # Gall is staged but never indexed


def _staged_bigwig_tracks():
    """The 165-file coolbox bigWig inventory: every merged20 track except
    Rad21, plus every AcMe3 track."""
    tracks = [
        t for t in SS.tracks
        if naming.parse_track(t)["mark"] not in STAGING_EXCLUDE_MARKS
    ]
    tracks += SS.acme3_tracks()
    return sorted(tracks)


def _staged_boundary_tracks():
    """Valley calls are H3K27me3-only, and coolbox reads the chrX set."""
    return sorted(SS.tracks_of(mark="H3K27me3"))


def _staged_indexed_tracks():
    """Xa / Xi only. `Gall` valley BEDs are staged but never tabix-indexed."""
    return [
        t for t in _staged_boundary_tracks()
        if naming.parse_track(t)["allele"] in STAGING_INDEXED_ALLELES
    ]


def _staged_bigwig_source(wildcards):
    """AcMe3 lives in its own tree; every other mark under merged20/{mark}/."""
    mark = wildcards.track.split("_", 1)[0]
    if mark == "AcMe3":
        return P.acme3_bw20(wildcards.track)
    return P.merged20(mark, wildcards.track)


STAGED_BIGWIG_TRACKS = _staged_bigwig_tracks()
STAGED_BOUNDARY_TRACKS = _staged_boundary_tracks()
STAGED_INDEXED_TRACKS = _staged_indexed_tracks()

if len(STAGED_BIGWIG_TRACKS) != STAGING_EXPECTED_BIGWIGS:
    # Not fatal: a user running their own clones legitimately has a different
    # count. But the paper run must be 165, so say so loudly rather than
    # discovering it in a browser panel three hours later.
    sys.stderr.write(
        f"NOTE  coolbox bigWig staging set is {len(STAGED_BIGWIG_TRACKS)} tracks, "
        f"not the paper's {STAGING_EXPECTED_BIGWIGS}. Expected only when the "
        f"sample sheet differs from the published one.\n"
    )


# ---------------------------------------------------------------------
# coolbox staging
# ---------------------------------------------------------------------
rule stage_coolbox_bigwigs:
    """Symlink a 20 bp merged bigWig into the coolbox staging tree.

    `ln -sfn` to the *resolved* source, so the link survives the source
    being itself a symlink and so `realpath` on the staged file lands on
    the real bigWig -- that identity is how the testing team checks this
    rule (realpath equality, never content comparison).
    """
    input:
        _staged_bigwig_source,
    output:
        P.staged_bigwig("{track}"),
    log:
        P.log("stage_coolbox_bigwigs", "{track}"),
    resources:
        mem_mb=1000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        ln -sfn "$(readlink -f {input})" {output} 2> {log}
        echo "staged {output} -> $(readlink -f {input})" >> {log}
        """


rule stage_coolbox_boundaries:
    """Symlink a chrX valley BED + bigWig into the coolbox staging tree.

    `01_08_coolbox_diff_hic_only.py` derives the BED path from the bigWig
    path by stem substitution, so both must land side by side under the
    same stem. Keeping them in one rule makes that pairing structural
    rather than incidental.
    """
    input:
        bed=P.valleys("chrX", "{track}"),
        bw=P.valleys("chrX", "{track}", ext="bw"),
    output:
        bed=P.staged_boundary("{track}", ext="bed"),
        bw=P.staged_boundary("{track}", ext="bw"),
    log:
        P.log("stage_coolbox_boundaries", "{track}"),
    resources:
        mem_mb=1000,
    shell:
        r"""
        mkdir -p "$(dirname {output.bed})" "$(dirname {log})"
        ln -sfn "$(readlink -f {input.bed})" {output.bed} 2> {log}
        ln -sfn "$(readlink -f {input.bw})"  {output.bw} 2>> {log}
        echo "staged {output.bed} -> $(readlink -f {input.bed})" >> {log}
        echo "staged {output.bw} -> $(readlink -f {input.bw})" >> {log}
        """


rule index_coolbox_boundaries:
    """bgzip + tabix the staged valley BED. NO ORIGINAL SCRIPT DOES THIS.

    The `.bed.bgz` / `.bed.bgz.tbi` pairs in the ground-truth
    `boundaries_for_coolbox` tree were produced by hand; the command line
    was never recorded. This is the reconstruction, from the file format:
    BGZF-compressed, tabix `-p bed` (0-based, chrom/start/end in columns
    1/2/3). Xa and Xi only -- `Gall` stems are staged unindexed, matching
    the 22-of-33 split on disk.

    The BED is compressed from the *resolved* symlink target rather than
    through the link, so the index offsets refer to a file written by this
    rule and never to a stale staged copy.
    """
    input:
        bed=P.staged_boundary("{track}", ext="bed"),
    output:
        bgz=P.staged_boundary("{track}", ext="bed.bgz"),
        tbi=P.staged_boundary("{track}", ext="bed.bgz.tbi"),
    wildcard_constraints:
        track=r"[^/]+_(Xa|Xi)",
    log:
        P.log("index_coolbox_boundaries", "{track}"),
    resources:
        mem_mb=1000,
    shell:
        r"""
        mkdir -p "$(dirname {output.bgz})" "$(dirname {log})"
        sort -k1,1 -k2,2n "$(readlink -f {input.bed})" \
            | bgzip -c > {output.bgz} 2> {log}
        tabix -f -p bed {output.bgz} 2>> {log}
        """


# ---------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------
rule check_fixtures:
    """Verify every human-judgement input against `fixtures.sha256`.

    Runs from the repository root because the manifest stores repo-relative
    paths. `verify_checksums: false` degrades this to presence-only, which
    is what a user bringing their own curated loops wants; the paper run
    leaves it on.
    """
    input:
        manifest=P.fixture("fixtures.sha256"),
    output:
        touch(P.fixtures_verified()),
    params:
        repo=P.repo,
        verify=bool(config["fixtures"].get("verify_checksums", True)),
    log:
        P.log("check_fixtures"),
    resources:
        mem_mb=1000,
    shell:
        r"""
        mkdir -p "$(dirname {output})" "$(dirname {log})"
        LOG="$(readlink -f {log} 2>/dev/null || echo "$PWD/{log}")"
        cd {params.repo}
        if [ "{params.verify}" = "True" ]; then
            sha256sum -c "{input.manifest}" > "$LOG" 2>&1 \
                || {{ echo "FIXTURE CHECKSUM MISMATCH -- see $LOG" >&2; exit 1; }}
            echo "verified $(wc -l < "{input.manifest}") fixtures" >> "$LOG"
        else
            echo "fixtures.verify_checksums=false: presence-only check" > "$LOG"
            cut -d' ' -f3- "{input.manifest}" | while read -r f; do
                [ -e "$f" ] || {{ echo "missing fixture: $f" >&2; exit 1; }}
            done
        fi
        """


# ---------------------------------------------------------------------
# Convenience target
# ---------------------------------------------------------------------
rule staging:
    """Everything 90_visualise's coolbox panels need staged, plus the
    fixture gate. Not part of `rule all` -- the panel rules pull the
    individual files they use."""
    input:
        P.fixtures_verified(),
        [P.staged_bigwig(t) for t in STAGED_BIGWIG_TRACKS],
        [P.staged_boundary(t, ext=e) for t in STAGED_BOUNDARY_TRACKS for e in ("bed", "bw")],
        [
            P.staged_boundary(t, ext=e)
            for t in STAGED_INDEXED_TRACKS
            for e in ("bed.bgz", "bed.bgz.tbi")
        ],
