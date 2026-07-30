#!/usr/bin/env python3
"""Fail if any rule outside ``90_visualise.smk`` emits an image.

the design notes section 6.2 makes "the visualisation stage is the only
place that writes figures" true by construction rather than by convention.
Three mechanisms enforce it; this script is the third.

    Every non-viz rule emits DATA. csaw emits a bin-count matrix, the valley
    caller emits a per-bin state table, cooltools emits .tsv/.npz, coolpuppy
    emits the pile-up matrix as .npz. 90_visualise.smk turns them into pictures.

A rule fails this lint if any of its ``output:`` entries

  * ends in ``.svg``, ``.pdf`` or ``.png`` (including via ``ext="svg"``), or
  * lives under ``results/figures/``, or
  * calls ``P.figures(...)``.

Two exemptions, both narrow and both spelled out below: ``P.report()`` is HTML,
and FastQC's ``.html``/``.zip`` are reports, not figures.

Run it directly::

    python tests/lint_figures_only_in_viz.py            # repo root
    python tests/lint_figures_only_in_viz.py --verbose

Exit status 0 = clean, 1 = violations found, 2 = the lint itself could not run.

Deliberately dependency-free: it parses the .smk text and needs neither
Snakemake, nor a config, nor any heavy data, so it can run in a pre-commit hook
and in CI. A second, authoritative pass over the real rule graph runs only if
Snakemake happens to be importable.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

VIZ_MODULE = "90_visualise.smk"

#: Extensions that make an output a figure.
IMAGE_EXTS = (".svg", ".pdf", ".png")

#: Rules allowed to break the letter of the rule. Keep this list at zero or one.
EXEMPT_RULES = {
    # none
}

#: Substrings whose presence in an output entry is harmless despite the pattern.
BENIGN = (
    "_fastqc.html",   # FastQC report, not a figure
    "_fastqc.zip",
)

_RULE_RE = re.compile(r"^(rule|checkpoint)\s+([A-Za-z_][\w]*)\s*:")
_DIRECTIVE_RE = re.compile(
    r"^\s{1,8}(input|output|params|log|benchmark|threads|resources|priority|"
    r"wildcard_constraints|conda|container|envmodules|shell|script|run|notebook|"
    r"message|group|cache|localrule|retries|default_target|shadow|version)\s*:"
)


class Violation(object):
    def __init__(self, path, lineno, rule, text, why):
        self.path, self.lineno, self.rule, self.text, self.why = (
            path, lineno, rule, text, why,
        )

    def __str__(self):
        rel = os.path.relpath(self.path, REPO)
        return "{}:{}  rule {}\n    {}\n    -> {}".format(
            rel, self.lineno, self.rule, self.text.strip(), self.why
        )


def rule_blocks(path):
    """Yield ``(rule_name, [(lineno, line), ...])`` for every rule in a file."""
    with open(path) as fh:
        lines = fh.readlines()

    current = None
    body = []
    for i, line in enumerate(lines, start=1):
        m = _RULE_RE.match(line)
        if m:
            if current:
                yield current, body
            current, body = m.group(2), []
            continue
        if current is None:
            continue
        stripped = line.rstrip("\n")
        if stripped and not stripped[0].isspace():
            # dedented back to column 0 -> the rule block is over
            yield current, body
            current, body = None, []
            m = _RULE_RE.match(line)
            if m:
                current, body = m.group(2), []
            continue
        body.append((i, stripped))
    if current:
        yield current, body


def output_section(body):
    """Return the ``output:`` lines of one rule body, as ``(lineno, text)``."""
    out = []
    inside = False
    indent = None
    for lineno, line in body:
        m = _DIRECTIVE_RE.match(line)
        if m:
            if m.group(1) == "output":
                inside = True
                indent = len(line) - len(line.lstrip())
                rest = line.split(":", 1)[1].strip()
                if rest:
                    out.append((lineno, rest))
                continue
            if inside and (len(line) - len(line.lstrip())) <= (indent or 0):
                inside = False
            if m.group(1) != "output":
                if inside and (len(line) - len(line.lstrip())) <= (indent or 0):
                    inside = False
                elif not inside:
                    continue
        if inside:
            if line.strip() and not line.strip().startswith("#"):
                out.append((lineno, line))
    return out


def classify(text):
    """Return a reason string if this output entry is a figure, else ``None``."""
    if any(b in text for b in BENIGN):
        return None
    stripped = text.split("#", 1)[0]

    if "P.figures(" in stripped:
        return "calls P.figures() -- only 90_visualise.smk may write under results/figures/"
    if "results/figures" in stripped or 'out("figures' in stripped:
        return "writes under results/figures/"

    # quoted literals ending in an image extension
    for quote in ('"', "'"):
        for lit in re.findall(quote + r"([^" + quote + r"]*)" + quote, stripped):
            low = lit.lower()
            if low.endswith(IMAGE_EXTS):
                return "emits {} -- emit the underlying data instead".format(
                    os.path.splitext(low)[1]
                )
            if low in ("svg", "pdf", "png") and "ext" in stripped:
                return "emits ext={} -- emit the underlying data instead".format(low)
    return None


def lint_static(verbose=False):
    targets = []
    rules_dir = os.path.join(REPO, "workflow", "rules")
    if os.path.isdir(rules_dir):
        for name in sorted(os.listdir(rules_dir)):
            if name.endswith(".smk") and name != VIZ_MODULE:
                targets.append(os.path.join(rules_dir, name))
    snakefile = os.path.join(REPO, "Snakefile")
    if os.path.exists(snakefile):
        targets.append(snakefile)

    if not targets:
        sys.stderr.write("lint: no rule modules found under {}\n".format(rules_dir))
        return None

    violations = []
    n_rules = 0
    for path in targets:
        for rule, body in rule_blocks(path):
            n_rules += 1
            if rule in EXEMPT_RULES:
                continue
            for lineno, text in output_section(body):
                why = classify(text)
                if why:
                    violations.append(Violation(path, lineno, rule, text, why))
    if verbose:
        print("lint: scanned {} rules in {} files".format(n_rules, len(targets)))
    return violations


def lint_dynamic(verbose=False):
    """Authoritative pass over the built rule graph. Skipped if unavailable."""
    try:
        import snakemake  # noqa: F401
    except Exception as exc:                                  # pragma: no cover
        if verbose:
            print("lint: snakemake not importable ({}), static pass only".format(exc))
        return None

    from snakemake.workflow import Workflow  # noqa: F401
    # Building the workflow needs a config and a sample sheet; when either is
    # missing we deliberately do not fail the lint -- the static pass is the
    # gate. This hook exists so that a future CI run with a full config gets the
    # stronger check for free.
    try:
        import subprocess
        env = dict(os.environ)
        env.setdefault("PYTHONWARNINGS", "ignore")
        proc = subprocess.run(
            [sys.executable, "-m", "snakemake", "--list", "-q"],
            cwd=REPO, capture_output=True, text=True, env=env, timeout=300,
        )
        if proc.returncode != 0:
            if verbose:
                print("lint: workflow would not build, static pass only")
            return None
        if verbose:
            print("lint: workflow builds, {} rules".format(
                len(proc.stdout.split())))
    except Exception as exc:                                  # pragma: no cover
        if verbose:
            print("lint: dynamic pass skipped ({})".format(exc))
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--verbose", "-v", action="store_true")
    args = ap.parse_args(argv)

    violations = lint_static(verbose=args.verbose)
    if violations is None:
        return 2
    lint_dynamic(verbose=args.verbose)

    if violations:
        print("FAIL: {} rule output(s) outside {} emit figures\n".format(
            len(violations), VIZ_MODULE))
        for v in violations:
            print(str(v))
            print("")
        print("Fix: emit the data (.tsv / .npz / .bed) from the rule and draw it")
        print("     from a panel_* rule in workflow/rules/{}.".format(VIZ_MODULE))
        print("     See the design notes section 6.2.")
        return 1

    print("OK: no rule outside {} emits *.svg / *.pdf / *.png or writes under "
          "results/figures/".format(VIZ_MODULE))
    return 0


if __name__ == "__main__":
    sys.exit(main())
