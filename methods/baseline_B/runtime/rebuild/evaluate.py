#!/usr/bin/env python3
from __future__ import annotations

import copy
import gzip
import hashlib
import json
import platform
import sys
import time
from pathlib import Path
from typing import Any

import numpy as np
import scipy
from scipy.stats import binomtest


ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / "restored/candidate"
SCRIPTS = CANDIDATE / "scripts"
sys.path.insert(0, str(SCRIPTS))
import evaluate_d5_successor as old

LABELS = ("FULL", "HALF_A", "HALF_B")
COMPARISONS = (
    "FULL_BASIS_ON_FULL",
    "HALF_A_BASIS_ON_HALF_A",
    "HALF_B_BASIS_ON_HALF_B",
    "HALF_A_BASIS_ON_HALF_B_HELD_OUT",
    "HALF_B_BASIS_ON_HALF_A_HELD_OUT",
)


def canon(x: Any) -> bytes:
    return (json.dumps(x, sort_keys=True, separators=(",", ":")) + "\n").encode()


def sem(x: Any) -> str:
    return hashlib.sha256(canon(x)).hexdigest()


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        while b := f.read(8 << 20):
            h.update(b)
    return h.hexdigest()


def load(path: Path):
    if path.suffix == ".gz":
        with gzip.open(path, "rb") as f:
            return json.loads(f.read())
    return json.loads(path.read_text())


def basis(fit):
    return np.asarray(fit["basis"], float) if fit.get("candidate_returned") else None


def bases(result):
    out = {k: basis(result["candidate_fits"][k]) for k in LABELS}
    return None if any(v is None for v in out.values()) else out


def truth_free(case):
    c = copy.deepcopy(case)
    c["truth"] = {"basis": None, "used_by_candidate_generation_or_decision": False}
    return c


def support_mask(q, b, cutoff=0.003):
    n = np.rint(q @ np.linalg.inv(b).T)
    return np.linalg.norm(q - n @ b.T, axis=1) <= cutoff


def paired(q, parent, fine, maximum_p):
    p, f = support_mask(q, parent), support_mask(q, fine)
    plus = int(np.sum(f & ~p))
    minus = int(np.sum(p & ~f))
    discord = plus + minus
    pv = float(binomtest(plus, discord, 0.5, alternative="greater").pvalue) if discord else 1.0
    gain = float(np.mean(f) - np.mean(p))
    return {
        "parent_support": float(np.mean(p)),
        "completion_support": float(np.mean(f)),
        "gain": gain,
        "completion_only": plus,
        "parent_only": minus,
        "p": pv,
        "passes": gain > 0 and pv <= maximum_p,
    }


def equivalence(a, b, tol=0.03, det=1):
    return old.primitive_equivalence(a, b, tol, det)


def consensus(bs, tol):
    pairs = {
        "FULL_TO_HALF_A": equivalence(bs["FULL"], bs["HALF_A"], tol),
        "FULL_TO_HALF_B": equivalence(bs["FULL"], bs["HALF_B"], tol),
        "HALF_A_TO_HALF_B": equivalence(bs["HALF_A"], bs["HALF_B"], tol),
    }
    return {"passes": all(x["primitive_lattice_equivalent"] for x in pairs.values()),
            "pairwise": pairs}


def phases(q, bs, successor):
    own = {k: old.phase_specificity(q[k], bs[k], successor) for k in LABELS}
    held = {
        "HALF_A_BASIS_ON_HALF_B_HELD_OUT": old.phase_specificity(q["HALF_B"], bs["HALF_A"], successor),
        "HALF_B_BASIS_ON_HALF_A_HELD_OUT": old.phase_specificity(q["HALF_A"], bs["HALF_B"], successor),
    }
    return own, held


def evidence(q, parent, fine, successor, rule):
    con = consensus(fine, rule["maximum_basis_integer_deviation"])
    own, held = phases(q, fine, successor)
    phase_pass = all(x["passes"] for x in [*own.values(), *held.values()])
    pairs = {
        "FULL_BASIS_ON_FULL": paired(q["FULL"], parent["FULL"], fine["FULL"], rule["paired_support_maximum_p"]),
        "HALF_A_BASIS_ON_HALF_A": paired(q["HALF_A"], parent["HALF_A"], fine["HALF_A"], rule["paired_support_maximum_p"]),
        "HALF_B_BASIS_ON_HALF_B": paired(q["HALF_B"], parent["HALF_B"], fine["HALF_B"], rule["paired_support_maximum_p"]),
        "HALF_A_BASIS_ON_HALF_B_HELD_OUT": paired(q["HALF_B"], parent["HALF_A"], fine["HALF_A"], rule["paired_support_maximum_p"]),
        "HALF_B_BASIS_ON_HALF_A_HELD_OUT": paired(q["HALF_A"], parent["HALF_B"], fine["HALF_B"], rule["paired_support_maximum_p"]),
    }
    gains = [pairs[k]["gain"] for k in COMPARISONS]
    all_paired = all(pairs[k]["passes"] for k in COMPARISONS)
    complexity = old.primitive_complexity_audit(q["FULL"], fine["FULL"], own["FULL"], successor)
    strong = (con["passes"] and phase_pass and all_paired
              and min(gains) >= rule["strong_completion_minimum_gain"]
              and complexity["passes"])
    weak = (con["passes"] and all_paired
            and min(gains) >= rule["weak_completion_minimum_gain"]
            and not strong)
    return {
        "consensus": con,
        "all_phase_gates_pass": phase_pass,
        "own_phase": own,
        "held_phase": held,
        "paired": pairs,
        "minimum_gain": float(min(gains)),
        "mean_gain": float(np.mean(gains)),
        "complexity": complexity,
        "strong": bool(strong),
        "weak": bool(weak),
    }


def member(kind, scale, bs, branch=None, parent=None, ev=None):
    return {
        "kind": kind,
        "scale": scale,
        "branch": branch,
        "bases": {k: bs[k].tolist() for k in LABELS},
        "parent_bases": None if parent is None else {k: parent[k].tolist() for k in LABELS},
        "evidence": ev,
    }


def full_basis(m):
    return np.asarray(m["bases"]["FULL"], float)


def cluster(members, tol):
    fams = []
    for m in sorted(members, key=lambda x: (-x["scale"], x["kind"], -1 if x["branch"] is None else x["branch"])):
        found = None
        for f in fams:
            if equivalence(np.asarray(f["basis"]), full_basis(m), tol)["primitive_lattice_equivalent"]:
                found = f
                break
        if found is None:
            found = {"index": len(fams), "basis": full_basis(m).tolist(), "members": []}
            fams.append(found)
        found["members"].append(m)
    for f in fams:
        f["scales"] = sorted({x["scale"] for x in f["members"]}, reverse=True)
        f["kinds"] = sorted({x["kind"] for x in f["members"]})
    return fams


def compact_direct(d):
    return {
        "decision": d["decision"],
        "reason": d["decision_reason"],
        "fits": {k: {
            "candidate_returned": d["candidate_fits"][k].get("candidate_returned"),
            "basis": d["candidate_fits"][k].get("basis"),
            "seed": d["candidate_fits"][k].get("seed"),
        } for k in LABELS},
        "consensus": d["primitive_lattice_consensus"],
        "own_phase": d["own_feed_phase_specificity"],
        "held_phase": d["held_out_phase_specificity"],
        "complexity": d["primitive_complexity_audit"],
    }


def completion_branches(scale, direct, case, portability, frozen, successor, rule):
    parent = bases(direct)
    if parent is None or not direct["primitive_lattice_consensus"]["passes"]:
        return []
    if not direct["own_feed_phase_specificity"].get("FULL", {}).get("passes", False):
        return []
    q = {k: old.unpack_feed(case["feeds"][k])[0] for k in LABELS}
    rows = []
    for i, h in enumerate(old.upper_hnf_coarsenings(2)):
        seeds = {k: parent[k] @ np.linalg.inv(h) for k in LABELS}
        seed_ev = evidence(q, parent, seeds, successor, rule)
        refined = {}
        failed = None
        for k in LABELS:
            try:
                refined[k] = portability.FROZEN_D5.refine_lattice(q[k], seeds[k], frozen)[0]
            except Exception as exc:
                failed = f"{type(exc).__name__}: {exc}"
                break
        refined_ev = None if failed else evidence(q, parent, refined, successor, rule)
        strong = bool(refined_ev and refined_ev["strong"])
        refined_weak = bool(refined_ev and refined_ev["weak"])
        seed_weak = bool(seed_ev["weak"] and not strong and not refined_weak)
        chosen_bs = refined if strong or refined_weak else seeds
        chosen_ev = refined_ev if strong or refined_weak else seed_ev
        kind = ("STRONG_COMPLETION" if strong else
                "WEAK_COMPLETION" if refined_weak or seed_weak else None)
        rows.append({
            "branch": i,
            "hnf": h.tolist(),
            "failure": failed,
            "seed_evidence": seed_ev,
            "refined_evidence": refined_ev,
            "provisional_kind": kind,
            "provisional_member": None if kind is None else member(kind, scale, chosen_bs, i, parent, chosen_ev),
        })
    if len(rows) == 7:
        winners = []
        for name in COMPARISONS:
            values = [r["seed_evidence"]["paired"][name]["gain"] for r in rows]
            mx = max(values)
            ws = [i for i, v in enumerate(values) if abs(v - mx) <= 1e-12]
            winners.append(ws)
        unique = winners[0][0] if all(len(w) == 1 for w in winners) and len({w[0] for w in winners}) == 1 else None
        for r in rows:
            if r["provisional_kind"] == "WEAK_COMPLETION" and r["branch"] != unique:
                r["provisional_kind"] = None
                r["provisional_member"] = None
                r["weak_rejection"] = "NONUNIQUE_EQUAL_INDEX_GAIN"
            r["unanimous_seed_gain_winner"] = unique
    return rows


def finite_index(fine, coarse, tol):
    t = np.linalg.inv(fine) @ coarse
    n = np.rint(t).astype(int)
    dev = float(np.max(np.abs(t - n)))
    det = int(round(np.linalg.det(n)))
    return {"deviation": dev, "determinant": det,
            "passes": dev <= tol and abs(det) == 2}


def evaluate_case(case, truth, portability, frozen, successor, rule):
    start = time.monotonic()
    c = truth_free(case)
    direct_members, completion_members, scales = [], [], []
    for floor in rule["seed_floors_cycles_per_angstrom"]:
        fr = copy.deepcopy(frozen)
        fr["lattice_search"]["difference_norm_min_cycles_per_angstrom"] = floor
        t = time.monotonic()
        d = old.evaluate_case(c, portability, fr, successor)
        bs = bases(d)
        if d["decision"] == "LATTICE_RECOVERED" and bs is not None:
            direct_members.append(member("DIRECT", floor, bs))
        branches = completion_branches(floor, d, c, portability, fr, successor, rule)
        completion_members += [r["provisional_member"] for r in branches if r["provisional_member"]]
        scales.append({"floor": floor, "elapsed_seconds": time.monotonic() - t,
                       "direct": compact_direct(d), "completions": branches})

    tol = rule["maximum_basis_integer_deviation"]
    minimum = rule["minimum_persistent_scales"]
    df = cluster(direct_members, tol)
    cf = cluster(completion_members, tol)
    persistent_direct = [f for f in df if len(f["scales"]) >= minimum]
    single_direct = [f for f in df if len(f["scales"]) < minimum]
    strong = [f for f in cf if sum(m["kind"] == "STRONG_COMPLETION" for m in f["members"]) >= minimum]
    weak = [f for f in cf if f not in strong and sum(m["kind"] == "WEAK_COMPLETION" for m in f["members"]) >= minimum]
    one_strong = [f for f in cf if f not in strong and any(m["kind"] == "STRONG_COMPLETION" for m in f["members"])]
    admitted = cluster(
        [m for f in persistent_direct for m in f["members"]] +
        [m for f in strong for m in f["members"] if m["kind"] == "STRONG_COMPLETION"], tol)
    admitted = [f for f in admitted if len(f["scales"]) >= minimum]
    suppressed = set()
    suppressions = []
    for fine in admitted:
        if "STRONG_COMPLETION" not in fine["kinds"]:
            continue
        for coarse in admitted:
            if fine is coarse or "DIRECT" not in coarse["kinds"]:
                continue
            rel = finite_index(np.asarray(fine["basis"]), np.asarray(coarse["basis"]), tol)
            if rel["passes"]:
                suppressed.add(coarse["index"])
                suppressions.append({"fine": fine["index"], "coarse": coarse["index"], "relation": rel})
    survivors = [f for f in admitted if f["index"] not in suppressed]

    def equiv_fam(a, b):
        return equivalence(np.asarray(a["basis"]), np.asarray(b["basis"]), tol)["primitive_lattice_equivalent"]

    weak_conflict = [f for f in weak if not any(equiv_fam(f, s) for s in survivors)]
    strong_conflict = [f for f in one_strong if survivors and all(not equiv_fam(f, s) for s in survivors)]
    direct_conflict = [f for f in single_direct if survivors and all(not equiv_fam(f, s) for s in survivors)]
    if direct_conflict:
        decision, reason = "AMBIGUOUS_LATTICE", "SINGLE_SCALE_DIRECT_CONFLICT_WITNESS"
    elif strong_conflict:
        decision, reason = "AMBIGUOUS_LATTICE", "SINGLE_SCALE_STRONG_COMPLETION_CONFLICT_WITNESS"
    elif weak_conflict:
        decision, reason = "AMBIGUOUS_LATTICE", "PERSISTENT_UNIQUE_WEAK_COMPLETION_EVIDENCE"
    elif len(survivors) > 1:
        decision, reason = "AMBIGUOUS_LATTICE", "MULTIPLE_PERSISTENT_FAMILIES"
    elif len(survivors) == 1:
        decision, reason = "LATTICE_RECOVERED", "ONE_PERSISTENT_FAMILY_SURVIVES"
    else:
        decision, reason = "INSUFFICIENT_SIGNAL", "NO_PERSISTENT_FAMILY"

    decision_record = {"decision": decision, "reason": reason, "truth_consulted": False}
    decision_record["semantic_sha256"] = sem(decision_record)
    scoring = {"truth_available": truth is not None,
               "after_decision_sha256": decision_record["semantic_sha256"]}
    if truth is not None and decision == "LATTICE_RECOVERED":
        f = survivors[0]
        m = sorted(f["members"], key=lambda x: -x["scale"])[0]
        per = {k: equivalence(np.asarray(truth), np.asarray(m["bases"][k]), tol) for k in LABELS}
        scoring["per_feed"] = per
        scoring["correct"] = all(x["primitive_lattice_equivalent"] for x in per.values())
    elif truth is not None:
        scoring["correct"] = False
    required = case["required_outcome"]
    satisfied = bool(scoring.get("correct")) if required == "RECOVER_PRIMITIVE_LATTICE" else decision != "LATTICE_RECOVERED"
    return {
        "case_id": case["case_id"],
        "required_outcome": required,
        "decision": decision_record,
        "scales": scales,
        "families": {
            "direct": df, "completion": cf, "admitted": admitted,
            "survivors": survivors, "weak_conflicts": weak_conflict,
            "strong_single_scale_conflicts": strong_conflict,
            "direct_single_scale_conflicts": direct_conflict,
            "suppressions": suppressions,
        },
        "truth_scoring": scoring,
        "required_outcome_satisfied": satisfied,
        "incorrect_recovery": decision == "LATTICE_RECOVERED" and not satisfied,
        "elapsed_seconds": time.monotonic() - start,
    }


def main():
    work = ROOT / "rebuild/work"
    inherited_path = CANDIDATE / "inputs/NFC_CRYST_D45_SUCCESSOR_FIXED_CONTROL_FEEDS_0_1_0.json.gz"
    six_path = ROOT / "restored/diag/inputs/NFC_CRYST_6MFU_FIXED_D45_CASE_0_1_0.json.gz"
    score_path = ROOT / "restored/diag/prior_result/NFC_CRYST_6MFU_ARCHIVED_D45_D5_POST_COMMITMENT_SCORING_RESULT_0_1_0.json"
    new_path = work / "NFC_CRYST_D5_MULTISCALE_ALIAS_CONTROLS_0_1_0.json.gz"
    rule_path = ROOT / "rebuild/d5_multiscale_bidirectional_rule_0_1_0.json"
    frozen_path = SCRIPTS / "frozen_core/d5_reciprocal_rule_freeze_0_1_0.json"
    successor_path = SCRIPTS / "d5_successor_rule_0_1_0.json"
    inherited, six, new = load(inherited_path), load(six_path), load(new_path)
    score, rule, frozen, successor = load(score_path), load(rule_path), load(frozen_path), load(successor_path)
    portability = old.load_portability(SCRIPTS / "pilatus_portability.py")
    items = []
    for c in inherited["real_cases"]:
        items.append(("INHERITED_REAL", "REAL", c, c["truth"]["basis"]))
    for model, cases in inherited["synthetic_models"].items():
        for c in cases:
            items.append(("INHERITED_SYNTHETIC", model, c, c["truth"]["basis"]))
    six = copy.deepcopy(six)
    six["historical_required_outcome"] = six["required_outcome"]
    six["required_outcome"] = "RECOVER_PRIMITIVE_LATTICE"
    items.append(("SIX_MFU_DEVELOPMENT", "REAL", six,
                  score["conventional_reference"]["truth_basis"]))
    for c in new["cases"]:
        items.append(("NEW_SYNTHETIC", "ALIAS", c, c["truth"]["basis"]))
    bindings = {
        "inherited": sha(inherited_path), "six_mfu": sha(six_path),
        "six_mfu_scoring": sha(score_path), "new_controls": sha(new_path),
        "rule": sha(rule_path), "frozen_rule": sha(frozen_path),
        "successor_rule": sha(successor_path),
    }
    checkpoints = work / "checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)
    results = []
    for group, model, case, truth in items:
        contract = sem({"bindings": bindings, "group": group, "model": model, "case": case})
        path = checkpoints / f"{case['case_id']}_{contract[:20]}.json"
        if path.exists():
            r = load(path)["result"]
            print("CACHE", group, model, case["case_id"], r["decision"]["decision"], flush=True)
        else:
            print("EVALUATE", group, model, case["case_id"], flush=True)
            r = evaluate_case(case, truth, portability, frozen, successor, rule)
            path.write_bytes(canon({"contract": contract, "result": r}))
            print("RESULT", case["case_id"], r["decision"]["decision"],
                  r["required_outcome_satisfied"], flush=True)
        r = copy.deepcopy(r)
        r["group"], r["model"] = group, model
        results.append(r)
    census = {
        "case_count": len(results),
        "required_satisfied": sum(r["required_outcome_satisfied"] for r in results),
        "incorrect_recoveries": sum(r["incorrect_recovery"] for r in results),
        "decisions": {d: sum(r["decision"]["decision"] == d for r in results)
                      for d in ("LATTICE_RECOVERED", "AMBIGUOUS_LATTICE", "INSUFFICIENT_SIGNAL")},
    }
    body = {
        "artifact_id": "NFC_CRYST_D5_MULTISCALE_BIDIRECTIONAL_INDEX_TWO_EVALUATION_0_1_0",
        "scientific_scope": "OPEN_EXPLORATORY_DEVELOPMENT_NOT_CONFIRMATORY",
        "principal_outcome": "PASSES_ALL_CURRENT_CONTROLS" if census["required_satisfied"] == len(results) else "DOES_NOT_PASS_ALL_CURRENT_CONTROLS",
        "bindings": bindings,
        "results": results,
        "census": census,
        "historical_boundary": rule["historical_boundary"],
        "environment": {"python": platform.python_version(), "numpy": np.__version__,
                        "scipy": scipy.__version__, "platform": platform.platform()},
    }
    out = body | {"semantic_sha256": sem(body)}
    path = work / "NFC_CRYST_D5_MULTISCALE_BIDIRECTIONAL_INDEX_TWO_EVALUATION_0_1_0.json"
    path.write_bytes(canon(out))
    print(json.dumps({"output": str(path), "semantic_sha256": out["semantic_sha256"],
                      "principal_outcome": out["principal_outcome"], "census": census},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
