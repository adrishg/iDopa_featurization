# -*- coding: utf-8 -*-
#!/usr/bin/env python3
import os, json, csv, argparse, glob
import numpy as np

REQ_KEYS_JSON = ("affinity_pred_value", "affinity_probability_binary")

# Heuristics for .npz: map possible NPZ keys -> our two outputs
# Add/adjust aliases here if your NPZ uses different names
NPZ_KEY_ALIASES = {
    "affinity_pred_value": ["affinity_pred_value", "pred", "y_pred", "affinity", "score", "value"],
    "affinity_probability_binary": ["affinity_probability_binary", "prob", "probability", "p_binary", "p"],
}

def pick_npz_value(npz, candidates):
    """Return first present key's scalar value from candidates, else None."""
    for k in candidates:
        if k in npz:
            arr = npz[k]
            # np.savez stores arrays; squeeze to scalar if needed
            try:
                val = np.array(arr).squeeze()
                # If it's still an array with >1 element, skip
                if val.shape == ():
                    return float(val)  # cast to float if numeric
                # Sometimes stored as 1-element array
                if val.size == 1:
                    return float(val.reshape(()))
            except Exception:
                # Last resort: try to decode if it’s bytes-like
                try:
                    return float(val)
                except Exception:
                    pass
    return None

def find_candidate_files(folder_path):
    # search depth 1 and 2 for json/npz
    patterns = [
        os.path.join(folder_path, "*.json"),
        os.path.join(folder_path, "*", "*.json"),
        os.path.join(folder_path, "*.npz"),
        os.path.join(folder_path, "*", "*.npz"),
    ]
    candidates = []
    for pat in patterns:
        candidates.extend(glob.glob(pat))
    # Prefer JSON over NPZ if both exist
    candidates.sort(key=lambda p: (not p.endswith(".json"), -os.path.getmtime(p)))
    return candidates

def try_from_json(path):
    try:
        with open(path, "r") as f:
            data = json.load(f)
        if all(k in data for k in REQ_KEYS_JSON):
            return {
                "affinity_pred_value": float(data["affinity_pred_value"]),
                "affinity_probability_binary": float(data["affinity_probability_binary"]),
            }, None
        return None, f"JSON missing keys {REQ_KEYS_JSON}"
    except Exception as e:
        return None, f"JSON parse error: {e}"

def try_from_npz(path):
    try:
        with np.load(path, allow_pickle=False) as npz:
            pred = pick_npz_value(npz, NPZ_KEY_ALIASES["affinity_pred_value"])
            prob = pick_npz_value(npz, NPZ_KEY_ALIASES["affinity_probability_binary"])
            if pred is None or prob is None:
                # Helpful debug: list keys
                return None, f"NPZ missing needed values; keys={list(npz.keys())}"
            return {
                "affinity_pred_value": float(pred),
                "affinity_probability_binary": float(prob),
            }, None
    except Exception as e:
        return None, f"NPZ read error: {e}"

def extract_affinity_values(input_dir, output_csv):
    rows, skipped, details = [], [], []
    folder_names = sorted(d for d in os.listdir(input_dir) if os.path.isdir(os.path.join(input_dir, d)))

    for folder_name in folder_names:
        folder_path = os.path.join(input_dir, folder_name)
        files = find_candidate_files(folder_path)
        if not files:
            skipped.append(folder_name)
            details.append(f"[MISS] {folder_name}: no JSON/NPZ found")
            continue

        used = None
        rec = None
        err = None
        for path in files:
            if path.endswith(".json"):
                rec, err = try_from_json(path)
                if rec: 
                    used = os.path.basename(path)
                    break
            elif path.endswith(".npz"):
                rec, err = try_from_npz(path)
                if rec:
                    used = os.path.basename(path)
                    break

        if rec:
            rows.append({"Tag": folder_name, **rec})
            details.append(f"[OK]   {folder_name}: {used}")
        else:
            skipped.append(folder_name)
            details.append(f"[SKIP] {folder_name}: {err}")

    # write CSV
    with open(output_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["Tag", "affinity_pred_value", "affinity_probability_binary"])
        w.writeheader()
        w.writerows(rows)

    # summary
    print(f"[INFO] CSV written to {output_csv}")
    print(f"[INFO] Folders scanned: {len(folder_names)}")
    print(f"[INFO] Records found:   {len(rows)}")
    print(f"[INFO] Skipped folders: {len(skipped)}")
    if skipped:
        print("[INFO] Skipped list:")
        for name in skipped:
            print(f"  - {name}")
    # Uncomment for detailed per-folder log:
    # for line in details: print(line)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Extract affinity values from JSON/NPZ files in subfolders.")
    ap.add_argument('--input-dir', required=True, help="Parent folder containing result subfolders.")
    ap.add_argument('--output-csv', required=True, help="Output CSV file path.")
    args = ap.parse_args()
    extract_affinity_values(args.input_dir, args.output_csv)

