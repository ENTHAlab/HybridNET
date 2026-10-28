from pathlib import Path
from datetime import datetime
import logging
import math
import numpy as np
import pandas as pd

# The locations of the input data and output data.
PARENT_FILE = r"..\parent.xlsx"
PRODUCT_FILE = r"..\product.xlsx"
OUTPUT_ROOT = r"..\result"
# Parameter settings
THRESHOLD = 0.7
MS2_TOLERANCE_DA = 0.01
LOW_INTENSITY_RATIO = 0.01
CHUNK_SIZE = 1_000_000
# Configure the log format (time + level + message)
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
LOG = logging.getLogger("all_five_similarity")
logging.getLogger("spec2vec").setLevel(logging.ERROR)
# Import the library
try:
    from matchms import Spectrum
    from matchms.similarity import CosineGreedy, ModifiedCosineGreedy, NeutralLossesCosine, FlashSimilarity
    from spec2vec import Spec2Vec
    from gensim.models import Word2Vec
except ImportError as exc:
    raise SystemExit("Missing dependencies: please instal pandas numpy openpyxl matchms spec2vec gensim") from exc
METHODS = ("cosine", "modifiedcosine", "neutralloss", "flashentropy", "spec2vec")
SPEC2VEC_MODEL = r"..\spec2vec_AllPositive_ratio05_filtered_201101_iter_15.model"
# Parse the data
def parse_spectrum(value):
    mzs, intensities = [], []
    if isinstance(value, str):
        for token in value.split():
            try:
                mz, intensity = (float(x) for x in token.split(":", 1))
                if np.isfinite(mz) and np.isfinite(intensity):
                    mzs.append(mz); intensities.append(intensity)
            except (ValueError, TypeError):
                continue
    return mzs, intensities
# export data
def export_spectra(path, spectra):
    rows = [{"ID":s["id"], "mz":s["precursor_mz"],
             "MS2":" ".join(f"{m}:{i}" for m,i in zip(s["mz_list"],s["intensity_list"])),
             "num_peaks":len(s["mz_list"])} for s in spectra]
    pd.DataFrame(rows, columns=["ID","mz","MS2","num_peaks"]).to_csv(path, index=False, encoding="utf-8-sig")
# save data
def save_csv_in_chunks(frame, round_dir, method):
    if frame.empty:
        frame.to_csv(round_dir / f"{method}.csv", index=False, encoding="utf-8-sig")
        return
    total = len(frame)
    if total <= CHUNK_SIZE:
        frame.to_csv(round_dir / f"{method}.csv", index=False, encoding="utf-8-sig")
        return
    num_chunks = math.ceil(total / CHUNK_SIZE)
    for idx in range(num_chunks):
        chunk = frame.iloc[idx * CHUNK_SIZE:(idx + 1) * CHUNK_SIZE]
        chunk.to_csv(round_dir / f"{method}-{idx + 1}.csv", index=False, encoding="utf-8-sig")
    LOG.info("%s: %d rows saved into %d files.", method, total, num_chunks)
# load data
def load_and_preprocess(file_path, kind, prep_dir):
    df = pd.read_excel(file_path)
    missing = {"ID", "mz", "MS2"} - set(df.columns)
    if missing:
        raise ValueError(f"{file_path} is missing columns: {sorted(missing)}")
    parsed = []
    for _, row in df.iterrows():
        mzs, ints = parse_spectrum(row["MS2"])
        if mzs:
            try: precursor = float(row["mz"])
            except (ValueError, TypeError): continue
            parsed.append({"id":row["ID"], "precursor_mz":precursor, "mz_list":mzs,
                           "intensity_list":ints, "type":kind})
    export_spectra(prep_dir / f"{kind}_1_parsed.csv", parsed)
    filtered = []
    for s in parsed:
        it = np.asarray(s["intensity_list"], dtype=float)
        if not it.size or np.max(it) <= 0: continue
        keep = it >= np.max(it) * LOW_INTENSITY_RATIO
        if keep.any():
            filtered.append({**s, "mz_list":np.asarray(s["mz_list"])[keep].tolist(),
                             "intensity_list":it[keep].tolist()})
    export_spectra(prep_dir / f"{kind}_2_filtered_below_1pct.csv", filtered)
    normalized = [{**s, "intensity_list":(np.asarray(s["intensity_list"], dtype=float) /
                   np.max(s["intensity_list"])).tolist()} for s in filtered]
    export_spectra(prep_dir / f"{kind}_3_normalized.csv", normalized)
    return normalized

def to_spectrum(spec):
    peaks = np.column_stack((spec["mz_list"], spec["intensity_list"])).astype(np.float64)
    valid = (peaks[:, 1] > 0) & np.isfinite(peaks[:, 0]) & np.isfinite(peaks[:, 1])
    peaks = peaks[valid]
    if not len(peaks): return None
    peaks = peaks[np.argsort(peaks[:, 0])]
    pmz = spec["precursor_mz"]
    return Spectrum(mz=peaks[:,0], intensities=peaks[:,1],
                    metadata={"precursor_mz":float(pmz) if pmz else 0.0})

def score_value(value):
    try: return float(value["score"])
    except (TypeError, ValueError, KeyError, IndexError):
        try: return float(value[0])
        except (TypeError, ValueError, IndexError): return 0.0
# Similarity calculation.
def build_scorers():
    model = Word2Vec.load(SPEC2VEC_MODEL)
    spec2vec = Spec2Vec(model=model, intensity_weighting_power=0.5,
                        allowed_missing_percentage=10, progress_bar=False)
    return {
        "cosine":CosineGreedy(tolerance=MS2_TOLERANCE_DA, mz_power=0.0, intensity_power=1),
        "modifiedcosine":ModifiedCosineGreedy(tolerance=MS2_TOLERANCE_DA, mz_power=0.0, intensity_power=1),
        "neutralloss":NeutralLossesCosine(tolerance=MS2_TOLERANCE_DA, mz_power=0.0, intensity_power=1,
                                           ignore_peaks_above_precursor=True),
        "flashentropy":FlashSimilarity(score_type="spectral_entropy", matching_mode="hybrid", tolerance=0.01,
            use_ppm=False, remove_precursor=False, noise_cutoff=0.01, normalize_to_half=True,
            merge_within=0.0, dtype=np.float64),
        "spec2vec":spec2vec,
    }

def pair_score(method, scorer, parent, product, parent_sp, product_sp):
    if method == "spec2vec":
        value = scorer.pair(parent_sp, product_sp)
        if not np.isfinite(value): return 0.0
        value = float(value)
        return 0.0 if -1e-6 < value < 0 else value
    result = scorer.pair(parent_sp, product_sp)
    return score_value(result)

def flash_matrix(scorer, p_specs, q_specs):
    n_p, n_q = len(p_specs), len(q_specs)
    if n_p == 0 or n_q == 0:
        return np.zeros((n_p, n_q), dtype=float)
    p_list = [sp for _, sp in p_specs]
    q_list = [sp for _, sp in q_specs]
    try:
        mat = scorer.matrix(p_list, q_list, progress_bar=False)
        mat = np.asarray(mat, dtype=float)
    except TypeError:
        mat = scorer.matrix(p_list, q_list)
        mat = np.asarray(mat, dtype=float)
    if mat.shape != (n_p, n_q):
        if mat.T.shape == (n_p, n_q):
            mat = mat.T
        else:
            mat = mat.reshape(n_p, n_q)
    mat = np.where(np.isfinite(mat), mat, 0.0)
    return mat

def run():
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")  # Generate a timestamp.
    root = Path(OUTPUT_ROOT) / f"analysis_{stamp}"  # Create the output directory structure.
    prep_dir, result_dir = root / "Data Preprocessing", root / "Similarity Results"
    prep_dir.mkdir(parents=True, exist_ok=True)
    result_dir.mkdir(parents=True, exist_ok=True)
    parents = load_and_preprocess(PARENT_FILE, "parent", prep_dir)
    products = load_and_preprocess(PRODUCT_FILE, "product", prep_dir)
    if not parents or not products:
        raise ValueError("No usable data in parent or product spectra; calculation stopped")
    scorers = build_scorers()
    current_parents, remaining = parents, products[:]
    round_no = 1
    pair_columns = ["parent_id", "product_id", "similarity", "parent_mz", "product_mz"]
    while current_parents and remaining:
        round_dir = result_dir / f"{round_no} round"
        round_dir.mkdir(parents=True, exist_ok=True)
        p_specs = [(s, to_spectrum(s)) for s in current_parents]
        q_specs = [(s, to_spectrum(s)) for s in remaining]
        p_specs = [(s, sp) for s,sp in p_specs if sp is not None]
        q_specs = [(s, sp) for s,sp in q_specs if sp is not None]

        flash_scores = flash_matrix(scorers["flashentropy"], p_specs, q_specs)
        LOG.info("Round %d: flashentropy matrix computed, shape=%s", round_no, flash_scores.shape)

        method_rows = {name:[] for name in METHODS}
        for pi, (parent, psp) in enumerate(p_specs, 1):
            if pi % 25 == 0 or pi == len(p_specs):
                LOG.info("%d round：parent data %d/%d", round_no, pi, len(p_specs))
            for qi, (product, qsp) in enumerate(q_specs):
                for method in METHODS:
                    if method == "flashentropy":
                        score = float(flash_scores[pi - 1, qi])
                    else:
                        score = pair_score(method, scorers[method], parent, product, psp, qsp)
                    if np.isfinite(score) and score >= THRESHOLD:
                        method_rows[method].append({"parent_id":parent["id"], "product_id":product["id"],
                            "similarity":round(float(score),6), "parent_mz":parent["precursor_mz"],
                            "product_mz":product["precursor_mz"]})
        combined = []
        for method in METHODS:
            frame = pd.DataFrame(method_rows[method], columns=pair_columns)
            if not frame.empty:
                combined.append(frame)
        if not combined:
            pd.DataFrame(columns=pair_columns).to_csv(
                round_dir / "result_test.csv", index=False, encoding="utf-8-sig")
            LOG.info("No matches above the threshold in round %d; iteration finished.", round_no)
            break
        hits = pd.concat(combined, ignore_index=True)
        best = (hits.sort_values("similarity", ascending=False, kind="mergesort")
                    .drop_duplicates(["parent_id", "product_id"], keep="first"))  # Keep the highest similarity result.
        save_csv_in_chunks(best, round_dir, "result_test")
        matched_ids = set(best["product_id"])
        current_parents = [s for s in remaining if s["id"] in matched_ids]
        remaining = [s for s in remaining if s["id"] not in matched_ids]
        LOG.info("Round %d: %d matched edges, %d product spectra hit, %d remaining.",
                 round_no, len(best), len(matched_ids), len(remaining))
        round_no += 1
    LOG.info("Done. Preprocessing files and per-round result directory：%s", root)

if __name__ == "__main__":
    run()
