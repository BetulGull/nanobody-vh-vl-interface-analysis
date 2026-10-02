"""Stage 1.5 — targeted mmCIF/entity verification of Stage-1 metadata anomalies.

Read-only diagnostic. No SASA/geometry/interaction/electrostatics/aggregation/
statistics performed. No cohort membership changed.
"""
import re
from pathlib import Path

import gemmi
import pandas as pd

PROJECT_ROOT = Path("/Users/betulgul/Desktop/staj")

CONFLICT_IDS = [
    "sabdab2_H00QEL0000", "sabdab2_H021EL0000", "sabdab2_H02DAL0000",
    "sabdab2_H02Y5L0000", "sabdab2_H06BKL0000",
]

MHC_H02Y5_ID = "sabdab2_H02Y5L0000"

# These files carry NO _entity.pdbx_description field at all (checked directly in the raw
# mmCIF: only _entity.id / _entity.type are present), so MHC/B2M identity must be read from
# the actual entity sequence, not a description field. B2M_SIGNATURE is a long, distinctive
# substring of human beta-2-microglobulin; HLA_SIGNATURE is the conserved HLA class I
# mature N-terminus.
B2M_SIGNATURE = "IQRTPKIQVYSRHPAENGKSNFLNCYVSGFHPSDIEVDLLKNGERIEKVEHSDLSFSKDWSFYLLYYTEFTPTEKDEYACRVNHVTLSQPKIVKWDRDM"
HLA_SIGNATURE = "SHSMRYF"
VHH_START_RE = re.compile(r"^(GS)?(QVQL|EVQL|QVKL|EVKL|QVEL|AVQL|EVMV|DVQL|QLQL|QVQV|EVQV)")


def read_structure(path: Path):
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    return st


def get_chain_polymer_info(st, chain_id: str):
    """Return (coord_seq_with_numbering, entity_full_seq, entity_description)."""
    model = st[0]
    chain = None
    for c in model:
        if c.name == chain_id:
            chain = c
            break
    if chain is None:
        return None, None, None

    poly = chain.get_polymer()
    coord_residues = [(r.seqid.num, r.seqid.icode.strip(), gemmi.find_tabulated_residue(r.name).one_letter_code.upper()
                        if gemmi.find_tabulated_residue(r.name) else "X")
                       for r in poly]
    coord_seq = "".join(x[2] for x in coord_residues)

    entity_full_seq = None
    entity_description = None
    subchain = chain[0].subchain if len(chain) else None
    # find entity via subchain -> entity mapping
    ent = None
    if len(poly) > 0:
        sub = poly[0].subchain
        for e in st.entities:
            if sub in e.subchains:
                ent = e
                break
    if ent is not None:
        entity_full_seq = gemmi.one_letter_code(ent.full_sequence).upper()
        entity_description = getattr(ent, "name", None) or None

    return coord_residues, entity_full_seq, entity_description


def get_all_protein_entities(st):
    """Return list of dicts: entity name/id, subchains, chains, full_seq, length."""
    model = st[0]
    subchain_to_chain = {}
    for chain in model:
        for res in chain:
            subchain_to_chain.setdefault(res.subchain, set()).add(chain.name)

    out = []
    for ent in st.entities:
        if ent.entity_type != gemmi.EntityType.Polymer:
            continue
        if ent.polymer_type not in (gemmi.PolymerType.PeptideL, gemmi.PolymerType.PeptideD):
            continue
        chains = set()
        for sub in ent.subchains:
            chains |= subchain_to_chain.get(sub, set())
        seq = gemmi.one_letter_code(ent.full_sequence).upper()
        out.append({
            "entity_name": ent.name,
            "chains": sorted(chains),
            "subchains": list(ent.subchains),
            "full_seq": seq,
            "length": len(seq),
        })
    return out


def extract_hchain_details(row):
    path = PROJECT_ROOT / row["extracted_structure_path"]
    if not path.exists():
        return {"error": f"file not found: {path}"}
    st = read_structure(path)
    coord_residues, entity_full_seq, entity_description = get_chain_polymer_info(st, row["Hchain"])
    if coord_residues is None:
        return {"error": f"Hchain {row['Hchain']} not found in structure"}

    coord_seq = "".join(x[2] for x in coord_residues)
    numbering = [(n, ic) for n, ic, _ in coord_residues]
    return {
        "coord_seq": coord_seq,
        "coord_len": len(coord_seq),
        "entity_full_seq": entity_full_seq,
        "entity_full_len": len(entity_full_seq) if entity_full_seq else None,
        "entity_description": entity_description,
        "numbering_min": min(n for n, _ in numbering) if numbering else None,
        "numbering_max": max(n for n, _ in numbering) if numbering else None,
        "has_icode": any(ic for _, ic in numbering),
        "numbering_list": numbering,
    }


def compare_sequences_by_position(instances_details):
    """Compare coord sequences across instances using shared IMGT residue numbers only."""
    # build num->aa maps
    maps = []
    for d in instances_details:
        m = {}
        for (num, ic), aa in zip(d["numbering_list"], d["coord_seq"]):
            m[(num, ic)] = aa
        maps.append(m)

    all_positions = set()
    for m in maps:
        all_positions |= set(m.keys())

    shared_positions = set(maps[0].keys())
    for m in maps[1:]:
        shared_positions &= set(m.keys())

    n_compared = len(shared_positions)
    n_mismatch = 0
    mismatches = []
    for pos in sorted(shared_positions):
        vals = set(m[pos] for m in maps)
        if len(vals) > 1:
            n_mismatch += 1
            mismatches.append((pos, [m[pos] for m in maps]))

    return {
        "n_shared_positions": n_compared,
        "n_mismatched_positions": n_mismatch,
        "mismatches": mismatches[:20],
        "identical_at_shared_positions": n_mismatch == 0,
    }


def main():
    df = pd.read_parquet(PROJECT_ROOT / "data" / "metadata" / "master_stage1.parquet")
    sdh = df[df["source_cohort"] == "SD-H"]

    # =========================================================================
    # 1 & 2: five conflicting SABDAB_IDs — full per-instance extraction + sequence comparison
    # =========================================================================
    conflict_detail_rows = []
    seq_identity_rows = []

    for sid in CONFLICT_IDS:
        grp = sdh[sdh["SABDAB_ID"] == sid].drop_duplicates("PDB")
        instance_details = []
        for _, row in grp.iterrows():
            det = extract_hchain_details(row)
            instance_details.append(det)
            conflict_detail_rows.append({
                "SABDAB_ID": sid,
                "PDB": row["PDB"],
                "INSTANCE": row["INSTANCE"],
                "Hchain": row["Hchain"],
                "Lchain": row["Lchain"],
                "heavy_species": row["heavy_species"],
                "heavy_taxid": row["heavy_taxid"],
                "organism": row["organism"],
                "organism_taxid": row["organism_taxid"],
                "compound": row["compound"],
                "short_header": row["short_header"],
                "coord_seq": det.get("coord_seq"),
                "coord_len": det.get("coord_len"),
                "entity_full_seq": det.get("entity_full_seq"),
                "entity_full_len": det.get("entity_full_len"),
                "entity_description": det.get("entity_description"),
                "numbering_min": det.get("numbering_min"),
                "numbering_max": det.get("numbering_max"),
                "has_icode": det.get("has_icode"),
                "error": det.get("error"),
            })

        valid_details = [d for d in instance_details if "error" not in d]
        n_instances = len(grp)

        if len(valid_details) >= 2:
            cmp_result = compare_sequences_by_position(valid_details)
            unique_full_seqs = set(d["entity_full_seq"] for d in valid_details if d["entity_full_seq"])
            unique_coord_seqs_exact = set(d["coord_seq"] for d in valid_details)
        else:
            cmp_result = {"n_shared_positions": None, "n_mismatched_positions": None, "identical_at_shared_positions": None}
            unique_full_seqs = set()
            unique_coord_seqs_exact = set()

        species_labels = sorted(set(grp["heavy_species"]))
        seq_identity_rows.append({
            "SABDAB_ID": sid,
            "n_instances": n_instances,
            "n_structures_extracted_ok": len(valid_details),
            "n_unique_entity_full_sequences": len(unique_full_seqs),
            "n_unique_coord_sequences_exact_string": len(unique_coord_seqs_exact),
            "n_shared_imgt_positions_compared": cmp_result["n_shared_positions"],
            "n_mismatched_positions_at_shared_sites": cmp_result["n_mismatched_positions"],
            "identical_at_all_shared_positions": cmp_result["identical_at_all_shared_positions"] if "identical_at_all_shared_positions" in cmp_result else cmp_result.get("identical_at_shared_positions"),
            "species_labels_present": ";".join(species_labels),
            "example_mismatches": str(cmp_result.get("mismatches", []))[:500],
        })

    conflict_detail_df = pd.DataFrame(conflict_detail_rows)
    seq_identity_df = pd.DataFrame(seq_identity_rows)
    seq_identity_df.to_csv(PROJECT_ROOT / "results" / "audits" / "stage1_5_sequence_identity_check.csv", index=False)

    # =========================================================================
    # 3: provenance conflict resolution table
    # =========================================================================
    provenance_rows = []
    for _, r in seq_identity_df.iterrows():
        sid = r["SABDAB_ID"]
        grp_all = sdh[sdh["SABDAB_ID"] == sid]  # NOT deduped by PDB — need every instance to catch within-PDB/within-row label disagreement
        species_text_set = set(grp_all["heavy_species"].str.lower())
        taxid_set = set(grp_all["heavy_taxid"])
        text_all = " ".join(grp_all["compound"].fillna("") + " " + grp_all["short_header"].fillna("")).lower()
        humanized_text = "humanized" in text_all or "humanised" in text_all

        # internal self-consistency: does heavy_species TEXT disagree with heavy_taxid NUMBER
        # within the very same row (e.g. species says "vicugna pacos" but taxid says human 9606)?
        camelid_taxids_local = {"9838", "9837", "419612", "9844", "30538", "9835", "9839", "1579311"}
        self_contradictory_rows = grp_all[
            (grp_all["heavy_species"].str.lower().isin(["vicugna pacos", "lama glama", "camelus dromedarius", "camelus bactrianus", "lama", "camelidae"]))
            & (~grp_all["heavy_taxid"].isin(camelid_taxids_local))
        ]

        identical = r["identical_at_all_shared_positions"]
        n_shared = r["n_shared_imgt_positions_compared"]

        if len(self_contradictory_rows) > 0:
            status = "POSSIBLE_DATABASE_INCONSISTENCY"
            bad_pdbs = sorted(set(self_contradictory_rows["PDB"]))
            note = (f"heavy_species TEXT and heavy_taxid NUMBER disagree within the same row for "
                    f"PDB(s) {bad_pdbs}: heavy_species names a camelid species but heavy_taxid points "
                    f"to a non-camelid taxid. This is a self-contradictory metadata record in the source "
                    f"summary CSV, not merely a cross-instance labeling difference.")
        elif identical is True and n_shared and n_shared > 50:
            if len(taxid_set) > 1 and "9606" in taxid_set:
                status = "MIXED_PROVENANCE_SAME_SEQUENCE"
                note = ("Identical variable-domain sequence at every shared IMGT position across "
                        "instances, but heavy_taxid annotation differs (includes both camelid and "
                        "Homo sapiens=9606) across PDB entries for the same SABDAB_ID. This is a "
                        "provenance/annotation inconsistency, not a genuine sequence difference.")
            else:
                status = "NATURAL_CAMELID_PROVENANCE_SUPPORTED"
                note = "Identical sequence at shared positions; species labels vary only among camelid genera/species (e.g. Lama vs Vicugna), consistent with a shared SAbDab2 ID under differing per-entry species annotation."
        elif identical is False:
            status = "POSSIBLE_DATABASE_INCONSISTENCY"
            note = f"{r['n_mismatched_positions_at_shared_sites']} mismatched residue(s) found at shared IMGT positions across instances of the same SABDAB_ID — needs manual sequence-level review; see example_mismatches column."
        else:
            status = "UNRESOLVED_PROVENANCE"
            note = "Insufficient overlapping structural data (too few shared resolved positions, or extraction error) to determine sequence identity directly."

        provenance_rows.append({
            "SABDAB_ID": sid,
            "n_instances": r["n_instances"],
            "species_text_labels_present": ";".join(sorted(species_text_set)),
            "heavy_taxid_values_present": ";".join(sorted(taxid_set)),
            "identical_at_shared_positions": identical,
            "n_shared_positions_compared": n_shared,
            "n_mismatched_positions": r["n_mismatched_positions_at_shared_sites"],
            "humanized_text_present": humanized_text,
            "self_contradictory_species_vs_taxid_rows": len(self_contradictory_rows),
            "recommended_status": status,
            "note": note,
            "sequence_provenance_ambiguity": status in ("MIXED_PROVENANCE_SAME_SEQUENCE", "UNRESOLVED_PROVENANCE", "POSSIBLE_DATABASE_INCONSISTENCY"),
        })

    provenance_df = pd.DataFrame(provenance_rows)
    provenance_df.to_csv(PROJECT_ROOT / "results" / "audits" / "stage1_5_sabdab_id_provenance_conflicts.csv", index=False)

    # =========================================================================
    # 4: H02Y5L0000 chain/entity verification (the "58 MHC rows" claim)
    # =========================================================================
    h02y5_rows = sdh[sdh["SABDAB_ID"] == MHC_H02Y5_ID].drop_duplicates("PDB")
    h02y5_verification = []
    mhc_keywords = ["hla", "mhc", "beta-2-microglobulin", "beta2-microglobulin", "b2-microglobulin", "microglobulin"]

    for _, row in h02y5_rows.iterrows():
        path = PROJECT_ROOT / row["extracted_structure_path"]
        if not path.exists():
            h02y5_verification.append({"PDB": row["PDB"], "error": "file not found"})
            continue
        st = read_structure(path)
        all_entities = get_all_protein_entities(st)
        hchain_det = extract_hchain_details(row)

        hchain_seq = hchain_det.get("entity_full_seq") or hchain_det.get("coord_seq") or ""
        hchain_has_b2m = B2M_SIGNATURE[:40] in hchain_seq
        hchain_has_hla = HLA_SIGNATURE in hchain_seq
        hchain_is_mhc = hchain_has_b2m or hchain_has_hla
        hchain_looks_like_vhh = bool(VHH_START_RE.match(hchain_seq)) and hchain_det.get("coord_len", 0) < 150

        other_entities_desc = []
        mhc_entity_found = None
        for e in all_entities:
            if row["Hchain"] in e["chains"]:
                continue
            seq = e["full_seq"]
            tag = ""
            if B2M_SIGNATURE[:40] in seq:
                tag = "[beta-2-microglobulin]"
            elif HLA_SIGNATURE in seq:
                tag = "[HLA class I heavy chain]"
            elif len(seq) < 15:
                tag = "[short peptide antigen]"
            other_entities_desc.append(f"{e['chains']}:len={e['length']}{tag}")
            if tag in ("[beta-2-microglobulin]", "[HLA class I heavy chain]"):
                mhc_entity_found = f"{e['chains']}:len={e['length']}{tag}"

        h02y5_verification.append({
            "PDB": row["PDB"],
            "INSTANCE": row["INSTANCE"],
            "Hchain": row["Hchain"],
            "Hchain_coord_len": hchain_det.get("coord_len"),
            "Hchain_entity_full_len": hchain_det.get("entity_full_len"),
            "Hchain_seq_preview": (hchain_seq or "")[:60],
            "Hchain_contains_B2M_signature": hchain_has_b2m,
            "Hchain_contains_HLA_signature": hchain_has_hla,
            "Hchain_is_MHC_construct": hchain_is_mhc,
            "Hchain_looks_like_VHH_by_sequence": hchain_looks_like_vhh,
            "Lchain_value": row["Lchain"],
            "Lchain_absent_as_expected_for_SDH": row["Lchain"] in ("NA", None),
            "other_protein_entities": " | ".join(other_entities_desc),
            "mhc_entity_identified_separately": mhc_entity_found,
            "compound_text": row["compound"],
        })

    h02y5_df = pd.DataFrame(h02y5_verification)
    h02y5_df.to_csv(PROJECT_ROOT / "results" / "audits" / "stage1_5_h02y5_chain_verification.csv", index=False)

    return conflict_detail_df, seq_identity_df, provenance_df, h02y5_df


if __name__ == "__main__":
    conflict_detail_df, seq_identity_df, provenance_df, h02y5_df = main()
    pd.set_option("display.width", 220)
    pd.set_option("display.max_colwidth", 60)
    print("=== seq_identity_df ===")
    print(seq_identity_df.to_string(index=False))
    print("\n=== provenance_df ===")
    print(provenance_df[["SABDAB_ID", "n_instances", "species_text_labels_present", "heavy_taxid_values_present", "identical_at_shared_positions", "self_contradictory_species_vs_taxid_rows", "recommended_status"]].to_string(index=False))
    print("\n=== h02y5_df (Hchain identity) ===")
    print(h02y5_df[["PDB", "Hchain", "Hchain_coord_len", "Hchain_is_MHC_construct", "Hchain_looks_like_VHH_by_sequence", "mhc_entity_identified_separately"]].to_string(index=False))
