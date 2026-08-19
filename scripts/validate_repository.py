from pathlib import Path
import csv
import re


ROOT = Path(__file__).resolve().parents[1]
DATASET = "Gene_Expression_Analysis_and_Disease_Relationship_Synthetic.csv"
MODELS = ["rf_celltype.rds", "rf_disease_basic.rds", "rf_disease_ct.rds"]
REQUIRED = [
    "README.md",
    "app.R",
    "01_train_models.R",
    DATASET,
    "gene_input_example.csv",
    "assets/app-preview.png",
    "www/TME_workflow.png",
    ".Renviron.example",
    *MODELS,
]

missing = [name for name in REQUIRED if not (ROOT / name).is_file()]
assert not missing, f"Missing runtime files: {missing}"

with (ROOT / DATASET).open(encoding="utf-8-sig", newline="") as handle:
    reader = csv.reader(handle)
    header = next(reader)
    row_count = sum(1 for _ in reader)

features = {
    "Gene_A_Oncogene",
    "Gene_B_Immune",
    "Gene_C_Stromal",
    "Gene_D_Therapy",
    "Gene_E_Housekeeping",
    "Pathway_Score_Inflam",
    "Cell_Type",
    "Disease_Status",
}
assert features.issubset(header), f"Dataset schema is missing: {sorted(features - set(header))}"
assert row_count == 3000, f"Expected 3,000 synthetic cells, found {row_count:,}"
assert all((ROOT / model).stat().st_size > 10_000 for model in MODELS), "Model asset is unexpectedly small"

app_text = (ROOT / "app.R").read_text(encoding="utf-8")
for filename in [DATASET, "gene_input_example.csv", *MODELS]:
    assert filename in app_text, f"app.R does not reference {filename}"
assert 'Sys.getenv("GEMINI_API_KEY")' in app_text, "AI key must come from the environment"

secret_pattern = re.compile(r"(?i)(api[_-]?key|token|password)\s*[:=]\s*['\"][^'\"]+['\"]")
assert not secret_pattern.search(app_text), "Possible hard-coded secret in app.R"
for forbidden in (".Renviron", ".Rhistory", ".RData", "rsconnect"):
    assert not (ROOT / forbidden).exists(), f"Remove local-only item: {forbidden}"

print(f"OK: TME Explorer complete; {row_count:,} synthetic cells and {len(MODELS)} models validated")
