import shutil
from pathlib import Path

from pkdb_models.models.data import collect_tsv_files

def collect_glimepiride_data():
    common_parent: Path = Path(__file__).parents[5]
    source_dir = common_parent / "pkdb_data" / "studies" / "glimepiride"
    target_dir = Path(__file__).parent / "glimepiride"

    # collect glimepiride
    collect_tsv_files(source_dir=source_dir, target_dir=target_dir)

    # collect excel files
    for filename in ["cyp2c9_data.xlsx", "dose_dependency_data.xlsx"]:
        shutil.copy2(source_dir / filename, target_dir / filename)

    # collect dapagliflozin
    def is_Kasichayanula2011c(study_name) -> bool:
        return study_name == "Kasichayanula2011c"

    collect_tsv_files(
        source_dir=common_parent / "pkdb_data" / "studies" / "dapagliflozin",
        target_dir=Path(__file__).parent / "dapagliflozin",
        filter_study=is_Kasichayanula2011c,
    )


if __name__ == "__main__":
    collect_glimepiride_data()

