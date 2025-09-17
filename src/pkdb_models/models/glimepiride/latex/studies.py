from pathlib import Path
import pandas as pd


def create_latex_table(df: pd.DataFrame) -> str:
    """Transform DataFrame into LaTeX table."""

    # Preamble:
    latex_header = r"""
\begin{table}[H]
\centering
\tabcolsep=3.5pt\relax
\scriptsize
\begin{threeparttable} 

\caption{\textbf{Summary of studies for modeling.} 
Overview of study identifiers, PK-DB IDs, administration routes, dosing regimens, doses [mg], co-administered drugs 
(\emph{Co-admin.}), and participant characteristics, including health status, renal impairment (\emph{Ren. imp.}), 
type 2 diabetes mellitus (\emph{T2DM}), and the studied genotypes\alleles (\emph{Alleles}).}
\label{table:curated_data_overview}

\begin{tabularx}{\textwidth}{
  p{2.8cm}  % Study
  p{1.7cm}  % PK-DB ID
  p{1.5cm}  % Route
  p{1.0cm}  % Dosing
  p{1.0cm}  % Dose [mg]
  p{1.8cm}  % Co-Admin
  p{1.1cm}  % Healthy
  p{0.8cm}  % Ren. imp.
  p{1.0cm}  % T2DM
  p{0.8cm}  % Alleles
}
\toprule
""".strip('\n')

    # Selecting needed columns
    df = df.loc[:, [
                       "Study identifier",
                       "PK-DB identifier",
                       "dose route",
                       "dosing",
                       "dose [mg]",
                       "co-administration",
                       "healthy",
                       "renal impairment",
                       "diabetes mellitus 2",
                       "Genotype"
                   ]].copy()

    # Renaming columns
    df = df.rename(columns={
        "Study identifier": "Study",
        "PK-DB identifier": "PK-DB ID",
        "dose route": "Route",
        "dosing": "Dosing",
        "dose [mg]": "Dose [mg]",
        "co-administration": "Co-admin.",
        "healthy": "Healthy",
        "renal impairment": "Ren. imp.",
        "diabetes mellitus 2": "T2DM",
        "Genotype": "Alleles"
    })

    # Limit rows
    df = df.drop(df.index[4]).iloc[1:20].copy()

    # Convert True/False to checkmarks
    checkmark_columns = ["Healthy", "Ren. imp.", "T2DM"]
    for col in checkmark_columns:
        df[col] = df[col].apply(lambda x: r"\checkmark" if str(x).strip().upper() == "TRUE" else "")

    # Column headers row
    column_names = df.columns.tolist()
    column_headers = " & ".join([f"\\textbf{{{col}}}" for col in column_names]) + r" \\"

    # Genotypes columns "*" -> "\textasteriskcentered"
    if "Alleles" in df.columns:
        df["Alleles"] = df["Alleles"].apply(
            lambda x: x.replace("*", r"\textasteriskcentered") if isinstance(x, str) else x
        )

    # LaTeX body
    latex_body = f"{column_headers}\n\\midrule\n"

    # "Badian1996" footnote

    footnote_marker = r"\tnote{1}"

    for _, row in df.iterrows():
        values = list(row.astype(str).values)

        if "Badian1996" in values[0]:
            values[0] = f"{values[0]} \\cite{{{values[0]}}}{footnote_marker}"
        else:
            values[0] = f"{values[0]} \\cite{{{values[0]}}}"

        # PK-DB ID to clickable link
        values[1] = f"\\href{{https://identifiers.org/pkdb:{values[1]}}}{{{values[1]}}}"

        # Converting leftover "TRUE" or "-" strings
        row_str = " & ".join(v.replace("TRUE", r"\checkmark").replace("-", "")
                             for v in values) + r" \\"
        latex_body += row_str + "\n"

    # End of the tabular portion
    latex_footer = r"""
\bottomrule
\end{tabularx}

\begin{tablenotes} 
\item [1] M1 metabolite was administered. 
\end{tablenotes}

\end{threeparttable} 
\end{table}
""".strip('\n')

    # Combine all parts
    full_latex = latex_header + "\n" + latex_body + latex_footer
    return full_latex


if __name__ == "__main__":
    # Input and output file paths
    tsv_path = Path(__file__).parent / 'glimepiride_studies.tsv'
    latex_path = Path(__file__).parent / 'glimepiride_studies.tex'

    # Load TSV file
    df = pd.read_csv(tsv_path, sep="\t")

    # Replace NaN values with an empty string
    df = df.fillna("")

    # Create LaTeX string
    latex_str = create_latex_table(df)

    # Save LaTeX table to file
    with open(latex_path, 'w') as f_tex:
        f_tex.write(latex_str)
