#!/usr/bin/env python3
"""
合并 featureCounts 输出为表达矩阵
"""
import pandas as pd
import glob
import sys

def main():
    count_files = sorted(glob.glob("results/counts/*_counts.txt"))
    if not count_files:
        print("ERROR: No count files found in results/counts/")
        sys.exit(1)

    all_counts = []
    for f in count_files:
        if "all_samples" in f:
            continue
        df = pd.read_csv(f, sep="\\t", comment="#", low_memory=False)
        sample = f.split("/")[-1].replace("_counts.txt", "")
        if "Geneid" in df.columns and len(df.columns) >= 7:
            count_col = df.columns[6]
            all_counts.append(df[["Geneid", count_col]].rename(
                columns={"Geneid": "gene_id", count_col: sample}
            ))

    merged = all_counts[0]
    for df in all_counts[1:]:
        merged = merged.merge(df, on="gene_id", how="outer")

    merged.fillna(0, inplace=True)
    merged.to_csv("results/counts/all_samples_counts.txt", sep="\\t", index=False)
    print(f"Merged counts saved: results/counts/all_samples_counts.txt")
    print(f"Samples: {', '.join(merged.columns[1:])}")

if __name__ == "__main__":
    main()
