import pandas as pd

df = pd.read_csv('./INPUT2_GO-mapping/output_GOALL_terms_PHILLIPP.log', sep="\t")
print("whole df: ")
print(df.shape)
df_no_duplicates = df.drop_duplicates()
print("drop_duplicates: ")
print(df_no_duplicates.shape)
df_only_duplicates = df[df.duplicated()]
print("duplicated: ")
print(df_only_duplicates.shape)
print(df_only_duplicates)

df_no_duplicates.to_csv("./INPUT2_GO-mapping/output_GOALL_terms_PHILLIPP_NoDupes.log", sep="\t", index=False) # Output durch pfad zur neuen datei ersetzen