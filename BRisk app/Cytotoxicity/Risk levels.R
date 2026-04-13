ref_df = read.csv("Cytotoxicity_data.csv")

Group_I_df = subset(ref_df, panC_Group == "Group_I")
Group_II_df = subset(ref_df, panC_Group == "Group_II")
Group_III_df = subset(ref_df, panC_Group == "Group_III")
Group_IV_df = subset(ref_df, panC_Group == "Group_IV ")
Group_V_df = subset(ref_df, panC_Group == "Group_V")
Group_VII_df = subset(ref_df, panC_Group == "Group_VII")

ref_ecdf <- ecdf(ref_df$Average_Cell_Viability_F)

Group_I_df$percentile_ref <- ref_ecdf(Group_I_df$Average_Cell_Viability_F) * 100
Group_II_df$percentile_ref <- ref_ecdf(Group_II_df$Average_Cell_Viability_F) * 100
Group_III_df$percentile_ref <- ref_ecdf(Group_III_df$Average_Cell_Viability_F) * 100
Group_IV_df$percentile_ref <- ref_ecdf(Group_IV_df$Average_Cell_Viability_F) * 100
Group_V_df$percentile_ref <- ref_ecdf(Group_V_df$Average_Cell_Viability_F) * 100
Group_VII_df$percentile_ref <- ref_ecdf(Group_VII_df$Average_Cell_Viability_F) * 100

Group_I_pc = median(Group_I_df$percentile_ref)
Group_II_pc = median(Group_II_df$percentile_ref)
Group_III_pc = median(Group_III_df$percentile_ref)
Group_IV_pc = median(Group_IV_df$percentile_ref)
Group_V_pc = median(Group_V_df$percentile_ref)
Group_VII_pc = median(Group_VII_df$percentile_ref)
