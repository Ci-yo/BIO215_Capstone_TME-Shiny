## 01_train_models.R
## 作用：读入合成基因表达数据 -> 训练 3 个 RF 模型 -> 生成 .rds 和 example CSV

## 0. 环境 & 包 ---------------------------------------------------------
rm(list = ls())

if (!requireNamespace("randomForest", quietly = TRUE)) {
  install.packages("randomForest")
}
library(randomForest)

## 1. 读数据 ------------------------------------------------------------
# 确保这个文件和脚本在同一个目录下
data_df <- read.csv(
  "Gene_Expression_Analysis_and_Disease_Relationship_Synthetic.csv",
  stringsAsFactors = FALSE
)

# 看一下前几行（可选）
# head(data_df)

## 2. 基本预处理 --------------------------------------------------------
# 设置因子水平（固定顺序，后面模型/画图都统一）
data_df$Cell_Type <- factor(
  data_df$Cell_Type,
  levels = c("Cancer", "T_Cell", "Fibroblast")
)

data_df$Disease_Status <- factor(
  data_df$Disease_Status,
  levels = c("Healthy_Control", "Tumor")
)

# 特征工程：三个 score（跟你之前分析保持一致）
data_df$Tumor_score   <- data_df$Gene_A_Oncogene + data_df$Gene_D_Therapy
data_df$Immune_score  <- data_df$Gene_B_Immune   + data_df$Pathway_Score_Inflam
data_df$Stromal_score <- data_df$Gene_C_Stromal

## 3. 划分 train / test -------------------------------------------------
set.seed(123)
n        <- nrow(data_df)
n_train  <- round(0.8 * n)
train_id <- sample.int(n, size = n_train)

train_df <- data_df[train_id, ]
test_df  <- data_df[-train_id, ]

## 4. 随机森林：预测 Cell_Type -----------------------------------------
set.seed(123)
rf_celltype <- randomForest(
  Cell_Type ~ Gene_E_Housekeeping +
    Gene_A_Oncogene +
    Gene_B_Immune +
    Gene_C_Stromal +
    Gene_D_Therapy +
    Pathway_Score_Inflam,
  data       = train_df,
  importance = TRUE,
  ntree      = 500
)

print(rf_celltype)  # 看一下 OOB 精度（可选）

## 5. 随机森林：预测 Disease_Status（不用 Cell_Type，带 scores）---------
set.seed(123)
rf_disease_basic <- randomForest(
  Disease_Status ~ Gene_E_Housekeeping +
    Gene_A_Oncogene +
    Gene_B_Immune +
    Gene_C_Stromal +
    Gene_D_Therapy +
    Pathway_Score_Inflam +
    Tumor_score +
    Immune_score +
    Stromal_score,
  data       = train_df,
  importance = TRUE,
  ntree      = 500
)

print(rf_disease_basic)

## 6. 随机森林：预测 Disease_Status（加入 Cell_Type） ---------------------
set.seed(123)
rf_disease_ct <- randomForest(
  Disease_Status ~ Gene_E_Housekeeping +
    Gene_A_Oncogene +
    Gene_B_Immune +
    Gene_C_Stromal +
    Gene_D_Therapy +
    Pathway_Score_Inflam +
    Cell_Type,
  data       = train_df,
  importance = TRUE,
  ntree      = 500
)

print(rf_disease_ct)

## 7. 保存模型到 .rds ---------------------------------------------------
saveRDS(rf_celltype,      file = "rf_celltype.rds")
saveRDS(rf_disease_basic, file = "rf_disease_basic.rds")
saveRDS(rf_disease_ct,    file = "rf_disease_ct.rds")

## 8. 生成 Shiny 用 example 输入 CSV（只有特征列，没有标签）------------
set.seed(123)
input_cols <- c(
  "Gene_E_Housekeeping",
  "Gene_A_Oncogene",
  "Gene_B_Immune",
  "Gene_C_Stromal",
  "Gene_D_Therapy",
  "Pathway_Score_Inflam"
)

example_idx <- sample.int(nrow(data_df), size = 100)
example_df  <- data_df[example_idx, input_cols]

write.csv(example_df, file = "gene_input_example.csv", row.names = FALSE)

cat("\nOK：已生成这些文件：\n",
    "- rf_celltype.rds\n",
    "- rf_disease_basic.rds\n",
    "- rf_disease_ct.rds\n",
    "- gene_input_example.csv\n")
