## ===== Global setup: 数据 & 模型 &工具函数（放在 ui 之前） =====

options(shiny.usecairo = TRUE)       # 文字更清晰
options(stringsAsFactors = FALSE)

library(shiny)
library(shinythemes)
library(DT)
library(dplyr)
library(ggplot2)
library(randomForest)
library(shinydashboard)
library(httr)
library(jsonlite)
if (file.exists(".Renviron")) {
  readRenviron(".Renviron")
}
# ====== Gemini API 调用函数：输入 prompt，返回一段英文文本 ======
call_gemini_tip <- function(prompt_text,
                            model = "gemini-2.5-flash") {
  api_key <- Sys.getenv("GEMINI_API_KEY")
  if (identical(api_key, "")) {
    stop("GEMINI_API_KEY is not set in your environment.")
  }
  
  url <- paste0(
    "https://generativelanguage.googleapis.com/v1beta/models/",
    model,
    ":generateContent?key=",
    api_key
  )
  
  body <- list(
    contents = list(
      list(
        role = "user",
        parts = list(
          list(text = prompt_text)
        )
      )
    )
  )
  
  # 发送请求
  res <- httr::POST(
    url,
    body   = jsonlite::toJSON(body, auto_unbox = TRUE),
    encode = "raw",
    httr::add_headers("Content-Type" = "application/json")
  )
  
  status  <- httr::status_code(res)
  raw_txt <- httr::content(res, as = "text", encoding = "UTF-8")
  
  message("Gemini HTTP status: ", status)
  # 调试用：可以先保留，看完一次结构后再注释掉
  message("Gemini raw JSON (first 500 chars):\n",
          substr(raw_txt, 1, 500), if (nchar(raw_txt) > 500) " ...")
  
  if (status >= 400) {
    stop("HTTP error from Gemini API: status = ", status, "\nBody: ", raw_txt)
  }
  
  # 这里用 simplifyVector = FALSE，避免自动变成 data.frame 导致结构奇怪
  j <- jsonlite::fromJSON(raw_txt, simplifyVector = FALSE)
  
  if (!is.null(j$error)) {
    stop("Gemini API error: ", j$error$message)
  }
  if (is.null(j$candidates) || length(j$candidates) == 0) {
    stop("Gemini returned no candidates.")
  }
  
  cand <- j$candidates[[1]]
  
  # 递归地把所有 text 字段拿出来
  extract_texts <- function(x) {
    out <- character()
    if (is.list(x)) {
      # 如果自己就有 text 字段
      if (!is.null(x$text) && is.character(x$text)) {
        out <- c(out, x$text)
      }
      # 继续往子节点递归
      for (el in x) {
        out <- c(out, extract_texts(el))
      }
    }
    out
  }
  
  texts <- extract_texts(cand)
  texts <- texts[!is.na(texts)]
  
  if (!length(texts)) {
    stop("Gemini returned candidates but no 'text' fields.")
  }
  
  # 通常第一个或者最长的就是我们要的内容，这里取最长那段
  txt <- texts[which.max(nchar(texts))]
  
  # 去掉两侧空白
  trimws(txt)
}


# ====== 把单样本预测结果整理成一段 summary（给 prompt 用） ======
make_summary_for_api <- function(df_row) {
  paste0(
    "Cell type: ", df_row$Pred_Cell_Type, ". ",
    "Disease_Status prediction: ", df_row$Pred_Disease,
    " with p(Tumor) = ", round(df_row$Prob_Tumor, 3), ". ",
    "Scores: Tumor_score = ", round(df_row$Tumor_score, 2), ", ",
    "Immune_score = ", round(df_row$Immune_score, 2), ", ",
    "Stromal_score = ", round(df_row$Stromal_score, 2), ". ",
    "This dataset is fully synthetic and used only for teaching in a BIO215 project."
  )
}

## 1) 读入完整数据集 ------------------------------------------------------

# 每行一个“细胞”
data_df <- read.csv(
  "Gene_Expression_Analysis_and_Disease_Relationship_Synthetic.csv",
  stringsAsFactors = FALSE
)

# 目标变量设定因子顺序（和 ML 结果一致）
data_df$Cell_Type <- factor(
  data_df$Cell_Type,
  levels = c("Cancer", "T_Cell", "Fibroblast")
)

data_df$Disease_Status <- factor(
  data_df$Disease_Status,
  levels = c("Healthy_Control", "Tumor")
)

## 2) 特征工程：三个 score（后面模型 & 可视化 & tip 都会用到） -----------

add_scores <- function(df) {
  df$Tumor_score   <- df$Gene_A_Oncogene + df$Gene_D_Therapy
  df$Immune_score  <- df$Gene_B_Immune   + df$Pathway_Score_Inflam
  df$Stromal_score <- df$Gene_C_Stromal
  df
}

data_df <- add_scores(data_df)

# 机器学习特征列（输入给模型的基础六个基因/通路）
feature_cols <- c(
  "Gene_E_Housekeeping",
  "Gene_A_Oncogene",
  "Gene_B_Immune",
  "Gene_C_Stromal",
  "Gene_D_Therapy",
  "Pathway_Score_Inflam"
)

# 在 data_df <- add_scores(data_df) 之后
num_cols_dt <- c(
  "Gene_E_Housekeeping",
  "Gene_A_Oncogene",
  "Gene_B_Immune",
  "Gene_C_Stromal",
  "Gene_D_Therapy",
  "Pathway_Score_Inflam",
  "Tumor_score",
  "Immune_score",
  "Stromal_score",
  "UMAP_1"
)

gene_vars <- c(
  feature_cols,            # 6 个原始特征
  "Tumor_score",
  "Immune_score",
  "Stromal_score"
)

## 3) 为 Gene sandbox 准备 slider 范围（整数+少量小数） -------------------

# 这里直接 floor / ceiling，避免一大串小数
gene_ranges <- data_df %>%
  summarise(
    Gene_A_min = floor(min(Gene_A_Oncogene,        na.rm = TRUE)),
    Gene_A_max = ceiling(max(Gene_A_Oncogene,      na.rm = TRUE)),
    Gene_B_min = floor(min(Gene_B_Immune,          na.rm = TRUE)),
    Gene_B_max = ceiling(max(Gene_B_Immune,        na.rm = TRUE)),
    Gene_C_min = floor(min(Gene_C_Stromal,         na.rm = TRUE)),
    Gene_C_max = ceiling(max(Gene_C_Stromal,       na.rm = TRUE)),
    Gene_D_min = floor(min(Gene_D_Therapy,         na.rm = TRUE)),
    Gene_D_max = ceiling(max(Gene_D_Therapy,       na.rm = TRUE)),
    Gene_E_min = floor(min(Gene_E_Housekeeping,    na.rm = TRUE)),
    Gene_E_max = ceiling(max(Gene_E_Housekeeping,  na.rm = TRUE)),
    Path_min   = floor(min(Pathway_Score_Inflam,   na.rm = TRUE)),
    Path_max   = ceiling(max(Pathway_Score_Inflam, na.rm = TRUE))
  )

# 之后在 ui 里可以这样用：
# sliderInput("geneA", ..., min = gene_ranges$Gene_A_min, max = gene_ranges$Gene_A_max, step = 0.1)

## 4) 读入训练好的 RF 模型 + 示例输入文件 -------------------------------

rf_celltype      <- readRDS("rf_celltype.rds")
rf_disease_basic <- readRDS("rf_disease_basic.rds")
rf_disease_ct    <- readRDS("rf_disease_ct.rds")

# 提供给用户下载 / 作为 Multiple prediction 的模板
example_input <- read.csv("gene_input_example.csv", stringsAsFactors = FALSE)

## 5) 公共小函数：预处理 & 预测 & 文本 tip -------------------------------

# 统一预处理：保证列存在、转 numeric、补充 scores
prep_new_samples <- function(df) {
  num_cols <- c(
    "Gene_E_Housekeeping",
    "Gene_A_Oncogene",
    "Gene_B_Immune",
    "Gene_C_Stromal",
    "Gene_D_Therapy",
    "Pathway_Score_Inflam"
  )
  missing <- setdiff(num_cols, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }
  df[num_cols] <- lapply(df[num_cols], as.numeric)
  df <- add_scores(df)
  df
}

# 多样本：预测 Cell_Type（返回原数据 + 概率）
predict_celltype_df <- function(df, model = rf_celltype) {
  X <- prep_new_samples(df)
  probs <- predict(model, newdata = X, type = "prob")
  preds <- predict(model, newdata = X, type = "class")
  out   <- cbind(X, Pred_Cell_Type = preds, probs)
  rownames(out) <- NULL
  out
}

# 多样本：疾病预测（只用基因+scores）
predict_disease_basic_df <- function(df, model = rf_disease_basic, cutoff = 0.5) {
  X <- prep_new_samples(df)
  probs     <- predict(model, newdata = X, type = "prob")
  tumor_prob <- probs[, "Tumor"]
  pred_label <- ifelse(tumor_prob >= cutoff, "Tumor", "Healthy_Control")
  out <- data.frame(
    Pred_Tumor_Prob_basic = tumor_prob,
    Pred_Label_basic      = factor(pred_label, levels = c("Healthy_Control", "Tumor"))
  )
  out
}

# 多样本：疾病预测（两步：用 Cell_Type + 基因+scores）
predict_disease_ct_df <- function(df,
                                  model_ct  = rf_celltype,
                                  model_dis = rf_disease_ct,
                                  cutoff = 0.5) {
  X  <- prep_new_samples(df)
  
  # 如果输入没有 Cell_Type，就先用 model_ct 预测一个
  if (!"Cell_Type" %in% names(X)) {
    ct_pred        <- predict(model_ct, newdata = X, type = "class")
    X$Cell_Type    <- factor(ct_pred,
                             levels = c("Cancer", "T_Cell", "Fibroblast"))
  }
  
  probs      <- predict(model_dis, newdata = X, type = "prob")
  tumor_prob <- probs[, "Tumor"]
  pred_label <- ifelse(tumor_prob >= cutoff, "Tumor", "Healthy_Control")
  
  out <- data.frame(
    Pred_Tumor_Prob_ct = tumor_prob,
    Pred_Label_ct      = factor(pred_label, levels = c("Healthy_Control", "Tumor"))
  )
  out
}

# 总控：多样本预测（给 Shiny multiple upload 用）
# feature_df: 至少包含六个基因列
# dz_model : "basic"（只用基因+scores）或 "ct"（两步）
prediction_multiple <- function(feature_df,
                                dz_model = c("basic", "ct"),
                                positive_threshold = 0.5) {
  dz_model <- match.arg(dz_model)
  
  # 基础预处理 + scores
  X <- prep_new_samples(feature_df)
  
  # Step 1: 总是先预测 Cell_Type
  ct_probs <- predict(rf_celltype, newdata = X, type = "prob")
  ct_lab   <- colnames(ct_probs)[max.col(ct_probs, ties.method = "first")]
  
  # Step 2: 预测 Disease_Status
  if (dz_model == "basic") {
    dz_probs <- predict(rf_disease_basic, newdata = X, type = "prob")
    model_used <- "rf_disease_basic (genes + scores)"
  } else {
    X_ct             <- X
    X_ct$Cell_Type   <- factor(ct_lab, levels = levels(data_df$Cell_Type))
    dz_probs         <- predict(rf_disease_ct, newdata = X_ct, type = "prob")
    model_used       <- "rf_disease_ct (genes + scores + predicted Cell_Type)"
  }
  
  p_tumor  <- dz_probs[, "Tumor"]
  pred_dz  <- ifelse(p_tumor >= positive_threshold, "Tumor", "Healthy_Control")
  
  out_df <- cbind(
    feature_df,
    Tumor_score   = X$Tumor_score,
    Immune_score  = X$Immune_score,
    Stromal_score = X$Stromal_score,
    Pred_Cell_Type = ct_lab,
    Pred_Disease   = factor(pred_dz, levels = c("Healthy_Control", "Tumor")),
    Prob_Tumor     = p_tumor,
    Model_Disease_Used = model_used,
    stringsAsFactors = FALSE
  )
  
  rownames(out_df) <- NULL
  out_df
}

# 单样本预测（给 Shiny 单个 slider / numericInput 用）
prediction_single <- function(gene_E, gene_A, gene_B, gene_C, gene_D, path_inflam,
                              dz_model = c("basic", "ct"),
                              positive_threshold = 0.5) {
  df1 <- data.frame(
    Gene_E_Housekeeping = gene_E,
    Gene_A_Oncogene     = gene_A,
    Gene_B_Immune       = gene_B,
    Gene_C_Stromal      = gene_C,
    Gene_D_Therapy      = gene_D,
    Pathway_Score_Inflam = path_inflam
  )
  prediction_multiple(
    feature_df        = df1,
    dz_model          = dz_model,
    positive_threshold = positive_threshold
  )
}

# 风险分级（给 UI 上的小标签用：High / Intermediate / Low）
risk_label <- function(prob) {
  if (is.na(prob)) return(NA_character_)
  if (prob >= 0.7) {
    "High (simulated)"
  } else if (prob >= 0.3) {
    "Intermediate (simulated)"
  } else {
    "Low (simulated)"
  }
}

# 教学用 tip：根据疾病标签 + Immune_score + p_tumor 生成一段英文说明
generate_tip <- function(disease_label, immune_score, p_tumor) {
  immune_q75 <- quantile(data_df$Immune_score, 0.75)
  
  if (is.na(p_tumor)) {
    return("No valid probability is available for this sample.")
  }
  
  if (disease_label == "Tumor" && immune_score > immune_q75) {
    paste0(
      "Predicted Tumor-like expression with a relatively high Immune_score.\n",
      "In a real research setting, such profiles could motivate further study of ",
      "tumour–immune interactions and potential response to immune-modulating ",
      "strategies.\nAll results here are based on a fully synthetic dataset and ",
      "machine learning models for teaching only. They must not be used for any ",
      "real clinical decision making."
    )
  } else if (disease_label == "Tumor") {
    paste0(
      "Predicted Tumor-like expression pattern with p(Tumor) ≈ ",
      round(p_tumor, 3), ".\n",
      "In realistic applications, such a model would only provide a screening ",
      "signal and must be confirmed by standard clinical tests and specialist ",
      "interpretation.\n\nThis Shiny app is purely educational and is not a ",
      "medical device."
    )
  } else {
    paste0(
      "Predicted Healthy_Control-like pattern with p(Tumor) ≈ ",
      round(p_tumor, 3), ".\n",
      "Even if a sample is predicted as non-tumour by a model, real-world health ",
      "decisions require comprehensive clinical assessment and regular check-ups.\n\n",
      "The current dataset and models are synthetic BIO215 learning materials ",
      "only and should not be used to guide personal health behaviour."
    )
  }
}



# ===== UI 部分：导航栏 + 各个 tabPanel / navbarMenu =====
ui <- navbarPage(
  title = "TME Cancer Microenvironment Explorer",
  theme = shinytheme("flatly"),
  id    = "main_nav",
  header = tags$head(
    tags$link(rel="stylesheet", href="https://cdnjs.cloudflare.com/ajax/libs/admin-lte/2.4.3/css/AdminLTE.min.css"),
    
    tags$style(HTML("
      /* 顶部 navbar 颜色（深蓝） */
      .navbar-default {
        background-color: #2c3e50;
        border-color: #2c3e50;
      }
      .navbar-default .navbar-brand,
      .navbar-default .navbar-nav > li > a {
        color: #ecf0f1;
        font-weight: 500;
      }
      .navbar-default .navbar-brand:hover,
      .navbar-default .navbar-nav > li > a:hover {
        color: #ffffff;
      }

      .navbar-default .navbar-nav > .active > a,
      .navbar-default .navbar-nav > .active > a:focus,
      .navbar-default .navbar-nav > .active > a:hover,
      .navbar-default .navbar-nav > .open > a,
      .navbar-default .navbar-nav > .open > a:focus,
      .navbar-default .navbar-nav > .open > a:hover {
        background-color: #34495e;
        color: #ffffff;
      }

      .navbar-default .navbar-nav > li > a {
        width: 220px;
        text-align: center;
        white-space: nowrap;
      }

      .navbar-nav > li > .dropdown-menu {
        min-width: 220px;
        margin-top: 0;
        padding: 0;
        border-radius: 0 0 4px 4px;
      }
      .navbar-nav > li > .dropdown-menu > li > a {
        padding: 8px 16px;
        color: #2c3e50;
      }
      .navbar-nav > li > .dropdown-menu > li > a:hover {
        background-color: #ecf0f1;
        color: #2c3e50;
      }

      /* ---- 字号 & 标题 ---- */
      body {
        font-size: 15px;          /* 正文稍微大一点 */
      }
      h2 {
        margin-top: 6px;
        margin-bottom: 14px;
        font-weight: 800;
        font-size: 28px;          /* Welcome 这种标题更大 */
      }
      h3 {
        font-weight: 800;
        font-size: 22px;          /* Project workflow / Feature description */
      }
      h4 {
        font-weight: 700;
        font-size: 18px;
      }
      h5 {
        font-weight: 600;
        font-size: 15px;
      }

     
      /* 子标题：Background / Cell type / Disease */
      .overview-subtitle {
        font-weight: 700;
        font-size: 16px;
        margin-bottom: 6px;
      }

      /* Dataset overview 里的两个小表：统一白底 + 黑边框 */
      #summary_ct table,
      #summary_ct th,
      #summary_ct td,
      #summary_ds table,
      #summary_ds th,
      #summary_ds td {
      
        border: 1px solid #000000;
        border-collapse: collapse;
      }
      #summary_ct th,
      #summary_ct td,
      #summary_ds th,
      #summary_ds td {
        background-color: #ffffff;
        padding: 4px 8px;
        font-size: 12px;  
      }



      /* 下拉菜单里【当前选中的】那一项：深蓝背景 + 白字 */
      .navbar-nav > li > .dropdown-menu > li.active > a,
      .navbar-nav > li > .dropdown-menu > li.active > a:focus,
      .navbar-nav > li > .dropdown-menu > li.active > a:hover {
        background-color: #34495e;
        color: #ffffff;
      }
      
            .overview-subtitle {
        font-size: 17px;
        font-weight: 700;
        color: #111827;
        
        
              /* ===== valueBox 外观：仿 m6A 样式 ===== */
      .small-box {
        position: relative;
        display: block;
        border-radius: 4px;
        margin-bottom: 12px;
        box-shadow: 0 1px 2px rgba(15, 23, 42, 0.15);
        color: #ffffff;
      }

      .small-box > .inner {
        position: relative;
        z-index: 1;
        padding: 10px 12px;
      }

      .small-box h3 {
        font-size: 30px;
        font-weight: 700;
        margin: 0 0 4px 0;
        white-space: nowrap;
      }

      .small-box p {
        font-size: 14px;
        margin: 0;
        opacity: 0.98;
      }

      /* icon 放在右上角，稍微透明一点 */
      .small-box .icon {
        position: absolute;
        top: 6px;
        right: 12px;
        z-index: 0;
        font-size: 40px;
        opacity: 0.35;
      }

      /* 三种背景色：light-blue / teal / aqua */
      .small-box.bg-light-blue {
        background-color: #3B82F6;
      }
      .small-box.bg-teal {
        background-color: #14B8A6;
      }
      .small-box.bg-aqua {
        background-color: #06B6D4;
      }

      }
      
/* ===============================
   全站统一 hero 风格（覆盖所有标题）
   =============================== */

/* 全站背景字稍微亮一点 */

/* h1 = hero h1 */
h1 {
  font-size: 34px;
  font-weight: 900;
  letter-spacing: .02em;
  margin-top: 6px;
  margin-bottom: 10px;
}

/* h2 = hero h2 */
h2 {
  font-size: 28px;
  font-weight: 800;
  margin-top: 6px;
  margin-bottom: 12px;
}

/* h3 = hero h3 */
h3 {
  font-size: 22px;
  font-weight: 800;
  margin-top: 6px;
  margin-bottom: 10px;
}

/* h4 = hero h4 */
h4 {
  font-size: 18px;
  font-weight: 700;
  margin-top: 6px;
  margin-bottom: 8px;
}

/* h5 = hero h5 */
h5 {
  font-size: 16px;
  font-weight: 650;
  margin-top: 4px;
  margin-bottom: 6px;
}

/* p = hero p */
p, .well p, .box p {
  font-size: 16px;
  line-height: 1.45;
  margin-bottom: 10px;
}

/* Home：右侧“使用步骤”列表更紧凑 + 字号和 Background 一样偏小 */
.use-steps li {
  margin-bottom: 4px;
  font-size: 13.5px;   /* 和 Background 那段 p(style=13.5px) 对齐 */
  line-height: 1.4;
}


/* overview 字体也同步放大 */
.overview-subtitle {
  font-weight: 700;
  font-size: 17px;
  margin-bottom: 6px;
}


.ct-prob-row {
  margin-bottom: 6px;
}
.ct-prob-label {
  font-weight: 600;
  margin-bottom: 2px;
}

#home-page h4 {
  font-size: 15px;   /* 原来是 18，这样 Background 和右边小标题都会矮一点 */
  margin-bottom: 6px;
}
#home-page h5 {
  font-size: 12px;   /* 原来是 18，这样 Background 和右边小标题都会矮一点 */
  margin-bottom: 6px;
}

/* Home 右侧三个 wellPanel 之间的竖向间距缩小 */
#home-right-col .well {
  margin-bottom: 8px;   /* 默认大概 20px → 改小 */
  padding-top: 10px;    /* 也顺便让内部上下不要太厚 */
  padding-bottom: 10px;
}

/* 右侧 h4 标题也稍微紧凑一点（如果你还没加过） */
#home-right-col h4 {
  margin-top: 4px;
  margin-bottom: 6px;
}

    "))
  ),
  
  
  
  # ===== DATA & VISUALIZATION =====
  navbarMenu(
    title = tagList(icon("database"), "DATA & VISUALIZATION"),
    
    tabPanel(
      title = tagList(icon("home"), "Home"),
      fluidPage(
        div(
          id = "home-page",
          
          br(),
          h1("Welcome - TME Cancer Microenvironment Explorer"),
          p(
            "This Shiny application explores a synthetic single-cell RNA-seq dataset ",
            "that mimics the tumor microenvironment of a solid tumor. Each row ",
            "represents one cell with gene expression features, cell type annotations, ",
            "and a disease status label."
          ),
          
          hr(),
          h3("Project workflow"),
          p(
            "The figure on the left summarizes the overall analysis and app workflow. ",
            "On the right, you can see a brief dataset background, overview, and how to use this app."
          ),
          
          # ===== 中间：左图 + 右侧 overview + how-to =====
          fluidRow(
            # 左：流程图图片
            column(
              width = 8,
              wellPanel(
                imageOutput("workflow_img", height = "auto")
              )
            ),
            
            # 右：overview + background + how to
            column(
              width = 4,
              tags$div(
                id = "home-right-col", 
                # Part 1: Dataset overview
                wellPanel(
                  h4("Dataset overview"),
                  verbatimTextOutput("summary_ncells")
                ),
                
                # ⭐ 去掉这里原来的 br()，减少多余空白
                
                # Part 2: Background + Cell / Disease 分布 一块
                wellPanel(
                  h4("Background"),
                  p(
                    style = "font-size: 13.5px;",
                    "This app uses a synthetic single-cell RNA-seq-style dataset that mimics ",
                    "the tumour microenvironment (TME) of a solid tumour. It was generated ",
                    "for educational and benchmarking purposes, with clean gene-expression ",
                    "profiles, clear cell-type labels (Cancer, T_Cell, Fibroblast) and a ",
                    "Disease_Status label (Tumor vs Healthy_Control)."
                  ),
                  
                  tags$hr(style = "margin-top: 4px; margin-bottom: 4px;"),
                  
                  fluidRow(
                    column(
                      width = 6,
                      h5("Cell type distribution"),
                      tableOutput("summary_ct")
                    ),
                    column(
                      width = 6,
                      h5("Disease status distribution"),
                      tableOutput("summary_ds")
                    )
                  )
                ),
                
                # ⭐ 这里也把原来的 br() 删掉
                
                # Part 3: How to use
                wellPanel(
                  h4("How to use this app"),
                  tags$ol(
                    class = "use-steps",
                    tags$li("Start on this Home page to understand the dataset and key features."),
                    tags$li("Use the Data / Visualize tabs to explore the synthetic microenvironment."),
                    tags$li("Use the Prediction tabs for cell type and disease status predictions."),
                    tags$li("Try the Gene sandbox to perform simple 'what-if' simulations.")
                  )
                )
              )
            )
          ),
          
          
          
          # ===== 单独一行的 Feature description（整行 12 列） =====
          fluidRow(
            column(
              width = 12,
              wellPanel(
                h3("Feature description"),
                tags$ul(
                  tags$li(
                    tags$b("Gene_A_Oncogene: "),
                    "simulated oncogene-like marker; higher values tend to be associated with tumor-like cells."
                  ),
                  tags$li(
                    tags$b("Gene_B_Immune: "),
                    "immune-related marker; higher values resemble immune-infiltrating cells such as T cells."
                  ),
                  tags$li(
                    tags$b("Gene_C_Stromal: "),
                    "stromal marker; captures fibroblast-like microenvironmental cells."
                  ),
                  tags$li(
                    tags$b("Gene_D_Therapy: "),
                    "therapy-related marker; conceptually represents potential drug targets."
                  ),
                  tags$li(
                    tags$b("Gene_E_Housekeeping: "),
                    "housekeeping-like gene used as a reference for overall expression level."
                  ),
                  tags$li(
                    tags$b("Pathway_Score_Inflam: "),
                    "inflammatory pathway activity score summarizing multiple genes."
                  ),
                  tags$li(
                    tags$b("Tumor_score: "),
                    "composite score based on Gene_A_Oncogene and Gene_D_Therapy."
                  ),
                  tags$li(
                    tags$b("Immune_score: "),
                    "composite score based on Gene_B_Immune and Pathway_Score_Inflam."
                  ),
                  tags$li(
                    tags$b("Stromal_score: "),
                    "composite score based on Gene_C_Stromal."
                  )
                )
              )
            )
          ),
          
          hr(),
          em(
            "Note: This is a fully synthetic dataset created for BIO215 learning. ",
            "Any 'risk' information shown in the app is simulated and not intended for clinical use."
          )
        ) # div(id = "home-page") 结束
      )
    )
    ,
    
    # 其余 Data / Visualize / Download 先保持占位
    # ---- Data tab: interactive table + filters + download ----
    tabPanel(
      title = tagList(icon("table"), "Data"),
      fluidPage(
        br(),
        fluidRow(
          # 左侧：筛选 + 下载设置
          column(
            width = 3,
            h3("Filter & download"),
            p("Use the controls below to subset the synthetic scRNA-seq dataset, ",
              "inspect the table, and download subsets for further analysis."),
            
            h4("Filter dataset"),
            selectInput(
              "data_filter_disease", "Disease_Status:",
              choices  = c("All", levels(data_df$Disease_Status)),
              selected = "All"
            ),
            selectInput(
              "data_filter_celltype", "Cell_Type:",
              choices  = c("All", levels(data_df$Cell_Type)),
              selected = "All"
            ),
            
            hr(),
            h4("Download options"),
            radioButtons(
              "data_dl_scope", "Download scope:",
              choices  = c("Filtered rows" = "filtered",
                           "Selected rows" = "selected"),
              selected = "filtered",
              inline   = FALSE
            ),
            downloadButton("data_dl_csv", "Download CSV"),
            br(), br(),
            actionButton("data_clear_sel", "Clear selection")
          ),
          
          # 右侧：简介 + summary + DT 表
          column(
            width = 9,
            h2("TME mini bioinformatics database"),
            p(
              "This tab treats the synthetic tumour microenvironment scRNA-seq dataset as a ",
              "mini bioinformatics database. After applying the filters, you can inspect the ",
              "rows in the interactive table and optionally export them as CSV."
            ),
            
            # 三个 valueBox
            fluidRow(
              valueBoxOutput("dt_rows",    width = 4),
              valueBoxOutput("dt_cols",    width = 4),
              valueBoxOutput("dt_posrate", width = 4)
            ),
            
            br(),
            
            DTOutput("table_data")
          )
          
        )
      )
    ),
    
    tabPanel(
      title = tagList(icon("chart-pie"), "Visualize"),
      fluidPage(
        br(),
        h2("Visualize the synthetic TME dataset"),
        p(
          "This tab summarizes the synthetic tumour microenvironment dataset with ",
          "a donut-style pie chart, bar chart, boxplot and histogram."
        ),
        br(),
        
        # ===== Row 1: Donut + Bar =====
        fluidRow(
          box(
            width = 6, status = "primary", solidHeader = TRUE,
            title = tagList(icon("chart-pie"), "Label distribution"),
            div(
              style = "display:flex; align-items:center; gap:12px; margin-bottom:8px;",
              tags$label(
                "Categorical variable:",
                `for` = "pie_col",
                style = "margin:0; font-weight:600;"
              ),
              selectInput(
                "pie_col", NULL,
                choices  = c("Cell type" = "Cell_Type",
                             "Disease status" = "Disease_Status"),
                selected = "Cell_Type",
                width    = "220px"
              )
            ),
            plotOutput("pie_plot", height = 320)
          ),
          
          box(
            width = 6, status = "info", solidHeader = TRUE,
            title = tagList(icon("chart-bar"), "Gene / score summary by group"),
            fluidRow(
              column(
                6,
                selectInput(
                  "bar_gene", "Gene / score:",
                  choices  = gene_vars,
                  selected = "Gene_A_Oncogene"
                )
              ),
              column(
                6,
                selectInput(
                  "bar_group", "Group by:",
                  choices  = c("Cell_Type", "Disease_Status"),
                  selected = "Cell_Type"
                )
              )
            ),
            radioButtons(
              "bar_stat", "Statistic:",
              choices  = c("Mean" = "mean", "Median" = "median"),
              selected = "mean",
              inline   = TRUE
            ),
            plotOutput("bar_plot", height = 320)
          )
        ),
        
        # ===== Row 2: Boxplot + Histogram =====
        fluidRow(
          box(
            width = 6, status = "warning", solidHeader = TRUE,
            title = tagList(icon("box-open"), "Distribution by group"),
            fluidRow(
              column(
                6,
                selectInput(
                  "box_gene", "Gene / score:",
                  choices  = gene_vars,
                  selected = "Gene_A_Oncogene"
                )
              ),
              column(
                6,
                selectInput(
                  "box_group", "Group by:",
                  choices  = c("Cell_Type", "Disease_Status"),
                  selected = "Disease_Status"
                )
              )
            ),
            plotOutput("box_plot", height = 320)
          ),
          
          box(
            width = 6, status = "success", solidHeader = TRUE,
            title = tagList(icon("chart-area"), "Histogram of gene / score"),
            fluidRow(
              column(
                6,
                selectInput(
                  "hist_gene", "Gene / score:",
                  choices  = gene_vars,
                  selected = "Gene_A_Oncogene"
                )
              ),
              column(
                3,
                selectInput(
                  "hist_disease", "Disease_Status:",
                  choices  = c("All", levels(data_df$Disease_Status)),
                  selected = "All"
                )
              ),
              column(
                3,
                selectInput(
                  "hist_celltype", "Cell_Type:",
                  choices  = c("All", levels(data_df$Cell_Type)),
                  selected = "All"
                )
              )
            ),
            sliderInput(
              "hist_bins", "Bins:",
              min = 10, max = 80, value = 30, step = 5
            ),
            plotOutput("hist_plot", height = 320)
          )
        )
      )
    ),
    
    tabPanel(
      title = tagList(icon("download"), "Download"),
      fluidPage(
        br(),
        h2("Download custom subsets"),
        p(
          "Use this tab to build and export custom subsets of the synthetic tumour ",
          "microenvironment (TME) dataset. You can filter by disease status and cell ",
          "type, choose which columns to keep, preview the subset, and download it ",
          "as a CSV file for downstream analysis."
        ),
        br(),
        
        fluidRow(
          # ===== 左侧：筛选 + 选列 + 下载 =====
          box(
            width = 3, status = "primary", solidHeader = TRUE,
            title = "Filters and columns",
            
            h4("Row filters"),
            checkboxGroupInput(
              "down_disease", "Disease_Status:",
              choices  = levels(data_df$Disease_Status),
              selected = levels(data_df$Disease_Status)
            ),
            checkboxGroupInput(
              "down_celltype", "Cell_Type:",
              choices  = levels(data_df$Cell_Type),
              selected = levels(data_df$Cell_Type)
            ),
            
            tags$hr(),
            h4("Columns to include"),
            checkboxGroupInput(
              "down_cols", NULL,
              choices  = names(data_df),
              selected = names(data_df)
            ),
            
            tags$hr(),
            h4("Download subset"),
            helpText(
              "The CSV file will contain the rows and columns defined above. ",
              "Preview on the right updates in real time."
            ),
            downloadButton("down_btn", "Download subset (.csv)")
          ),
          
          # ===== 右侧：summary + 预览表 =====
          box(
            width = 9, status = "info", solidHeader = TRUE,
            title = "Preview of current subset",
            
            fluidRow(
              valueBoxOutput("down_nrows",  width = 4),
              valueBoxOutput("down_ncols",  width = 4),
              valueBoxOutput("down_tumor",  width = 4)   # 新增
            ),
            br(),
            h4("First 50 rows"),
            DTOutput("down_table"),
            br(),
            h4("Download ML-ready feature table"),
            helpText("Rows follow the same filters, columns are fixed: feature genes, scores and labels."),
            downloadButton("down_btn_ml", "Download features + scores (.csv)")
          )
        )
      )
    )
  ),
  
  # ===== PREDICTION =====
  navbarMenu(
    title = tagList(icon("project-diagram"), "PREDICTION"),
    tabPanel(
      title = tagList(icon("th-list"), "Cell type prediction"),
      fluidPage(
        br(),
        h2("Cell type prediction (Random Forest)"),
        p(
          "This module deploys the trained Random Forest (rf_celltype) to predict ",
          "cell type (Cancer / T_Cell / Fibroblast) from six gene / pathway features."
        ),
        br(),
        
        tabsetPanel(
          id   = "ct_tabs",
          type = "tabs",
          
          # ---------- Multiple (CSV) ----------
          tabPanel(
            title = tagList(icon("layer-group"), "Multiple (CSV)"),
            br(),
            fluidRow(
              # 左：上传 + 设置
              column(
                width = 4,
                box(
                  width = 12, status = "primary", solidHeader = TRUE,
                  title = "Step 1. Upload CSV or use example",
                  
                  fileInput(
                    "ct_file",
                    "Upload CSV with feature columns",
                    accept = ".csv"
                  ),
                  checkboxInput(
                    "ct_use_example",
                    "Use example input from gene_input_example.csv",
                    value = FALSE
                  ),
                  tags$hr(),
                  h4("Required columns in the input CSV:"),
                  tags$code(
                    "Gene_E_Housekeeping, Gene_A_Oncogene, Gene_B_Immune, ",
                    "Gene_C_Stromal, Gene_D_Therapy, Pathway_Score_Inflam"
                  ),
                  tags$hr(),
                  actionButton(
                    "ct_run_multi",
                    "Run cell type prediction",
                    icon = icon("play")
                  ),
                  br(), br(),
                  downloadButton(
                    "ct_download_multi",
                    "Download predictions (.csv)"
                  )
                )
              ),
              # 右：结果表
              column(
                width = 8,
                box(
                  width = 12, status = "info", solidHeader = TRUE,
                  title = "Prediction results",
                  p(
                    "The table below shows the original input features together with ",
                    "engineered scores, predicted cell type, and class probabilities."
                  ),
                  DTOutput("ct_multi_table")
                )
              )
            )
          ),
          
          # ---------- Single sample ----------
          tabPanel(
            title = tagList(icon("vial"), "Single sample"),
            br(),
            fluidRow(
              # 左：输入
              column(
                width = 4,
                box(
                  width = 12, status = "primary", solidHeader = TRUE,
                  title = "Input features",
                  
                  numericInput(
                    "ct_single_E", "Gene_E_Housekeeping",
                    value = round((gene_ranges$Gene_E_min + gene_ranges$Gene_E_max)/2, 1),
                    min   = gene_ranges$Gene_E_min,
                    max   = gene_ranges$Gene_E_max,
                    step  = 0.1
                  ),
                  numericInput(
                    "ct_single_A", "Gene_A_Oncogene",
                    value = round((gene_ranges$Gene_A_min + gene_ranges$Gene_A_max)/2, 1),
                    min   = gene_ranges$Gene_A_min,
                    max   = gene_ranges$Gene_A_max,
                    step  = 0.1
                  ),
                  numericInput(
                    "ct_single_B", "Gene_B_Immune",
                    value = round((gene_ranges$Gene_B_min + gene_ranges$Gene_B_max)/2, 1),
                    min   = gene_ranges$Gene_B_min,
                    max   = gene_ranges$Gene_B_max,
                    step  = 0.1
                  ),
                  numericInput(
                    "ct_single_C", "Gene_C_Stromal",
                    value = round((gene_ranges$Gene_C_min + gene_ranges$Gene_C_max)/2, 1),
                    min   = gene_ranges$Gene_C_min,
                    max   = gene_ranges$Gene_C_max,
                    step  = 0.1
                  ),
                  numericInput(
                    "ct_single_D", "Gene_D_Therapy",
                    value = round((gene_ranges$Gene_D_min + gene_ranges$Gene_D_max)/2, 1),
                    min   = gene_ranges$Gene_D_min,
                    max   = gene_ranges$Gene_D_max,
                    step  = 0.1
                  ),
                  numericInput(
                    "ct_single_path", "Pathway_Score_Inflam",
                    value = round((gene_ranges$Path_min + gene_ranges$Path_max)/2, 1),
                    min   = gene_ranges$Path_min,
                    max   = gene_ranges$Path_max,
                    step  = 0.1
                  ),
                  tags$hr(),
                  actionButton(
                    "ct_run_single",
                    "Predict cell type",
                    icon = icon("microscope")
                  )
                )
              ),
              
              # 右：输出
              column(
                width = 8,
                box(
                  width = 12, status = "info", solidHeader = TRUE,
                  title = "Single-sample prediction",
                  
                  h4("Text summary"),
                  verbatimTextOutput("ct_single_text"),
                  tags$hr(),
                  
                  h4("Class probabilities"),
                  uiOutput("ct_single_probbar")
                  ,
                  
                  tags$hr(),
                  h4("Compare multiple single-sample trials"),
                  p(
                    "You can keep a small list of single-sample tries to compare ",
                    "different feature combinations."
                  ),
                  
                  div(
                    style = "display:flex; gap:8px; margin-bottom:8px;",
                    actionButton("ct_single_add",   "Add to compare", icon = icon("plus")),
                    actionButton("ct_single_clear", "Clear list",     icon = icon("trash")),
                    downloadButton("ct_single_dl",  "Download compare (.csv)")
                  ),
                  
                  DTOutput("ct_single_compare")
                )
              )
            )
          )
        )
      )
    ),
    
    tabPanel(
      title = tagList(icon("heartbeat"), "Disease prediction"),
      fluidPage(
        br(),
        h2("Disease status prediction (Tumor vs Healthy_Control)"),
        p(
          "This module deploys two Random Forest models for Disease_Status: ",
          "a basic model using genes and pathway scores only, and a ct model that ",
          "additionally incorporates cell type information (either true labels or ",
          "predicted by the cell type classifier)."
        ),
        br(),
        
        tabsetPanel(
          id = "dz_tabs",
          type = "tabs",
          
          ## ---------- A. Multiple (CSV) ----------
          tabPanel(
            title = tagList(icon("layer-group"), "Multiple (CSV)"),
            br(),
            fluidRow(
              # 左：上传 + 参数
              column(
                width = 4,
                box(
                  width = 12, status = "primary", solidHeader = TRUE,
                  title = "Step 1. Upload CSV or use example",
                  
                  fileInput(
                    "dis_multi_file",
                    "Upload CSV with gene features",
                    accept = ".csv"
                  ),
                  checkboxInput(
                    "dis_multi_use_example",
                    "Use example input from gene_input_example.csv",
                    value = FALSE
                  ),
                  tags$hr(),
                  h4("Tumor probability cutoff"),
                  numericInput(
                    "dis_multi_cutoff",
                    label = "Cutoff for p(Tumor):",
                    value = 0.5,
                    min   = 0.01,
                    max   = 0.99,
                    step  = 0.01
                  ),
                  helpText(
                    "Cells with p(Tumor) ≥ cutoff are classified as Tumor, ",
                    "otherwise Healthy_Control."
                  ),
                  tags$hr(),
                  actionButton(
                    "dis_multi_run",
                    "Run disease prediction",
                    icon = icon("play")
                  ),
                  br(), br(),
                  downloadButton(
                    "dis_multi_download",
                    "Download predictions (.csv)"
                  )
                )
              ),
              
              # 右：结果表 + 说明
              column(
                width = 8,
                box(
                  width = 12, status = "info", solidHeader = TRUE,
                  title = "Model comparison: basic vs ct",
                  p(
                    "Both models use the same gene features and engineered scores. ",
                    "The basic model only uses genes and scores. ",
                    "The ct model additionally incorporates cell type information ",
                    "(true labels in the dataset or predicted by the cell type ",
                    "classifier in user-uploaded data)."
                  ),
                  br(),
                  DTOutput("dis_multi_table")
                )
              )
            )
          ),
          
          ## ---------- B. Single sample ----------
          tabPanel(
            title = tagList(icon("vial"), "Single sample + Tip"),
            br(),
            fluidRow(
              # 左：输入区
              column(
                width = 4,
                box(
                  width = 12, status = "primary", solidHeader = TRUE,
                  title = "Input features",
                  
                  numericInput(
                    "dis_single_E", "Gene_E_Housekeeping",
                    value = round((gene_ranges$Gene_E_min + gene_ranges$Gene_E_max)/2, 1),
                    min   = gene_ranges$Gene_E_min,
                    max   = gene_ranges$Gene_E_max,
                    step  = 0.1
                  ),
                  numericInput(
                    "dis_single_A", "Gene_A_Oncogene",
                    value = round((gene_ranges$Gene_A_min + gene_ranges$Gene_A_max)/2, 1),
                    min   = gene_ranges$Gene_A_min,
                    max   = gene_ranges$Gene_A_max,
                    step  = 0.1
                  ),
                  numericInput(
                    "dis_single_B", "Gene_B_Immune",
                    value = round((gene_ranges$Gene_B_min + gene_ranges$Gene_B_max)/2, 1),
                    min   = gene_ranges$Gene_B_min,
                    max   = gene_ranges$Gene_B_max,
                    step  = 0.1
                  ),
                  numericInput(
                    "dis_single_C", "Gene_C_Stromal",
                    value = round((gene_ranges$Gene_C_min + gene_ranges$Gene_C_max)/2, 1),
                    min   = gene_ranges$Gene_C_min,
                    max   = gene_ranges$Gene_C_max,
                    step  = 0.1
                  ),
                  numericInput(
                    "dis_single_D", "Gene_D_Therapy",
                    value = round((gene_ranges$Gene_D_min + gene_ranges$Gene_D_max)/2, 1),
                    min   = gene_ranges$Gene_D_min,
                    max   = gene_ranges$Gene_D_max,
                    step  = 0.1
                  ),
                  numericInput(
                    "dis_single_path", "Pathway_Score_Inflam",
                    value = round((gene_ranges$Path_min + gene_ranges$Path_max)/2, 1),
                    min   = gene_ranges$Path_min,
                    max   = gene_ranges$Path_max,
                    step  = 0.1
                  ),
                  
                  tags$hr(),
                  h4("Disease model"),
                  radioButtons(
                    "dis_single_model",
                    label   = "Choose model:",
                    choices = c("Genes only (basic)" = "basic",
                                "Genes + Cell_Type (ct)" = "ct"),
                    selected = "ct"
                  ),
                  
                  numericInput(
                    "dis_single_thresh",
                    label = "Tumor probability cutoff:",
                    value = 0.5,
                    min   = 0.01,
                    max   = 0.99,
                    step  = 0.01
                  ),
                  
                  tags$hr(),
                  actionButton(
                    "dis_single_run",
                    "Run single-sample prediction",
                    icon = icon("heartbeat")
                  )
                )
              ),
              
              # 右：结果区
              column(
                width = 8,
                box(
                  width = 12, status = "info", solidHeader = TRUE,
                  title = "Model predictions",
                  h4("Text summary"),
                  verbatimTextOutput("dis_single_summary"),
                  tags$hr(),
                  
                  h4("Tumor probability"),
                  uiOutput("dis_single_probbar")
                ),
                br(),
                box(
                  width = 12, status = "warning", solidHeader = TRUE,
                  title = "Model-based educational tip (not real medical advice)",
                  verbatimTextOutput("dis_single_tip"),
                  tags$hr(),
                  # ===== 新增：AI 按钮 + 说明文字 =====
                  div(
                    style = "display:flex; align-items:center; gap:10px; margin-bottom:6px;",
                    actionButton(
                      "dis_ai_btn",
                      label = "Generate AI tip (optional)",
                      icon  = icon("robot"),
                      class = "btn btn-sm btn-info"
                    ),
                    span(
                      "(Gemini-based explanation in simple academic English)",
                      style = "font-size: 12px; color:#4b5563;"
                    )
                  ),
                  
                  # ===== 新增：AI 输出，不用 verbatim，避免滚动条 =====
                  textOutput("dis_ai_text")
                )
              )
            )
          )
        )
      )
    )
    
  ),
  
  # ===== Gene sandbox =====
  # ===== Gene sandbox / What-if experiment =====
  tabPanel(
    title = tagList(icon("sliders-h"), "Gene sandbox"),
    fluidPage(
      br(),
      h2("Gene sandbox / what-if experiment"),
      p(
        "This interactive module lets you pick a synthetic cell from the dataset, ",
        "tweak its gene / pathway features, and see how cell type and tumour ",
        "probabilities change. All results are based on synthetic BIO215 data ",
        "for teaching only."
      ),
      br(),
      
      fluidRow(
        # -------- 左：选择 baseline 细胞 + baseline 预测 --------
        column(
          width = 3,
          box(
            width = 12, status = "primary", solidHeader = TRUE,
            title = "Step 1. Choose a baseline cell",
            
            selectInput(
              "sandbox_cell",
              label   = "Cell_ID from dataset:",
              choices = data_df$Cell_ID,
              selected = data_df$Cell_ID[1]
            ),
            actionButton(
              "sandbox_load",
              "Load baseline cell",
              icon = icon("dna")
            ),
            tags$hr(),
            
            h4("Baseline cell (original)"),
            div(
              style = "max-height:210px; overflow-y:auto; overflow-x:auto; font-size:12px;",
              tableOutput("sandbox_before_info")
            ),
            
            tags$hr(),
            h4("Baseline predictions"),
            verbatimTextOutput("sandbox_before_pred")
          )
        ),
        
        # -------- 中：基因 / 通路滑块 --------
        column(
          width = 4,
          box(
            width = 12, status = "info", solidHeader = TRUE,
            title = "Step 2. Adjust gene / pathway features",
            
            helpText(
              "Start from the baseline values and move the sliders to perform a ",
              "simple what-if experiment. Then click 'Update prediction'."
            ),
            
            sliderInput(
              "sandbox_geneE", "Gene_E_Housekeeping",
              min   = gene_ranges$Gene_E_min,
              max   = gene_ranges$Gene_E_max,
              value = round((gene_ranges$Gene_E_min + gene_ranges$Gene_E_max) / 2, 1),
              step  = 0.1
            ),
            sliderInput(
              "sandbox_geneA", "Gene_A_Oncogene",
              min   = gene_ranges$Gene_A_min,
              max   = gene_ranges$Gene_A_max,
              value = round((gene_ranges$Gene_A_min + gene_ranges$Gene_A_max) / 2, 1),
              step  = 0.1
            ),
            sliderInput(
              "sandbox_geneB", "Gene_B_Immune",
              min   = gene_ranges$Gene_B_min,
              max   = gene_ranges$Gene_B_max,
              value = round((gene_ranges$Gene_B_min + gene_ranges$Gene_B_max) / 2, 1),
              step  = 0.1
            ),
            sliderInput(
              "sandbox_geneC", "Gene_C_Stromal",
              min   = gene_ranges$Gene_C_min,
              max   = gene_ranges$Gene_C_max,
              value = round((gene_ranges$Gene_C_min + gene_ranges$Gene_C_max) / 2, 1),
              step  = 0.1
            ),
            sliderInput(
              "sandbox_geneD", "Gene_D_Therapy",
              min   = gene_ranges$Gene_D_min,
              max   = gene_ranges$Gene_D_max,
              value = round((gene_ranges$Gene_D_min + gene_ranges$Gene_D_max) / 2, 1),
              step  = 0.1
            ),
            sliderInput(
              "sandbox_path", "Pathway_Score_Inflam",
              min   = gene_ranges$Path_min,
              max   = gene_ranges$Path_max,
              value = round((gene_ranges$Path_min + gene_ranges$Path_max) / 2, 1),
              step  = 0.1
            ),
            
            tags$hr(),
            actionButton(
              "sandbox_update",
              "Update prediction",
              icon = icon("play")
            )
          )
        ),
        
        # -------- 右：修改后的预测 + 图 + tip --------
        column(
          width = 5,
          
          box(
            width = 12, status = "success", solidHeader = TRUE,
            title = "Updated cell after slider changes",
            div(
              style = "max-height:210px; overflow-y:auto; overflow-x:auto; font-size:12px;",
              tableOutput("sandbox_after_info")
            ),
            tags$hr(),
            h4("Updated predictions"),
            verbatimTextOutput("sandbox_after_pred")
          ),
          br(),
          
          box(
            width = 12, status = "warning", solidHeader = TRUE,
            title = "Gene / score changes (before vs after)",
            plotOutput("sandbox_plot", height = 260)
          ),
          br(),
          
          box(
            width = 12, status = "danger", solidHeader = TRUE,
            title = "Educational what-if tip (not clinical advice)",
            
            # 原来的内置 tip
            verbatimTextOutput("sandbox_tip"),
            
            tags$hr(),
            
            # ===== 新增：AI 按钮 =====
            div(
              style = "display:flex; align-items:center; gap:10px; margin-bottom:6px;",
              actionButton(
                "sandbox_ai_btn",
                label = "Ask AI for extra explanation (optional)",
                icon  = icon("robot"),
                class = "btn btn-sm btn-danger"
              ),
              span(
                "(Gemini-based explanation for this what-if scenario)",
                style = "font-size: 12px; color:#4b5563;"
              )
            ),
            
            # ===== 新增：AI 输出 =====
            textOutput("sandbox_ai_text")
          )
          ,
          
          br(),
          
          box(
            width = 12, status = "info", solidHeader = TRUE,
            title = "Score-based explanation",
            p(
              "Click the button below to get a short explanation based on ",
              "changes in Tumor_score, Immune_score and Stromal_score."
            ),
            actionButton("sandbox_explain", "Explain changes"),
            br(), br(),
            verbatimTextOutput("sandbox_explain_text")
          )
        )
      ),
      
      hr(),
      em(
        "All data and models in this sandbox are synthetic BIO215 teaching materials. ",
        "They must not be used for personal health or any clinical decision making."
      )
    )
  ),
  
  # ===== INFO / Help =====
  # ===== INFO / Help =====
  tabPanel(
    title = tagList(icon("question-circle"), "INFO"),
    fluidPage(
      br(),
      h2("Help & documentation"),
      p(
        "This page explains the overall structure of the app, how to use each tab, ",
        "and important notes about the synthetic dataset and machine learning models."
      ),
      
      ## ---------- Row 1: Quick start + App structure ----------
      fluidRow(
        box(
          width = 6, solidHeader = TRUE, title = "Quick start",
          
          height = 385,  
          HTML("
        <ol>
          <li><b>Home</b>: Read the project overview, dataset background and workflow figure.</li>
          <li><b>Data</b>: Filter cells by <code>Disease_Status</code> or <code>Cell_Type</code>, inspect the interactive table and download subsets.</li>
          <li><b>Visualize</b>: Explore label distributions, bar charts, boxplots and histograms of genes or composite scores.</li>
          <li><b>Download</b>: Build custom subsets (rows + columns) and export CSV files, including an ML-ready feature table.</li>
          <li><b>Prediction → Cell type</b>: Upload a CSV (or use the example) to predict cell type for many cells, or try single-sample predictions with numeric inputs.</li>
          <li><b>Prediction → Disease</b>: Compare the basic and ct models on CSV input, or run a single-sample prediction with an educational risk tip.</li>
          <li><b>Gene sandbox</b>: Choose a real cell from the dataset, move sliders to perturb gene values and observe how predictions change (what-if simulation).</li>
        </ol>
        ")
        ),
        box(
          width = 6, solidHeader = TRUE, title = "App structure (tabs)",
          
          height = 385,  
          HTML("
        <ul>
          <li><b>DATA & VISUALIZATION</b>
            <ul>
              <li><b>Home</b>: Overall introduction, workflow diagram and dataset overview.</li>
              <li><b>Data</b>: \"Mini database\" view of the synthetic TME scRNA-seq dataset.</li>
              <li><b>Visualize</b>: Summary plots for labels, genes and scores.</li>
              <li><b>Download</b>: Subset builder and export.</li>
            </ul>
          </li>
          <li><b>PREDICTION</b>
            <ul>
              <li><b>Cell type prediction</b>: Random Forest model <code>rf_celltype</code>.</li>
              <li><b>Disease prediction</b>: Two Random Forest models (<code>rf_disease_basic</code> and <code>rf_disease_ct</code>).</li>
            </ul>
          </li>
          <li><b>GENE SANDBOX</b>: Interactive what-if module combining genes, scores and predictions.</li>
          <li><b>INFO</b>: The current help page.</li>
        </ul>
        <p>For the report you may include a screenshot of the navbar (top bar) to illustrate this structure.</p>
        ")
        )
      ),
      
      ## ---------- Row 2: Dataset background ----------
      fluidRow(
        box(
          width = 12, solidHeader = TRUE, title = "Dataset background",
          HTML("
        <h4>Context and inspiration</h4>
        <p>
        The app is built on a <b>synthetic single-cell RNA-seq</b> dataset that mimics 
        the <b>tumour microenvironment (TME)</b> of a solid tumour. Each row corresponds 
        to one simulated cell, with gene-expression features, a <code>Cell_Type</code> label 
        (Cancer, T_Cell, Fibroblast) and a <code>Disease_Status</code> label (Tumor vs Healthy_Control).
        </p>
        <p>
        The dataset is designed for teaching and benchmarking: it has clean signal, 
        balanced classes and no real patient information. It is suitable for practising 
        scRNA-seq style tasks such as cell type classification, feature exploration and 
        disease-status prediction.
        </p>
        
        <h4>Core analytical challenges</h4>
        <ul>
          <li><b>Cell type classification</b>: Predict Cancer / T_Cell / Fibroblast from gene and pathway features.</li>
          <li><b>Disease status prediction</b>: Use composite scores and cell-type information to estimate Tumor vs Healthy_Control.</li>
          <li><b>Biomarker exploration</b>: Inspect which genes or scores separate different cell types or disease labels.</li>
          <li><b>What-if simulations</b>: Use the Gene sandbox to see how changes in gene values might push model predictions.</li>
        </ul>
        
        <h4>Source and limitations</h4>
        <p>
        The data is <b>entirely synthetic</b> and was generated to resemble real scRNA-seq data 
        while avoiding privacy issues. It should <b>not</b> be used for any real biological or 
        clinical conclusions.
        </p>
        ")
        )
      ),
      
      ## ---------- Row 3: Data dictionary ----------
      fluidRow(
        box(
          width = 12, solidHeader = TRUE, title = "Data dictionary (main fields)",
          HTML("
        <div style='overflow:auto'>
        <table style='border-collapse:collapse;width:100%'>
          <thead>
            <tr>
              <th style='border:1px solid #e5e7eb;padding:6px'>Field</th>
              <th style='border:1px solid #e5e7eb;padding:6px'>Type</th>
              <th style='border:1px solid #e5e7eb;padding:6px'>Description</th>
            </tr>
          </thead>
          <tbody>
            <tr>
              <td style='border:1px solid #e5e7eb;padding:6px'>Cell_ID</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>character / integer</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>Unique identifier for each simulated cell.</td>
            </tr>
            <tr>
              <td style='border:1px solid #e5e7eb;padding:6px'>Cell_Type</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>factor</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>True cell type label (Cancer, T_Cell, Fibroblast).</td>
            </tr>
            <tr>
              <td style='border:1px solid #e5e7eb;padding:6px'>Disease_Status</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>factor</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>Tumor vs Healthy_Control label.</td>
            </tr>
            <tr>
              <td style='border:1px solid #e5e7eb;padding:6px'>Gene_A_Oncogene</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>numeric</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>Simulated oncogene-like expression; higher values tend to be more tumour-like.</td>
            </tr>
            <tr>
              <td style='border:1px solid #e5e7eb;padding:6px'>Gene_B_Immune</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>numeric</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>Immune-related marker, enriched in immune-like cells (e.g. T_Cell).</td>
            </tr>
            <tr>
              <td style='border:1px solid #e5e7eb;padding:6px'>Gene_C_Stromal</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>numeric</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>Stromal marker associated with fibroblast-like microenvironment.</td>
            </tr>
            <tr>
              <td style='border:1px solid #e5e7eb;padding:6px'>Gene_D_Therapy</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>numeric</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>Conceptual therapy-related marker / potential drug target.</td>
            </tr>
            <tr>
              <td style='border:1px solid #e5e7eb;padding:6px'>Gene_E_Housekeeping</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>numeric</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>Housekeeping-like gene capturing overall expression level.</td>
            </tr>
            <tr>
              <td style='border:1px solid #e5e7eb;padding:6px'>Pathway_Score_Inflam</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>numeric</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>Inflammatory pathway activity score.</td>
            </tr>
            <tr>
              <td style='border:1px solid #e5e7eb;padding:6px'>Tumor_score</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>numeric (derived)</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>Composite of <code>Gene_A_Oncogene</code> and <code>Gene_D_Therapy</code>.</td>
            </tr>
            <tr>
              <td style='border:1px solid #e5e7eb;padding:6px'>Immune_score</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>numeric (derived)</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>Composite of <code>Gene_B_Immune</code> and <code>Pathway_Score_Inflam</code>.</td>
            </tr>
            <tr>
              <td style='border:1px solid #e5e7eb;padding:6px'>Stromal_score</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>numeric (derived)</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>Composite score based on <code>Gene_C_Stromal</code>.</td>
            </tr>
            <tr>
              <td style='border:1px solid #e5e7eb;padding:6px'>UMAP_1</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>numeric</td>
              <td style='border:1px solid #e5e7eb;padding:6px'>First UMAP embedding dimension (for visualisation only).</td>
            </tr>
          </tbody>
        </table>
        </div>
        ")
        )
      ),
      
      ## ---------- Row 4: Input CSV format + prediction outputs ----------
      fluidRow(
        box(
          width = 6, solidHeader = TRUE, title = "Input file format (CSV for prediction tabs)",
          
          height = 370,  
          HTML("
        <p><b>Required feature columns</b> for both cell type and disease prediction:</p>
        <p>
          <code>Gene_E_Housekeeping, Gene_A_Oncogene, Gene_B_Immune, 
          Gene_C_Stromal, Gene_D_Therapy, Pathway_Score_Inflam</code>
        </p>
        <ul>
          <li>Column names must match <b>exactly</b> (case sensitive).</li>
          <li>All six feature columns must be numeric (decimal point as dot).</li>
          <li>Additional columns (e.g. <code>Cell_ID</code>, <code>Cell_Type</code>, <code>Disease_Status</code>) are allowed and will be ignored by the models.</li>
          <li>Use UTF-8 encoding and a comma <code>,</code> as the separator.</li>
        </ul>
        <p>It is strongly recommended to start from the provided <i>gene_input_example.csv</i> to avoid column name or ordering issues.</p>
        ")
        ),
        box(
          width = 6, solidHeader = TRUE, title = "Prediction outputs & thresholds",
          
          height = 370,  
          HTML("
        <ul>
          <li><b>Cell type prediction (multiple CSV)</b>
            <ul>
              <li><code>Pred_Cell_Type</code>: predicted Cancer / T_Cell / Fibroblast.</li>
              <li><code>Cancer, T_Cell, Fibroblast</code>: class probabilities (0–1).</li>
            </ul>
          </li>
          <li><b>Disease prediction (multiple CSV)</b>
            <ul>
              <li><code>Pred_Label_basic</code>, <code>Pred_Tumor_Prob_basic</code>, <code>Risk_basic</code>: results from the basic model.</li>
              <li><code>Pred_Label_ct</code>, <code>Pred_Tumor_Prob_ct</code>, <code>Risk_ct</code>: results from the ct model (uses cell-type information).</li>
            </ul>
          </li>
          <li><b>Single-sample prediction</b>
            <ul>
              <li>Shows engineered scores, predicted <code>Cell_Type</code> and <code>Pred_Disease</code>.</li>
              <li>Displays <code>p(Tumor)</code> with a progress bar and a qualitative risk label
                  (High / Intermediate / Low, simulated).</li>
            </ul>
          </li>
        </ul>
        <p>The Tumor cut-off used in disease prediction can be adjusted in the UI (e.g. 0.5 by default).</p>
        ")
        )
      ),
      
      ## ---------- Row 5: Gene sandbox explanation ----------
      fluidRow(
        box(
          width = 12, solidHeader = TRUE, title = "Gene sandbox (what-if module)",
          HTML("
        <p>
        The <b>Gene sandbox</b> module demonstrates a simple what-if analysis:
        </p>
        <ol>
          <li>Select a real cell from the dataset. The app loads its original gene and score values.</li>
          <li>Use sliders to simulate up- or down-regulation of specific genes or the inflammation score.</li>
          <li>Recalculate cell type and disease predictions for the modified profile.</li>
          <li>Compare \"before\" vs \"after\" values in tables, bar plots and risk tips.</li>
        </ol>
        <p>
        This module is intended to show how small changes in input features can influence 
        Random Forest predictions. It is <b>not</b> a mechanistic model of real biology.
        </p>
        ")
        )
      ),
      
      ## ---------- Row 6: Disclaimer / FAQ / browser ----------
      fluidRow(
        box(
          width = 6, solidHeader = TRUE, title = "Disclaimer & intended use",
          
          height = 225,  
          HTML("
        <p>
        All data and models in this project are <b>synthetic</b> and created for the BIO215 capstone project. 
        They are designed for teaching machine learning and data visualisation, not for clinical decision-making.
        </p>
        <ul>
          <li>Predictions should not be interpreted as diagnosis or treatment recommendations.</li>
          <li>Any \"risk\" labels are simulated and only illustrate how one might communicate model output.</li>
          <li>For real health concerns, always consult qualified medical professionals.</li>
        </ul>
        ")
        ),
        box(
          width = 6, solidHeader = TRUE, title = "FAQ / troubleshooting & browser",
          
          height = 225,  
          HTML("
        <ul>
          <li><b>Missing required columns</b>: Re-create the file starting from <i>gene_input_example.csv</i>.</li>
          <li><b>Non-numeric values</b>: Check for text strings, extra symbols or NA values in the gene columns.</li>
          <li><b>Prediction tab does nothing</b>: Make sure you clicked the <i>Run</i> button after selecting example or uploading a file.</li>
          <li><b>Best browser</b>: The app was tested mainly in Google Chrome; other modern browsers usually work as well.</li>
        </ul>
        ")
        )
      )
    )
  )
  
)

server <- function(input, output, session) {
  
  ## ---- Home: Dataset overview ----
  output$summary_ncells <- renderPrint({
    n_cells    <- nrow(data_df)
    n_cols     <- ncol(data_df)
    n_features <- length(feature_cols)
    cat(
      "Number of cells  :", n_cells, "\n",
      "Number of columns:", n_cols,  "\n",
      "Expression features (genes / pathway):", n_features, "\n",
      "Engineered scores:", 3, "\n"
    )
  })
  
  output$summary_ct <- renderTable({
    as.data.frame(table(data_df$Cell_Type)) |>
      dplyr::rename(Cell_Type = Var1, N = Freq)
  }, striped = TRUE, bordered = TRUE, digits = 0)
  
  output$summary_ds <- renderTable({
    as.data.frame(table(data_df$Disease_Status)) |>
      dplyr::rename(Disease_Status = Var1, N = Freq)
  }, striped = TRUE, bordered = TRUE, digits = 0)
  
  ## ---- Home: Workflow 图 ----
  output$workflow_img <- renderImage({
    # 你把图片放在 app 根目录或 www/ 下，文件名改成你想要的
    paths <- c(
      "www/TME_workflow.png",
      "www/TME_workflow.jpg",
      "TME_workflow.png",
      "TME_workflow.jpg"
    )
    p <- paths[file.exists(paths)][1]
    
    shiny::validate(
      shiny::need(length(p) == 1,
                  "Please put workflow image (TME_workflow.png or .jpg) into app/ or www/.")
    )
    
    list(
      src = p,
      contentType = if (grepl("\\.png$", p, ignore.case = TRUE)) "image/png" else "image/jpeg",
      alt = "Workflow of the Tumor Microenvironment Explorer",
      width = "100%"  # 自动铺满左侧列
    )
  }, deleteFile = FALSE)
  ## =========================================================
  ## Data tab: filters + interactive DT + download
  ## =========================================================
  
  # 1) 先做“外部筛选”：Disease_Status / Cell_Type
  data_filtered_base <- reactive({
    df <- data_df
    
    if (!is.null(input$data_filter_disease) &&
        input$data_filter_disease != "All") {
      df <- df[df$Disease_Status == input$data_filter_disease, , drop = FALSE]
    }
    if (!is.null(input$data_filter_celltype) &&
        input$data_filter_celltype != "All") {
      df <- df[df$Cell_Type == input$data_filter_celltype, , drop = FALSE]
    }
    df
  })
  # 2) DataTable 过滤之后的“当前数据” （≈ m6A 里的 data_filtered_now）
  data_filtered_now <- reactive({
    base <- data_filtered_base()
    idx  <- input$table_data_rows_all   # after DT internal filters
    if (is.null(idx)) base else base[idx, , drop = FALSE]
  })
  # 3) 再叠加 DataTable 自带的列过滤（header 下面那一行 filter）
  data_selected_now <- reactive({
    base <- data_filtered_now()
    sel  <- input$table_data_rows_selected
    if (is.null(sel) || !length(sel)) base[0, , drop = FALSE] else base[sel, , drop = FALSE]
  })
  
  
  # 4) 渲染 DataTable
  output$table_data <- DT::renderDT({
    df <- data_filtered_base()
    
    # 当前数据中实际存在的数值列
    num_cols <- intersect(num_cols_dt, names(df))
    
    DT::datatable(
      df,
      filter     = "top",              # 列头下面有过滤框
      rownames   = FALSE,
      selection  = "multiple",
      class      = "display nowrap compact",
      extensions = c("Buttons", "ColReorder", "FixedHeader", "Scroller"),
      options    = list(
        dom         = "Bfrtip",
        buttons     = list(
          list(extend = "colvis", text = "Columns"),
          list(extend = "copy",   text = "Copy"),
          list(extend = "csv",    text = "CSV"),
          list(extend = "excel",  text = "Excel"),
          list(extend = "print",  text = "Print")
        ),
        colReorder  = TRUE,
        fixedHeader = TRUE,
        deferRender = TRUE,
        scrollX     = TRUE,
        scrollY     = 400,
        scroller    = TRUE,
        pageLength  = 25,
        lengthMenu  = c(10, 25, 50, 100),
        autoWidth   = TRUE
      )
    ) |>
      # 数值列保留两位小数
      DT::formatRound(num_cols, digits = 2) |>
      
      # 疾病状态：Tumor 红色，Healthy_Control 绿色，加粗
      DT::formatStyle(
        "Disease_Status",
        color      = DT::styleEqual(
          c("Tumor",        "Healthy_Control"),
          c("#DC2626",      "#0D9488")        # 红 + 绿
        ),
        fontWeight = "bold"
      ) |>
      
      # 细胞类型：背景上色区分三类细胞
      DT::formatStyle(
        "Cell_Type",
        backgroundColor = DT::styleEqual(
          c("Cancer",   "T_Cell",   "Fibroblast"),
          c("#FCA5A5",  "#BFDBFE",  "#FEF3C7")   # 粉红/浅蓝/浅黄
        ),
        fontWeight = "bold"
      )
  }, server = TRUE)
  
  
  # 3) 三个 valueBox（行数 / 列数 / Tumor 比例）
  output$dt_rows <- renderValueBox({
    valueBox(
      scales::comma(nrow(data_filtered_now())),
      "Rows (after filters)",
      
      color = "light-blue"
    )
  })
  
  output$dt_cols <- renderValueBox({
    valueBox(
      ncol(data_df),
      "Columns",
      
      color = "teal"
    )
  })
  
  output$dt_posrate <- renderValueBox({
    d    <- data_filtered_now()
    rate <- if (nrow(d)) mean(d$Disease_Status == "Tumor") else 0
    valueBox(
      scales::percent(rate, accuracy = 0.1),
      "Tumor proportion (filtered)",
      
      color = "aqua"
    )
  })
  
  # 4) 清除选中行按钮（注意这里用的是你 UI 里定义的 id）
  observeEvent(input$data_clear_sel, {
    DT::selectRows(DT::dataTableProxy("table_data"), NULL)
  })
  
  # 5) 下载：Filtered rows / Selected rows
  output$data_dl_csv <- downloadHandler(
    filename = function() {
      scope <- if (identical(input$data_dl_scope, "selected")) "selected" else "filtered"
      paste0("TME_", scope, ".csv")
    },
    content = function(file) {
      out <- if (identical(input$data_dl_scope, "selected")) {
        x <- data_selected_now()
        if (!nrow(x)) data_filtered_now() else x
      } else {
        data_filtered_now()
      }
      write.csv(out, file, row.names = FALSE)
    }
  )
  ## =========================================================
  ## Visualize tab: donut pie + bar + box + histogram
  ## =========================================================
  
  # ---- 1) Donut pie chart ----
  output$pie_plot <- renderPlot({
    req(input$pie_col)
    var <- input$pie_col
    
    label_threshold <- 0.20   # 占比 < 20% 的扇区不显示文字
    
    df <- data_df %>%
      count(level = .data[[var]], name = "n") %>%
      arrange(desc(n)) %>%
      mutate(
        prop = n / sum(n),
        lab  = ifelse(
          prop >= label_threshold,
          paste0(level, "\n", scales::percent(prop, accuracy = 0.1)),
          NA_character_
        ),
        emph = ifelse(n == max(n), "max", "other")
      )
    
    lvls <- df$level
    
    # 不同变量用不同配色
    pal <- if (var == "Disease_Status") {
      setNames(
        c("#EF4444", "#0D9488"),        # 红 / 绿
        c("Tumor",   "Healthy_Control")
      )
    } else if (var == "Cell_Type") {
      setNames(
        c("#F97373", "#60A5FA", "#FACC15"),  # 粉 / 蓝 / 黄
        c("Cancer", "T_Cell", "Fibroblast")
      )
    } else {
      setNames(scales::hue_pal()(length(lvls)), lvls)
    }
    
    ggplot(df, aes(x = 2, y = n, fill = level, alpha = emph)) +
      geom_col(color = "white", width = 1) +
      coord_polar(theta = "y", clip = "off") +
      xlim(0.5, 2.5) +                      # 中间挖空 → 甜甜圈效果
      geom_text(
        aes(label = lab),
        position  = position_stack(vjust = 0.5),
        size      = 4,
        lineheight = 1.1
      ) +
      scale_fill_manual(values = pal, guide = guide_legend(title = var)) +
      scale_alpha_manual(values = c(max = 1, other = 0.6), guide = "none") +
      theme_void(base_size = 14) +
      theme(
        plot.title  = element_text(face = "bold", size = 16, hjust = 0.5),
        legend.title = element_text(size = 11),
        legend.text  = element_text(size = 10)
      ) +
      labs(title = paste("Distribution of", var))
  })
  
  # ---- 2) Bar chart: 按组的 mean/median ----
  output$bar_plot <- renderPlot({
    req(input$bar_gene, input$bar_group, input$bar_stat)
    var <- input$bar_gene
    grp <- input$bar_group
    
    stat_fun  <- if (identical(input$bar_stat, "median")) median else mean
    stat_name <- if (identical(input$bar_stat, "median")) "Median" else "Mean"
    
    df <- data_df %>%
      group_by(.data[[grp]]) %>%
      summarise(
        stat_value = stat_fun(.data[[var]], na.rm = TRUE),
        .groups    = "drop"
      )
    
    ggplot(df, aes(x = .data[[grp]], y = stat_value, fill = .data[[grp]])) +
      geom_col(color = "white") +
      scale_fill_brewer(palette = "Set2", guide = "none") +
      theme_minimal(base_size = 13) +
      labs(
        title = paste(stat_name, var, "by", grp),
        x     = grp,
        y     = paste(stat_name, var)
      )
  })
  
  
  
  # ---- 3) Boxplot: 分组分布 ----
  output$box_plot <- renderPlot({
    req(input$box_gene, input$box_group)
    var <- input$box_gene
    grp <- input$box_group
    
    ggplot(data_df,
           aes(x = .data[[grp]], y = .data[[var]], fill = .data[[grp]])) +
      geom_boxplot(alpha = 0.8) +
      theme_minimal(base_size = 13) +
      theme(legend.position = "none") +
      labs(
        title = paste(var, "by", grp),
        x     = grp,
        y     = var
      )
  })
  
  # ---- 4) Histogram: 下面一整块 ----
  hist_filtered <- reactive({
    d <- data_df
    if (!is.null(input$hist_disease) && input$hist_disease != "All") {
      d <- d[d$Disease_Status == input$hist_disease, , drop = FALSE]
    }
    if (!is.null(input$hist_celltype) && input$hist_celltype != "All") {
      d <- d[d$Cell_Type == input$hist_celltype, , drop = FALSE]
    }
    d
  })
  
  output$hist_plot <- renderPlot({
    req(input$hist_gene, input$hist_bins)
    d <- hist_filtered()
    shiny::validate(shiny::need(nrow(d) > 0, "No data for current filters."))
    
    var  <- input$hist_gene
    bins <- input$hist_bins
    
    x   <- d[[var]]
    mu  <- mean(x, na.rm = TRUE)
    md  <- median(x, na.rm = TRUE)
    q05 <- quantile(x, 0.05, na.rm = TRUE)
    q95 <- quantile(x, 0.95, na.rm = TRUE)
    
    ggplot(d, aes(x = .data[[var]])) +
      # 5–95% 范围阴影
      annotate(
        "rect",
        xmin = q05, xmax = q95,
        ymin = 0, ymax = Inf,
        fill = "#DBEAFE", alpha = 0.3
      ) +
      geom_histogram(
        bins  = bins,
        fill  = "#60A5FA",
        color = "white",
        linewidth = 0.3,
        alpha = 0.9
      ) +
      geom_vline(xintercept = mu, color = "#1D4ED8", linewidth = 0.8) +
      geom_vline(
        xintercept = md,
        color = "#EF4444",
        linewidth = 0.8,
        linetype = "dashed"
      ) +
      annotate(
        "text",
        x = mu, y = Inf, vjust = -0.6,
        label = sprintf("mean = %.2f", mu),
        color = "#1D4ED8",
        size  = 3.2,
        fontface = "bold"
      ) +
      annotate(
        "text",
        x = md, y = Inf, vjust = -2.0,
        label = sprintf("median = %.2f", md),
        color = "#EF4444",
        size  = 3.2,
        fontface = "bold"
      ) +
      labs(
        x = var,
        y = "Count"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.minor  = element_blank(),
        panel.grid.major.x = element_blank(),
        axis.title.x      = element_text(face = "bold"),
        axis.title.y      = element_text(face = "bold")
      )
  })
  
  ## =========================================================
  ## Download tab: filters + column selection + CSV export
  ## =========================================================
  
  # 1) 先按 Disease_Status / Cell_Type 做行筛选
  down_filtered <- reactive({
    df <- data_df
    
    # Disease_Status 多选（如果全选，就不筛）
    if (!is.null(input$down_disease) &&
        length(input$down_disease) > 0 &&
        length(input$down_disease) < length(levels(data_df$Disease_Status))) {
      df <- df[df$Disease_Status %in% input$down_disease, , drop = FALSE]
    }
    
    # Cell_Type 多选
    if (!is.null(input$down_celltype) &&
        length(input$down_celltype) > 0 &&
        length(input$down_celltype) < length(levels(data_df$Cell_Type))) {
      df <- df[df$Cell_Type %in% input$down_celltype, , drop = FALSE]
    }
    
    df
  })
  
  # 2) 再按选中的列来做列子集
  down_subset <- reactive({
    df   <- down_filtered()
    cols <- input$down_cols
    
    # 若用户把所有列都取消选中，就返回空表，避免报错
    if (is.null(cols) || !length(cols)) {
      return(df[0, , drop = FALSE])
    }
    
    cols <- intersect(cols, names(df))
    df[, cols, drop = FALSE]
  })
  
  # 3) 两个 valueBox：当前子集的行数 / 列数
  output$down_nrows <- renderValueBox({
    valueBox(
      value    = scales::comma(nrow(down_subset())),
      subtitle = "Rows in subset",
      
      color    = "light-blue"
    )
  })
  
  output$down_ncols <- renderValueBox({
    # 用列数，而不是 checkbox 长度（防止全部取消时 NA）
    df   <- down_subset()
    ncol_now <- ncol(df)
    
    valueBox(
      value    = ncol_now,
      subtitle = "Columns selected",
      
      color    = "teal"
    )
  })
  
  output$down_tumor <- renderValueBox({
    df <- down_subset()
    rate <- if (nrow(df) == 0 || !"Disease_Status" %in% names(df)) 0 else
      mean(df$Disease_Status == "Tumor")
    
    valueBox(
      value    = scales::percent(rate, accuracy = 0.1),
      subtitle = "Tumor proportion",
      
      color    = "aqua"
    )
  })
  down_subset_ml <- reactive({
    df <- down_filtered()  # 行：只用筛选，不受 down_cols 影响
    
    keep_cols <- c(
      feature_cols,                # 六个原始特征
      "Tumor_score",
      "Immune_score",
      "Stromal_score",
      "Cell_Type",
      "Disease_Status"
    )
    keep_cols <- intersect(keep_cols, names(df))
    
    df[, keep_cols, drop = FALSE]
  })
  
  output$down_btn_ml <- downloadHandler(
    filename = function() {
      paste0("TME_featurescores_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(down_subset_ml(), file, row.names = FALSE)
    }
  )
  
  # 4) 预览表：只显示前 50 行，避免太卡
  output$down_table <- DT::renderDT({
    df <- down_subset()
    DT::datatable(
      head(df, 50),
      rownames = FALSE,
      options  = list(
        pageLength = 10,
        lengthMenu = c(10, 25, 50),
        scrollX    = TRUE
      )
    )
  })
  
  # 5) 下载按钮：导出完整子集（所有行）
  output$down_btn <- downloadHandler(
    filename = function() {
      paste0("TME_subset_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(down_subset(), file, row.names = FALSE)
    }
  )
  
  ## =========================================================
  ## Cell type prediction tab
  ## =========================================================
  
  # ---------- Multiple (CSV) ----------
  ct_multi_input <- eventReactive(input$ct_run_multi, {
    if (isTRUE(input$ct_use_example)) {
      # 直接用示例输入
      example_input
    } else {
      file <- input$ct_file
      
      if (is.null(file)) {
        # 不再在页面上塞一行文字，而是弹一个小 notification
        showNotification(
          "Please upload a CSV file or tick 'Use example input from gene_input_example.csv'.",
          type     = "error",
          duration = 5
        )
        return(NULL)
      }
      
      read.csv(file$datapath, stringsAsFactors = FALSE)
    }
  })
  
  
  # 跑 rf_celltype，返回原始特征 + scores + Pred_Cell_Type + 三个概率
  ct_multi_pred <- eventReactive(input$ct_run_multi, {
    df <- ct_multi_input()
    req(!is.null(df)) 
    predict_celltype_df(df)
  })
  
  # 结果表
  output$ct_multi_table <- DT::renderDT({
    res <- ct_multi_pred()
    req(res)
    
    # 列顺序：预测在最前，然后特征 + scores + 概率
    ct_cols <- c(
      "Pred_Cell_Type",
      feature_cols,
      "Tumor_score", "Immune_score", "Stromal_score",
      "Cancer", "T_Cell", "Fibroblast"
    )
    ct_cols <- intersect(ct_cols, colnames(res))
    res <- res[, ct_cols, drop = FALSE]
    
    dt <- DT::datatable(
      res,
      rownames = FALSE,
      class    = "display nowrap compact",
      options  = list(
        pageLength = 10,
        lengthMenu = c(5, 10, 20, 50),
        scrollX    = TRUE
      )
    )
    
    # 能找到这些列再格式化，防止 example 空表时报错
    num_to_round <- intersect(
      c("Gene_E_Housekeeping", "Gene_A_Oncogene", "Gene_B_Immune",
        "Gene_C_Stromal", "Gene_D_Therapy", "Pathway_Score_Inflam",
        "Tumor_score", "Immune_score", "Stromal_score",
        "Cancer", "T_Cell", "Fibroblast"),
      colnames(res)
    )
    if (length(num_to_round)) {
      dt <- dt |>
        DT::formatRound(num_to_round, digits = 3)
    }
    
    # ⭐ 数字列全部居中
    num_cols_center <- names(res)[sapply(res, is.numeric)]
    if (length(num_cols_center)) {
      dt <- dt |>
        DT::formatStyle(
          num_cols_center,
          `text-align` = "center"
        )
    }
    
    # Pred_Cell_Type 上色，配色跟 Data 页相同
    if ("Pred_Cell_Type" %in% colnames(res)) {
      dt <- dt |>
        DT::formatStyle(
          "Pred_Cell_Type",
          backgroundColor = DT::styleEqual(
            c("Cancer",   "T_Cell",   "Fibroblast"),
            c("#FCA5A5",  "#BFDBFE",  "#FEF3C7")  # 粉/蓝/黄
          ),
          fontWeight = "bold",
          color      = "#000000"
        )
    }
    
    dt
  })
  
  
  # 下载预测结果
  output$ct_download_multi <- downloadHandler(
    filename = function() {
      paste0("CellType_predictions_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(ct_multi_pred(), file, row.names = FALSE)
    }
  )
  
  # ---------- Single sample ----------
  
  ct_single_pred <- eventReactive(input$ct_run_single, {
    df1 <- data.frame(
      Gene_E_Housekeeping = input$ct_single_E,
      Gene_A_Oncogene     = input$ct_single_A,
      Gene_B_Immune       = input$ct_single_B,
      Gene_C_Stromal      = input$ct_single_C,
      Gene_D_Therapy      = input$ct_single_D,
      Pathway_Score_Inflam = input$ct_single_path
    )
    predict_celltype_df(df1)
  })
  
  # 文本 summary
  output$ct_single_text <- renderPrint({
    res <- ct_single_pred()
    req(nrow(res) == 1)
    
    ct   <- as.character(res$Pred_Cell_Type[1])
    p_ca <- round(res$Cancer[1],     3)
    p_tc <- round(res$T_Cell[1],     3)
    p_fb <- round(res$Fibroblast[1], 3)
    
    cat("Predicted cell type:\n")
    cat("  ", ct, "\n\n", sep = "")
    
    cat("Class probabilities:\n")
    cat("  Cancer    :", p_ca, "\n")
    cat("  T_Cell    :", p_tc, "\n")
    cat("  Fibroblast:", p_fb, "\n\n")
    
    cat("Note: The model was trained on a fully synthetic TME dataset ",
        "and achieves near-perfect performance on that data. ",
        "All predictions here are for teaching purposes only.\n", sep = "")
  })
  
  output$ct_single_probbar <- renderUI({
    res <- ct_single_pred()
    req(nrow(res) == 1)
    
    probs <- c(
      Cancer     = res$Cancer[1],
      T_Cell     = res$T_Cell[1],
      Fibroblast = res$Fibroblast[1]
    )
    pred <- as.character(res$Pred_Cell_Type[1])
    
    # 为了排版更好看，少于 5% 的条也给个最小长度
    make_bar <- function(label, p, is_pred) {
      pct  <- round(p * 100)
      width_pct <- max(pct, 5)
      
      bar_class <- if (is_pred) {
        "progress-bar progress-bar-success"
      } else {
        "progress-bar"
      }
      
      tags$div(
        class = "ct-prob-row",
        tags$div(
          class = "ct-prob-label",
          sprintf("%s (%.1f%%)", label, pct)
        ),
        tags$div(
          class = "progress",
          tags$div(
            class = bar_class,
            role  = "progressbar",
            `aria-valuenow` = pct,
            `aria-valuemin` = "0",
            `aria-valuemax` = "100",
            style = paste0("width:", width_pct, "%;"),
            paste0(pct, "%")
          )
        )
      )
    }
    
    tagList(
      make_bar("Cancer",     probs["Cancer"],     pred == "Cancer"),
      make_bar("T_Cell",     probs["T_Cell"],     pred == "T_Cell"),
      make_bar("Fibroblast", probs["Fibroblast"], pred == "Fibroblast")
    )
  })
  
  ## ---------- Single sample: compare list ----------
  
  # 用来存多次尝试的结果
  ct_single_store <- reactiveVal(data.frame())
  
  # 点击 "Add to compare"：根据当前输入重新算一遍，然后追加一行
  observeEvent(input$ct_single_add, {
    # 按当前输入构造一行，直接调用你写好的 predict_celltype_df()
    df1 <- data.frame(
      Gene_E_Housekeeping = input$ct_single_E,
      Gene_A_Oncogene     = input$ct_single_A,
      Gene_B_Immune       = input$ct_single_B,
      Gene_C_Stromal      = input$ct_single_C,
      Gene_D_Therapy      = input$ct_single_D,
      Pathway_Score_Inflam = input$ct_single_path
    )
    res <- predict_celltype_df(df1)
    
    row <- data.frame(
      Gene_E_Housekeeping = res$Gene_E_Housekeeping,
      Gene_A_Oncogene     = res$Gene_A_Oncogene,
      Gene_B_Immune       = res$Gene_B_Immune,
      Gene_C_Stromal      = res$Gene_C_Stromal,
      Gene_D_Therapy      = res$Gene_D_Therapy,
      Pathway_Score_Inflam = res$Pathway_Score_Inflam,
      Tumor_score         = res$Tumor_score,
      Immune_score        = res$Immune_score,
      Stromal_score       = res$Stromal_score,
      Pred_Cell_Type      = as.character(res$Pred_Cell_Type),
      Cancer              = as.numeric(res$Cancer),
      T_Cell              = as.numeric(res$T_Cell),
      Fibroblast          = as.numeric(res$Fibroblast),
      stringsAsFactors    = FALSE
    )
    
    cur <- ct_single_store()
    ct_single_store(rbind(cur, row))
  })
  
  # 清空列表
  observeEvent(input$ct_single_clear, {
    ct_single_store(data.frame())
  })
  
  # 渲染 compare 表
  output$ct_single_compare <- DT::renderDT({
    dat <- ct_single_store()
    
    if (is.null(dat) || nrow(dat) == 0) {
      dat <- data.frame(
        Message = "List is empty. Click 'Add to compare' after running a single-sample prediction."
      )
      return(DT::datatable(dat, rownames = FALSE, options = list(dom = "t")))
    }
    
    # 列顺序：特征 + score + 预测
    cols <- c(
      feature_cols,
      "Tumor_score", "Immune_score", "Stromal_score",
      "Pred_Cell_Type", "Cancer", "T_Cell", "Fibroblast"
    )
    cols <- intersect(cols, colnames(dat))
    dat  <- dat[, cols, drop = FALSE]
    
    dt <- DT::datatable(
      dat,
      rownames = FALSE,
      options  = list(
        pageLength = 6,
        scrollX    = TRUE
      )
    )
    
    num_to_round <- intersect(
      c("Gene_E_Housekeeping", "Gene_A_Oncogene", "Gene_B_Immune",
        "Gene_C_Stromal", "Gene_D_Therapy", "Pathway_Score_Inflam",
        "Tumor_score", "Immune_score", "Stromal_score",
        "Cancer", "T_Cell", "Fibroblast"),
      colnames(dat)
    )
    if (length(num_to_round)) {
      dt <- dt |>
        DT::formatRound(num_to_round, digits = 3)
    }
    # ⭐ 数字列全部居中
    num_cols_center <- names(dat)[sapply(dat, is.numeric)]
    if (length(num_cols_center)) {
      dt <- dt |>
        DT::formatStyle(
          num_cols_center,
          `text-align` = "center"
        )
    }
    
    if ("Pred_Cell_Type" %in% colnames(dat)) {
      dt <- dt |>
        DT::formatStyle(
          "Pred_Cell_Type",
          backgroundColor = DT::styleEqual(
            c("Cancer",   "T_Cell",   "Fibroblast"),
            c("#F97373",  "#60A5FA",  "#FACC15")
          ),
          fontWeight = "bold",
          color      = "#000000"
        )
    }
    
    dt
  })
  
  # 下载 compare 列表
  output$ct_single_dl <- downloadHandler(
    filename = function() {
      paste0("CellType_single_compare_", Sys.Date(), ".csv")
    },
    content = function(file) {
      dat <- ct_single_store()
      if (is.null(dat) || nrow(dat) == 0) {
        dat <- data.frame(
          Message = "No entries in compare list when downloading."
        )
      }
      write.csv(dat, file, row.names = FALSE)
    }
  )
  ## =========================================================
  ## Disease prediction tab
  ## =========================================================
  
  ## ---------- Multiple (CSV): basic vs ct 比较 ----------
  
  # 读取上传 / 示例 CSV
  dis_multi_input <- eventReactive(input$dis_multi_run, {
    if (isTRUE(input$dis_multi_use_example)) {
      example_input
    } else {
      file <- input$dis_multi_file
      if (is.null(file)) {
        showNotification(
          "Please upload a CSV file or tick 'Use example input from gene_input_example.csv'.",
          type     = "error",
          duration = 5
        )
        return(NULL)
      }
      read.csv(file$datapath, stringsAsFactors = FALSE)
    }
  })
  
  # 调用 basic 和 ct 两个模型，并把结果拼在一张表里
  dis_multi_pred <- eventReactive(input$dis_multi_run, {
    df <- dis_multi_input()
    req(!is.null(df))
    
    cutoff <- input$dis_multi_cutoff
    
    # 只用原始特征 + scores 做一份“干净的基表”，方便展示
    base <- prep_new_samples(df)
    base_show <- base[, c(feature_cols,
                          "Tumor_score", "Immune_score", "Stromal_score"),
                      drop = FALSE]
    
    # 两个疾病模型
    basic_res <- predict_disease_basic_df(df, cutoff = cutoff)
    ct_res    <- predict_disease_ct_df(df, cutoff = cutoff)
    
    out <- cbind(
      base_show,
      Pred_Tumor_Prob_basic = basic_res$Pred_Tumor_Prob_basic,
      Pred_Label_basic      = basic_res$Pred_Label_basic,
      Risk_basic            = vapply(
        basic_res$Pred_Tumor_Prob_basic,
        risk_label,
        character(1)
      ),
      Pred_Tumor_Prob_ct    = ct_res$Pred_Tumor_Prob_ct,
      Pred_Label_ct         = ct_res$Pred_Label_ct,
      Risk_ct               = vapply(
        ct_res$Pred_Tumor_Prob_ct,
        risk_label,
        character(1)
      )
    )
    
    rownames(out) <- NULL
    out
  })
  
  # 结果表
  output$dis_multi_table <- DT::renderDT({
    res <- dis_multi_pred()
    req(res)
    
    # 调整一下列顺序：先 basic/ct 预测，再特征
    front_cols <- c(
      "Pred_Label_basic", "Pred_Tumor_Prob_basic", "Risk_basic",
      "Pred_Label_ct",    "Pred_Tumor_Prob_ct",    "Risk_ct"
    )
    front_cols <- intersect(front_cols, colnames(res))
    other_cols <- setdiff(colnames(res), front_cols)
    res <- res[, c(front_cols, other_cols), drop = FALSE]
    
    dt <- DT::datatable(
      res,
      rownames = FALSE,
      class    = "display nowrap compact",
      options  = list(
        pageLength = 10,
        lengthMenu = c(5, 10, 20, 50),
        scrollX    = TRUE
      )
    )
    
    # 数值列保留三位小数并居中
    num_cols <- names(res)[sapply(res, is.numeric)]
    if (length(num_cols)) {
      dt <- dt |>
        DT::formatRound(num_cols, digits = 3) |>
        DT::formatStyle(num_cols, `text-align` = "center")
    }
    
    # 预测标签：Tumor 红、Healthy_Control 绿
    if ("Pred_Label_basic" %in% colnames(res)) {
      dt <- dt |>
        DT::formatStyle(
          "Pred_Label_basic",
          backgroundColor = DT::styleEqual(
            c("Tumor",        "Healthy_Control"),
            c("#FCA5A5",      "#A7F3D0")
          ),
          fontWeight = "bold"
        )
    }
    if ("Pred_Label_ct" %in% colnames(res)) {
      dt <- dt |>
        DT::formatStyle(
          "Pred_Label_ct",
          backgroundColor = DT::styleEqual(
            c("Tumor",        "Healthy_Control"),
            c("#F97373",      "#6EE7B7")
          ),
          fontWeight = "bold"
        )
    }
    
    dt
  })
  
  # 下载预测结果
  output$dis_multi_download <- downloadHandler(
    filename = function() {
      paste0("Disease_predictions_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(dis_multi_pred(), file, row.names = FALSE)
    }
  )
  
  ## ---------- Single sample + tip ----------
  
  dis_single_result <- eventReactive(input$dis_single_run, {
    prediction_single(
      gene_E        = input$dis_single_E,
      gene_A        = input$dis_single_A,
      gene_B        = input$dis_single_B,
      gene_C        = input$dis_single_C,
      gene_D        = input$dis_single_D,
      path_inflam   = input$dis_single_path,
      dz_model      = input$dis_single_model,
      positive_threshold = input$dis_single_thresh
    )
  })
  
  # 文本 summary：特征、scores、预测标签、概率、模型
  output$dis_single_summary <- renderText({
    df <- dis_single_result()
    req(nrow(df) == 1)
    
    paste0(
      "Input features:\n",
      "  Gene_E_Housekeeping = ", round(df$Gene_E_Housekeeping, 3), "\n",
      "  Gene_A_Oncogene     = ", round(df$Gene_A_Oncogene, 3), "\n",
      "  Gene_B_Immune       = ", round(df$Gene_B_Immune, 3), "\n",
      "  Gene_C_Stromal      = ", round(df$Gene_C_Stromal, 3), "\n",
      "  Gene_D_Therapy      = ", round(df$Gene_D_Therapy, 3), "\n",
      "  Pathway_Score_Inflam= ", round(df$Pathway_Score_Inflam, 3), "\n\n",
      "Engineered scores:\n",
      "  Tumor_score   = ", round(df$Tumor_score, 3), "\n",
      "  Immune_score  = ", round(df$Immune_score, 3), "\n",
      "  Stromal_score = ", round(df$Stromal_score, 3), "\n\n",
      "Model predictions:\n",
      "  Predicted Cell_Type   = ", df$Pred_Cell_Type, "\n",
      "  Predicted Disease     = ", df$Pred_Disease, "\n",
      "  p(Tumor)              = ", round(df$Prob_Tumor, 3), "\n",
      "  Disease model used    = ", df$Model_Disease_Used, "\n"
    )
  })
  
  # 小进度条：显示 p(Tumor) 大小，并根据 cutoff 决定红/绿
  output$dis_single_probbar <- renderUI({
    df <- dis_single_result()
    req(nrow(df) == 1)
    
    p   <- as.numeric(df$Prob_Tumor[1])
    pct <- round(p * 100)
    width_pct <- max(pct, 5)  # 最小宽度，避免 1–2% 太细
    
    cutoff <- input$dis_single_thresh
    is_tumor <- p >= cutoff
    
    bar_class <- if (is_tumor) {
      "progress-bar progress-bar-danger"
    } else {
      "progress-bar progress-bar-success"
    }
    
    risk_txt <- risk_label(p)
    tagList(
      tags$p(
        sprintf("p(Tumor) = %.3f (cutoff = %.2f), risk level: %s",
                p, cutoff, risk_txt)
      ),
      tags$div(
        class = "progress",
        tags$div(
          class = bar_class,
          role  = "progressbar",
          `aria-valuenow` = pct,
          `aria-valuemin` = "0",
          `aria-valuemax` = "100",
          style = paste0("width:", width_pct, "%;"),
          paste0(pct, "%")
        )
      )
    )
  })
  
  # 教学 tip：用你写好的 generate_tip()
  output$dis_single_tip <- renderText({
    df <- dis_single_result()
    req(nrow(df) == 1)
    
    generate_tip(
      disease_label = as.character(df$Pred_Disease),
      immune_score  = df$Immune_score,
      p_tumor       = df$Prob_Tumor
    )
  })
  ## ---------- Disease 单样本：AI 额外 tip ----------
  # 先给一个空的默认值（可选）
  output$dis_ai_text <- renderText({ "" })
  
  observeEvent(input$dis_ai_btn, {
    df <- dis_single_result()
    req(nrow(df) == 1)   # 确保已经跑过单样本预测
    
    # 先显示“加载中”，避免用户以为没反应
    output$dis_ai_text <- renderText("Contacting AI (Gemini)...")
    
    # 整理一行供 make_summary_for_api 使用
    row <- df[1, ]
    
    summary_txt <- make_summary_for_api(row)
    
    prompt <- paste0(
      "You are helping a university student explain predictions from a synthetic tumour microenvironment (TME) dataset.\n",
      "The dataset is fully synthetic and for teaching only. It must never be used for medical decisions.\n\n",
      "Here is a short summary of the model output:\n",
      summary_txt, "\n\n",
      "Task: Write 2–3 short sentences in simple academic English that:\n",
      "1) Briefly describe what this pattern may mean for disease risk in this synthetic dataset.\n",
      "2) Clearly state that this is NOT medical advice and real decisions must be made by doctors.\n",
      "3) Optionally mention very general healthy habits (balanced diet, exercise, regular check-ups).\n",
      "Do NOT mention specific drugs, medical procedures, or dosages.\n",
      "Keep the answer within 80–120 words."
    )
    
    tip_text <- tryCatch(
      {
        call_gemini_tip(prompt)
      },
      error = function(e) {
        paste(
          "Gemini API error or no free credits available.\n",
          "Please refer to the built-in educational tip above.\n\n",
          "Technical message: ", e$message
        )
      }
    )
    
    # 最终把 AI 文本渲染出来（正常换行，没有滚动条）
    output$dis_ai_text <- renderText(tip_text)
  })
  
  
  
  ## =========================================================
  ## Gene sandbox tab: what-if experiment
  ## =========================================================
  
  # 存放 baseline 细胞
  sandbox_base <- reactiveVal(NULL)
  
  # Step 1: 根据选中的 Cell_ID 载入 baseline 细胞，并同步 slider
  observeEvent(input$sandbox_load, {
    row <- data_df[data_df$Cell_ID == input$sandbox_cell, , drop = FALSE]
    
    shiny::validate(
      shiny::need(nrow(row) == 1, "Selected cell not found in dataset.")
    )
    
    # 确保 scores 和当前规则一致（data_df 本身已算过，再跑一次也没关系）
    row <- add_scores(row)
    sandbox_base(row)
    
    # 用该细胞的真实值更新滑块
    updateSliderInput(session, "sandbox_geneE",
                      value = row$Gene_E_Housekeeping)
    updateSliderInput(session, "sandbox_geneA",
                      value = row$Gene_A_Oncogene)
    updateSliderInput(session, "sandbox_geneB",
                      value = row$Gene_B_Immune)
    updateSliderInput(session, "sandbox_geneC",
                      value = row$Gene_C_Stromal)
    updateSliderInput(session, "sandbox_geneD",
                      value = row$Gene_D_Therapy)
    updateSliderInput(session, "sandbox_path",
                      value = row$Pathway_Score_Inflam)
  })
  
  # Step 2: 根据 slider 修改基因，得到 after 版本
  sandbox_after_df <- eventReactive(input$sandbox_update, {
    base <- sandbox_base()
    shiny::validate(
      shiny::need(!is.null(base), "Please load a baseline cell first.")
    )
    
    modified <- base
    modified$Gene_E_Housekeeping  <- input$sandbox_geneE
    modified$Gene_A_Oncogene      <- input$sandbox_geneA
    modified$Gene_B_Immune        <- input$sandbox_geneB
    modified$Gene_C_Stromal       <- input$sandbox_geneC
    modified$Gene_D_Therapy       <- input$sandbox_geneD
    modified$Pathway_Score_Inflam <- input$sandbox_path
    
    modified <- add_scores(modified)
    modified
  }, ignoreNULL = FALSE)
  
  # ---- baseline 信息表 ----
  output$sandbox_before_info <- renderTable({
    base <- sandbox_base()
    shiny::validate(
      shiny::need(!is.null(base), "Load a baseline cell to see details.")
    )
    base[, c(
      "Cell_ID", "Cell_Type", "Disease_Status",
      "Gene_E_Housekeeping", "Gene_A_Oncogene", "Gene_B_Immune",
      "Gene_C_Stromal", "Gene_D_Therapy", "Pathway_Score_Inflam",
      "Tumor_score", "Immune_score", "Stromal_score"
    ), drop = FALSE]
  }, rownames = FALSE)
  
  # ---- after 信息表 ----
  output$sandbox_after_info <- renderTable({
    after <- sandbox_after_df()
    shiny::validate(
      shiny::need(!is.null(after), "Adjust genes and click 'Update prediction'.")
    )
    after[, c(
      "Cell_ID", "Cell_Type", "Disease_Status",
      "Gene_E_Housekeeping", "Gene_A_Oncogene", "Gene_B_Immune",
      "Gene_C_Stromal", "Gene_D_Therapy", "Pathway_Score_Inflam",
      "Tumor_score", "Immune_score", "Stromal_score"
    ), drop = FALSE]
  }, rownames = FALSE)
  
  # ---- baseline 预测 ----
  output$sandbox_before_pred <- renderPrint({
    base <- sandbox_base()
    shiny::validate(
      shiny::need(!is.null(base), "Load a baseline cell to see predictions.")
    )
    
    ct_res  <- predict_celltype_df(base)[1, ]
    dis_res <- predict_disease_basic_df(base)[1, ]
    
    cat("Baseline predictions (original cell)\n")
    cat("----------------------------------\n")
    cat("Cell type (rf_celltype): ",
        as.character(ct_res$Pred_Cell_Type), "\n", sep = "")
    cat("Tumor probability (basic disease model): ",
        round(dis_res$Pred_Tumor_Prob_basic, 3), "\n", sep = "")
    cat("Risk category: ",
        risk_label(dis_res$Pred_Tumor_Prob_basic), "\n")
  })
  
  # ---- after 预测 ----
  output$sandbox_after_pred <- renderPrint({
    after <- sandbox_after_df()
    shiny::validate(
      shiny::need(!is.null(after), "Adjust genes and click 'Update prediction' first.")
    )
    
    ct_res  <- predict_celltype_df(after)[1, ]
    dis_res <- predict_disease_basic_df(after)[1, ]
    
    cat("Updated predictions (after slider changes)\n")
    cat("-----------------------------------------\n")
    cat("Cell type (rf_celltype): ",
        as.character(ct_res$Pred_Cell_Type), "\n", sep = "")
    cat("Tumor probability (basic disease model): ",
        round(dis_res$Pred_Tumor_Prob_basic, 3), "\n", sep = "")
    cat("Risk category: ",
        risk_label(dis_res$Pred_Tumor_Prob_basic), "\n")
  })
  
  # ---- Before vs After 的柱状图 ----
  output$sandbox_plot <- renderPlot({
    base  <- sandbox_base()
    after <- sandbox_after_df()
    shiny::validate(
      shiny::need(!is.null(base),  "Load a baseline cell first."),
      shiny::need(!is.null(after), "Adjust genes and click 'Update prediction'.")
    )
    
    df_plot <- dplyr::bind_rows(
      data.frame(
        State = "Before",
        Gene  = gene_vars,
        Value = c(
          base$Gene_E_Housekeeping,
          base$Gene_A_Oncogene,
          base$Gene_B_Immune,
          base$Gene_C_Stromal,
          base$Gene_D_Therapy,
          base$Pathway_Score_Inflam,
          base$Tumor_score,
          base$Immune_score,
          base$Stromal_score
        )
      ),
      data.frame(
        State = "After",
        Gene  = gene_vars,
        Value = c(
          after$Gene_E_Housekeeping,
          after$Gene_A_Oncogene,
          after$Gene_B_Immune,
          after$Gene_C_Stromal,
          after$Gene_D_Therapy,
          after$Pathway_Score_Inflam,
          after$Tumor_score,
          after$Immune_score,
          after$Stromal_score
        )
      )
    )
    
    ggplot(df_plot, aes(x = Gene, y = Value, fill = State)) +
      geom_col(position = "dodge") +
      theme_minimal(base_size = 13) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
      labs(
        title = "Gene expression and scores: before vs after adjustment",
        x     = NULL,
        y     = "Expression / score"
      )
  })
  
  # ---- 教学 tip（基于 after）----
  output$sandbox_tip <- renderText({
    after <- sandbox_after_df()
    shiny::validate(
      shiny::need(!is.null(after), "Adjust genes and click 'Update prediction' first.")
    )
    
    dis_res <- predict_disease_basic_df(after)[1, ]
    p_tumor <- dis_res$Pred_Tumor_Prob_basic
    label   <- as.character(dis_res$Pred_Label_basic)
    
    generate_tip(
      disease_label = label,
      immune_score  = after$Immune_score,
      p_tumor       = p_tumor
    )
  })
  ## ---------- Sandbox：AI 额外解释 ----------
  output$sandbox_ai_text <- renderText({ "" })
  
  observeEvent(input$sandbox_ai_btn, {
    after <- sandbox_after_df()
    shiny::validate(
      shiny::need(!is.null(after),
                  "Please adjust sliders and click 'Update prediction' before asking AI.")
    )
    
    # 用 basic disease 模型 + cell type 模型算一下当前 after 状态
    ct_res  <- predict_celltype_df(after)[1, ]
    dis_res <- predict_disease_basic_df(after)[1, ]
    
    # 构造一行，字段名和 make_summary_for_api 需要的一致
    df_row <- data.frame(
      Pred_Cell_Type = as.character(ct_res$Pred_Cell_Type),
      Pred_Disease   = as.character(dis_res$Pred_Label_basic),
      Prob_Tumor     = as.numeric(dis_res$Pred_Tumor_Prob_basic),
      Tumor_score    = after$Tumor_score,
      Immune_score   = after$Immune_score,
      Stromal_score  = after$Stromal_score
    )
    
    output$sandbox_ai_text <- renderText("Contacting AI (Gemini)...")
    
    summary_txt <- make_summary_for_api(df_row[1, ])
    
    prompt <- paste0(
      "You are helping a university student interpret a what-if experiment in a synthetic tumour microenvironment (TME) sandbox.\n",
      "The dataset and models are fully synthetic and used only for teaching. They must never be used for real medical decisions.\n\n",
      "Here is a summary of the updated model output after slider changes:\n",
      summary_txt, "\n\n",
      "Task: Write 2–3 short sentences in simple academic English that:\n",
      "1) Explain qualitatively how this what-if change may have influenced the modelled tumour risk.\n",
      "2) Emphasise that this is not medical advice and that real decisions require doctors and real clinical tests.\n",
      "3) Optionally suggest very general healthy habits (balanced meals, regular exercise, no smoking, regular check-ups).\n",
      "Do NOT mention specific drugs, medical procedures, or dosages.\n",
      "Keep the answer within 80–120 words."
    )
    
    tip_text <- tryCatch(
      {
        call_gemini_tip(prompt)
      },
      error = function(e) {
        paste(
          "Gemini API error or no free credits available.\n",
          "Please refer to the built-in sandbox tip and score-based explanation above.\n\n",
          "Technical message: ", e$message
        )
      }
    )
    
    output$sandbox_ai_text <- renderText(tip_text)
  })
  
  
  # ---- 简单“解释”按钮：根据 score 变化给文字说明 ----
  sandbox_explain_res <- eventReactive(input$sandbox_explain, {
    base  <- sandbox_base()
    after <- sandbox_after_df()
    shiny::validate(
      shiny::need(!is.null(base),  "Load a baseline cell first."),
      shiny::need(!is.null(after), "Adjust genes and click 'Update prediction'.")
    )
    
    delta_tumor   <- after$Tumor_score   - base$Tumor_score
    delta_immune  <- after$Immune_score  - base$Immune_score
    delta_stromal <- after$Stromal_score - base$Stromal_score
    
    dis_base  <- predict_disease_basic_df(base)[1, ]
    dis_after <- predict_disease_basic_df(after)[1, ]
    
    list(
      delta_tumor   = delta_tumor,
      delta_immune  = delta_immune,
      delta_stromal = delta_stromal,
      prob_before   = dis_base$Pred_Tumor_Prob_basic,
      prob_after    = dis_after$Pred_Tumor_Prob_basic
    )
  })
  
  output$sandbox_explain_text <- renderPrint({
    x <- sandbox_explain_res()
    req(x)
    
    cat("Baseline p(Tumor):", round(x$prob_before, 3), "\n")
    cat("After adjustment p(Tumor):", round(x$prob_after, 3), "\n\n")
    
    if (x$delta_tumor < 0) {
      cat("- Tumor_score decreased. Lower oncogene / therapy-marker activity ",
          "likely reduces the predicted tumour probability.\n", sep = "")
    } else if (x$delta_tumor > 0) {
      cat("- Tumor_score increased. Higher oncogene / therapy-marker activity ",
          "may push the model towards the Tumor class.\n", sep = "")
    } else {
      cat("- Tumor_score did not change substantially.\n")
    }
    
    if (x$delta_immune > 0) {
      cat("- Immune_score increased. Stronger immune-related signalling ",
          "may partially counteract tumour-like patterns.\n", sep = "")
    } else if (x$delta_immune < 0) {
      cat("- Immune_score decreased. Weaker immune-related activity ",
          "may reduce protection against tumour-like cells.\n", sep = "")
    } else {
      cat("- Immune_score did not change substantially.\n")
    }
    
    if (x$delta_stromal != 0) {
      cat("- Stromal_score changed, indicating a shift in the stromal component ",
          "of the microenvironment.\n", sep = "")
    } else {
      cat("- Stromal_score did not change substantially.\n")
    }
    
    cat("\nThis explanation is heuristic and based on synthetic scores, ",
        "serving as an intuitive example of local model interpretation.\n",
        sep = "")
  })
}

shinyApp(ui, server)
