# 读取数据
data <- read.table("gwas.1e-6.genotype.tab", header = TRUE, na.strings = "NA", check.names = FALSE)

# 提取样本组别信息（D=处理组=1，S=对照组=0）
samples <- colnames(data)[3:ncol(data)]
group <- ifelse(grepl("^D", samples), 1, 0)  # 处理组=1，对照组=0

# 准备结果数据框
results <- data.frame(SNP = character(),
                     Chromosome = integer(),
                     Position = integer(),
                     p_value = numeric(),
                     Treatment_Homo = integer(),  # 纯合型(0+2)在处理组的样本数
                     Treatment_Het = integer(),   # 杂合型(1)在处理组的样本数
                     Control_Homo = integer(),    # 纯合型(0+2)在对照组的样本数
                     Control_Het = integer(),     # 杂合型(1)在对照组的样本数
                     stringsAsFactors = FALSE)

# 对每个SNP进行分析
for(i in 1:nrow(data)) {
  snp <- paste0(data[i,1], ":", data[i,2])
  chr <- data[i,1]
  pos <- data[i,2]
  
  # 提取基因型数据
  genotypes <- unlist(data[i, 3:ncol(data)])
  
  # 创建数据框用于分析
  df <- data.frame(genotype = genotypes, group = group)
  df <- df[!is.na(df$genotype), ]  # 移除缺失值
  
  # 将基因型转换为显性模型编码：0和2 -> 0（纯合型），1 -> 1（杂合型）
  df$dominant <- ifelse(df$genotype == 1, 1, 0)
  
  # 计算各组基因型计数（显性模型）
  counts <- table(df$dominant, df$group)
  
  # 准备计数结果（确保0和1都存在）
  treatment_counts <- rep(0, 2)  # 处理组（group=1）的纯合型(0)和杂合型(1)计数
  control_counts <- rep(0, 2)    # 对照组（group=0）的纯合型(0)和杂合型(1)计数
  
  for(gt in 0:1) {
    if(as.character(gt) %in% rownames(counts)) {
      if("1" %in% colnames(counts)) {  # 处理组（group=1）
        treatment_counts[gt+1] <- counts[as.character(gt), "1"]
      }
      if("0" %in% colnames(counts)) {  # 对照组（group=0）
        control_counts[gt+1] <- counts[as.character(gt), "0"]
      }
    }
  }
  
  # 进行GLM分析（仅当基因型有变异时）
  if(length(unique(df$dominant)) > 1) {
    model <- glm(group ~ dominant, data = df, family = binomial())
    p_value <- coef(summary(model))["dominant", "Pr(>|z|)"]
  } else {
    p_value <- NA
  }
  
  # 添加到结果
  results <- rbind(results, data.frame(SNP = snp,
                                      Chromosome = chr,
                                      Position = pos,
                                      p_value = p_value,
                                      Treatment_Homo = treatment_counts[1],
                                      Treatment_Het = treatment_counts[2],
                                      Control_Homo = control_counts[1],
                                      Control_Het = control_counts[2]))
}

# 输出结果
write.csv(results, "dom.csv", row.names = FALSE)