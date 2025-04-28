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
                     Treatment_0 = integer(),
                     Treatment_1 = integer(),
                     Treatment_2 = integer(),
                     Control_0 = integer(),
                     Control_1 = integer(),
                     Control_2 = integer(),
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
  
  # 计算各组基因型计数
  counts <- table(df$genotype, df$group)
  
  # 准备计数结果（确保所有基因型0,1,2都存在）
  treatment_counts <- rep(0, 3)  # 处理组（group=1）的0,1,2计数
  control_counts <- rep(0, 3)    # 对照组（group=0）的0,1,2计数
  
  for(gt in 0:2) {
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
  if(length(unique(df$genotype)) > 1) {
    model <- glm(group ~ genotype, data = df, family = binomial())
    p_value <- coef(summary(model))["genotype", "Pr(>|z|)"]
  } else {
    p_value <- NA
  }
  
  # 添加到结果
  results <- rbind(results, data.frame(SNP = snp,
                                      Chromosome = chr,
                                      Position = pos,
                                      p_value = p_value,
                                      Treatment_0 = treatment_counts[1],
                                      Treatment_1 = treatment_counts[2],
                                      Treatment_2 = treatment_counts[3],
                                      Control_0 = control_counts[1],
                                      Control_1 = control_counts[2],
                                      Control_2 = control_counts[3]))
}

# 输出结果
write.csv(results, "add.csv", row.names = FALSE)