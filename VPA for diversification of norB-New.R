# 加载必要的R包
library(readxl)
library(vegan)
library(VennDiagram)
library(ggplot2)

# 读取Excel文件
data <- read_excel("C:/Users/PC/Desktop/norB_pairwise_distances.xlsx")

# 查看数据结构
print(head(data))

# 确保数据是数值型
Geographic_data <- as.matrix(data[, c("distance", "Elevation","Moisture")])
Soil_properties_data <- as.matrix(data[, c("Totalnitrogen", "Totalcarbon","pH","Organicmatter",
                                         "Aluminum","Calcium","Copper","Iron","Potassium","Magnesium",
                                         "Manganese","Molybdenum","Sodium","Phosphorus","Sulfur",
                                         "Zinc")])

Climate_data <- as.matrix(data[, c("Precipitation", "Maxtemperature","Mintemperature","Windspeed")])
Land_use_data <- as.matrix(data[, c("Openwater", "Developed1","Developed2","Barren",
                                    "Forest","Shrubland","Grassland","Cropland",
                                    "Pasture","Wetland")])# 提取每个维度的数据
# 确保所有数据都是数值型
Geographic_data <- apply(Geographic_data, 2, as.numeric)
Soil_properties_data <- apply(Soil_properties_data, 2, as.numeric)
Climate_data <- apply(Climate_data, 2, as.numeric)
Land_use_data <- apply(Land_use_data, 2, as.numeric)

# 检查是否有非数值型数据
if (any(sapply(Geographic_data, is.character)) || any(is.na(Geographic_data)) ||
    any(sapply(Soil_properties_data, is.character)) || any(is.na(Soil_properties_data)) ||
    any(sapply(Climate_data, is.character)) || any(is.na(Climate_data)) ||
    any(sapply(Land_use_data, is.character)) || any(is.na(Land_use_data))) {
  stop("There are non-numeric or NA values in the data. Please check and handle them.")
}

# 处理缺失值
Geographic_data <- na.omit(Geographic_data)
Soil_properties_data <- na.omit(Soil_properties_data)
Climate_data <- na.omit(Climate_data)
Land_use_data <- na.omit(Land_use_data)

# 计算每个维度的综合贡献
Geographic_dissimilarity <- vegdist(Geographic_data, method = "euclidean")
Soil_properties_dissimilarity <- vegdist(Soil_properties_data, method = "euclidean")
Climate_dissimilarity <- vegdist(Climate_data, method = "euclidean")
Land_use_dissimilarity <- vegdist(Land_use_data, method = "euclidean")

# 假设的基因序列差异矩阵（示例数据）
norB_dissimilarity <- vegdist(as.matrix(data[, "norB"]), method = "euclidean")

# 进行VPA分析
vpa_result <- varpart(norB_dissimilarity, 
                      ~ Geographic_dissimilarity + Soil_properties_dissimilarity + 
                        Climate_dissimilarity + Land_use_dissimilarity)

# 绘制VPA结果
plot(vpa_result)

# 提取VPA结果中的贡献比例
Geographic_contribution <- vpa_result$var.part[1]
Soil_properties_contribution <- vpa_result$var.part[2]
Climate_contribution <- vpa_result$var.part[3]
Land_use_contribution <- vpa_result$var.part[4]

# 绘制Venn图
venn.plot <- draw.quad.venn(
  area1 = Geographic_contribution,
  area2 = Soil_properties_contribution,
  area3 = Climate_contribution,
  area4 = Land_use_contribution,
  category = c("Geographic", "Soil properties", "Climate", "Land use"),
  fill = c("#916FA999", "#D8A9C899", "#F9F0A899", "#A8CB7299"),
  alpha = 0.5
)

# 显示Venn图
grid.draw(venn.plot)

