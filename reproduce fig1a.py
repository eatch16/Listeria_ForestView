# 导入所需库
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# 设置数据文件和输出路径
blast_sum_path = r"D:\Reproduce_Fig_1a\\"
out_path = r"D:\Reproduce_Fig_1a\\"

# 读取功能性 ARG 数据
df_functional = pd.read_csv(
    blast_sum_path + "ARG_functional_prop.csv", header=0, names=["gene", "proportion"]
)
df_functional = df_functional[df_functional["proportion"] != 0]
df_functional["group"] = "Functional"

# 读取存在但截断的 ARG 数据
df_present = pd.read_csv(
    blast_sum_path + "ARG_present_prop.csv", header=0, names=["gene", "proportion"]
)
df_present = df_present[df_present["proportion"] != 0]
df_present["group"] = "Truncated"

# 合并两个数据表
df_ARG = pd.concat([df_functional, df_present], axis=0)

# 替换基因名称为更具可读性的标签
df_ARG = df_ARG.replace("lmo0919", "lin")
df_ARG = df_ARG.replace("lmo1695", "mprF")

# 设置绘图参数
plt.rcParams["figure.figsize"] = (4, 4)
sns.set(style="whitegrid")

# 绘制条形图
sns.barplot(
    data=df_ARG,
    x="gene",
    y="proportion",
    hue="group",
    palette=["#D98880", "#A9CCE3"],
    edgecolor="black",
)

# 图标题和坐标轴设置
plt.title("Prevalence of ARGs", size=15)
plt.legend(title="", bbox_to_anchor=(1, 1))
plt.xticks(fontsize=13, rotation=45, ha="right")
plt.xlabel("")
plt.ylabel("Fraction among Listeria genomes", fontsize=13)

# 添加柱状图数值注释
for p in plt.gca().patches:
    height = p.get_height()
    if not pd.isna(height) and height != 0:
        plt.annotate(
            format(height, '.6f'),
            (p.get_x() + p.get_width() / 2., height),
            ha='center', va='center',
            xytext=(0, 35),
            textcoords='offset points',
            rotation=90
        )

# 保存并显示图像
plt.savefig(out_path + "ARG_proportion.png", bbox_inches="tight", dpi=300)
plt.show()

