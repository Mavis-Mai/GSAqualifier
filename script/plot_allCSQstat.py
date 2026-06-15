import matplotlib.pyplot as plt
import pandas as pd
import os
import sys
import matplotlib.patches as mpatches

def read_stat_files(file_list):
    """读取统计文件并转换为数据框"""
    data = []
    for file in file_list:
        with open(file, 'r') as f:
            lines = f.readlines()
            if len(lines) < 2:
                continue  # 跳过空文件

            stat_dict = {}
            header = lines[0].strip().split('\t')  # 提取标题
            for line in lines[1:]:
                parts = line.strip().split('\t')
                #if len(parts) == len(header):
                stat_dict[parts[0]] = dict(zip(header, parts))

            # 单文件内的总数计算，包含 Missing
            total = sum(
                int(stat_dict[label]["count"])
                for label in ["Good", "Moderate", "Poor"]
                if label in stat_dict
            ) + int(stat_dict.get("Missing", {}).get("count", 0))

#            print(stat_dict)
            # 基于该文件统计并存储数据
            data.append({
                "file": os.path.basename(file).replace("_stat.txt", ""),  # 提取文件名前缀作为标签
                "good_freq": float(stat_dict["Good"]["frequency(%)"]),
                "moderate_freq": float(stat_dict["Moderate"]["frequency(%)"]),
                "poor_freq": float(stat_dict["Poor"]["frequency(%)"]),
                "missing_pct": int(stat_dict.get("Missing", {}).get("count", 0)) / total * 100 if total > 0 else 0,  # 转换为百分比
            })
 #   print(data)
    return pd.DataFrame(data)


def plot_stacked_bar(data, output_prefix):
    """绘制横向堆叠柱状图"""
    # 准备数据
    y_labels = data["file"]
    good_freq = data["good_freq"]
    moderate_freq = data["moderate_freq"]
    poor_freq = data["poor_freq"]
    missing_pct = data["missing_pct"]

    # 创建图形
    fig, ax = plt.subplots(figsize=(6, 6))

    # 绘制横向堆叠柱状图
    ax.barh(y_labels, good_freq, label="Good", color="#EF888A")      # 红色
    ax.barh(y_labels, moderate_freq, left=good_freq, label="Moderate", color="#45465e")  # 灰色
    ax.barh(y_labels, poor_freq, left=good_freq + moderate_freq, label="Poor", color="#b8b8b8")  # 深灰蓝
    ax.barh(y_labels, missing_pct,
            left=good_freq + moderate_freq + poor_freq,
            label="Missing",
            color="black",
            hatch="//", edgecolor="gray")  # 斜纹填充+灰色

    for i, g in enumerate(good_freq):
        ax.text(g / 2, i, f"{g:.1f}%", va="center", ha="center", fontsize=8, color="white")

    # 添加图例
    good_patch = mpatches.Patch(color="#EF888A", label="Good")  # 自定义图例
    moderate_patch = mpatches.Patch(color="#45465e", label="Moderate")
    poor_patch = mpatches.Patch(color="#45465e", label="Poor")
    missing_patch = mpatches.Patch(color="black", label="Missing", hatch="//", edgecolor="gray")
    ax.legend(handles=[good_patch, moderate_patch, poor_patch, missing_patch], loc="upper left", title="Categories")

    # 设置轴标签和标题
    ax.set_xlabel("Percentage (%)")
    ax.set_title("Stacked Horizontal Bar Chart of Label Frequencies")

    # 设置X轴范围，确保Missing部分能够显示超出100%的区域
  #  print(max(missing_pct))
    ax.set_xlim(0, 100 + max(missing_pct))
    ax.axvline(x=100, color="gray", linestyle="--", linewidth=1)

    # 保存图像
    plt.tight_layout()
    plt.savefig(f"{output_prefix}.svg", format="svg", dpi=200)
    plt.savefig(f"{output_prefix}.png", format="png", dpi=200)
    print(f"Saved plots as {output_prefix}.svg and {output_prefix}.png")


def main(args_list=None):
    # 检查输入参数
    if len(sys.argv) < 2:
        print("Usage: python script.py file1.txt file2.txt ...")
        sys.exit(1)

    input_files = sys.argv[1:]  # 从命令行获取文件列表
    output_prefix = "stacked_horizontal_bar_chart"

    # 读取数据并生成堆叠柱状图
    data = read_stat_files(input_files)
    plot_stacked_bar(data, output_prefix)


if __name__ == "__main__":
    main()
