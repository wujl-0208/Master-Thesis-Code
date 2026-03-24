# 硕士毕业论文：基于多种抽样方案的临床评价指标统计性质研究与实证分析

本项目包含了 **中央民族大学 应用统计专业** 硕士学位论文的核心代码与实验过程。研究重点在于对比不同抽样设计（如 Case-Control, R-balanced）下，临床预测模型评价指标（AUC, IDI, NRI, SNB 等）的统计表现及其参数一致性。

## 👤 作者信息
* **作者**：吴俊良 (Wu Junliang)
* **学校**：中央民族大学 (Minzu University of China)
* **专业**：应用统计 (Applied Statistics)
* **项目地址**：[https://github.com/wujl-0208/Master-Thesis-Code](https://github.com/wujl-0208/Master-Thesis-Code)

---

## 📁 项目结构说明

本项目主要分为两个核心模块：**真实数据分析** 与 **数值模拟实验**。

### 1. 📈 真实数据分析 (Real Data Analysis)
基于 **NHANES (美国国家健康与营养调查)** 数据库，验证模型在临床场景下的评价效果。

* **数据预处理 (Python)**
    * `数据合并.ipynb`: 原始数据集的清洗、多表关联及缺失值处理。
    * `AUC增量.ipynb`: 计算模型引入新变量后的 AUC 提升情况。
* **统计分析与绘图 (R)**
    * `真实数据处理.R`: 核心逻辑回归建模及统计检验。
    * `calculate.R`: 辅助指标计算函数库。
* **数据与结果**
    * `nhanes_final_merged.csv`: 预处理后的最终分析数据集。
    * `校准曲线图_临床版.png` / `决策曲线分析_临床版.png`: 模型的临床有效性验证图表。

### 2. 🧪 数值模拟 (Numerical Simulation)
通过受控实验验证不同抽样方案下的指标性质。分为四个子文件夹：

* **`benchmark/`**: 基准对比模块，计算基础模型与全变量模型的指标差异。
* **`R-balanced/`**: 针对 R-balanced 抽样设计的模拟实验。
* **`casecontrol/`**: 针对个案对照抽样设计的模拟实验。
* **`beta/`**: 模拟结果汇总与可视化。
    * `箱线图.ipynb`: **核心绘图脚本**，用于读取各模拟文件夹数据，对比 $\beta$ 参数的理论值与模拟值并生成 `picture_beta_1.eps`。

**模拟核心脚本说明 (R-balanced & Case-control 内):**
* `Data_sampling.R`: 负责总体构建与样本抽取函数。
* `calculate.R`: 封装 AUC、IDI、NRI、SNB 等指标的计算逻辑。
* `main.R`: 主实验脚本。支持通过修改参数（如 `Prev`, `Rho`, `n`, `K` 等）调整模拟场景。

---

## 🛠 环境配置与依赖

为确保代码正常运行，请安装以下环境及扩展包：

### **Python (数据处理与可视化)**
* **版本**: 3.8+
* **核心库**:
    ```python
    import pandas as pd
    import numpy as np
    import seaborn as sns
    import matplotlib.pyplot as plt
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    ```

### **R (统计建模与复杂抽样)**
* **版本**: 4.0+
* **核心库**:
    ```R
    require(osDesign)   # 处理两阶段/复杂抽样设计
    require(mvtnorm)    # 多元正态分布处理
    require(plyr)       # 数据操纵
    require(MESS)       # 统计实用工具
    require(ggplot2)    # 绘图
    require(tidyr)      # 数据清洗
    ```

---

## 🚀 快速上手指南

1.  **数据分析**:
    * 首先运行 `数据合并.ipynb` 生成干净的数据集。
    * 在 R 环境中调用 `真实数据处理.R` 生成临床评价曲线。
2.  **执行模拟**:
    * 进入相应的抽样方案文件夹（如 `casecontrol`），运行 `main.R`。
    * 如需调整模拟设置（如修改患病率为 Rare Event），请在 `main.R` 的参数设定区修改 `Prev` 变量。
3.  **结果可视化**:
    * 确保各模拟文件夹生成了对应的 CSV 结果，运行 `beta/箱线图.ipynb` 生成论文所需的 $\beta$ 参数一致性对比图。

---

## 📄 许可说明
本项目代码仅供学术交流使用，引用请注明出处。
