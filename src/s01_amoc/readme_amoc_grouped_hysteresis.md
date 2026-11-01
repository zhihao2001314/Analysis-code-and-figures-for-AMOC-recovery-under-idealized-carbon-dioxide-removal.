# n20：八模式AMOC–CO₂滞回曲线

## 输入、维度与指标

计算读取 `output/amoc_s01_n01/l_40/CMIP6_AMOC_9models_annual_series.csv` 的 `Forced AMOC strength (Sv)`，不读取旧订正异常列，也不读取原始NetCDF。源CSV实际维度为 `Model × Simulation year`，透视后为 `Simulation year × Model`。指标为各模式最近40°N纬度、500 m及以下Atlantic–Arctic翻转流函数最大值，不是固定500 m的流量，也不是局地深对流强度。实际纬度、成员及msftmz/msftyz来源记录于端点表；每模式一个已有成员。未做piControl漂移订正。

仅取合并第1–280年。升CO₂为1–140年，CDR为141–280年；C0=284.7 ppm，峰值1146.515147 ppm。源表晚于280年的记录不用于计算、平滑或补长。

## 分组及名称

沿用既有分层分类，不从滞回图外观重新分组：

1. 高初始层结组：ACCESS-ESM1-5、CanESM5、MIROC-ES2L。初始层结超过观测参考均值1.5倍，优先归入该组；该组三个模式的末期恢复程度不同。
2. AMOC恢复组：GFDL-ESM4、NorESM2-LM。通过初始层结筛选，且末期/初期AMOC≥0.8。
3. AMOC未恢复组：CESM2、CNRM-ESM2-1、UKESM1-0-LL。通过初始层结筛选，且末期/初期AMOC<0.8。

用户此次要求八模式，因此UKESM恢复其原始未恢复类别，仅用于本图；此前GFDL/Nor与CESM2/CNRM四模式大气复合不变。CAS不进入本图。

参考用户示意图，可在口头或正文中使用下列描述性名称，但不将其当作新统计分类，不在图中添加长解释：

| 模式 | 描述性名称 | 末期/初期（%） | 11年平滑面积A（Sv ppm） | 归一化面积A* |
| --- | --- | ---: | ---: | ---: |
| ACCESS-ESM1-5 | 端点近回归型滞回 | 100.5 | 724.3 | 0.0450 |
| CanESM5 | 端点近回归型滞回 | 102.0 | 1750.5 | 0.1554 |
| MIROC-ES2L | 持续偏弱型滞回 | 55.7 | 3057.9 | 0.2464 |
| GFDL-ESM4 | 端点近回归型滞回 | 98.0 | 2548.8 | 0.1607 |
| NorESM2-LM | 部分回归型滞回 | 88.1 | 4183.0 | 0.2168 |
| CESM2 | 持续偏弱型滞回 | 32.3 | 6782.1 | 0.3558 |
| CNRM-ESM2-1 | 持续偏弱型滞回 | 38.8 | 3958.2 | 0.2281 |
| UKESM1-0-LL | 持续偏弱型滞回 | 64.0 | 2704.4 | 0.1894 |

“近回归”不是统计等效于初始状态；“部分回归”不替代原≥0.8恢复分组。未设置large/small A的任意二分阈值。特别是Nor已达到原恢复阈值，但并非100%恢复。

## 计算与绘图方法

- 曲线：先取原AMOC年均强度，对完整1–280年连续序列做11年居中滑动平均，`min_periods=1`。边界窗口缩短，第1和280年的标记各平均6年。平滑穿过第140年转折点，但不使用第281年之后的数据。
- 横轴为线性CO₂浓度（ppm），纵轴为绝对AMOC强度（Sv），所有八幅图统一范围。蓝实线表示CO₂增加，朱红虚线表示CO₂降低；箭头均沿时间前进方向。空心圆为初始点，菱形为CO₂峰值，方形为末点，均是曲线上对应年份的平滑值。
- 为让两阶段曲线相连，降低分支从同一个第140年转折点开始绘制；这是共用已有转折点，不增加独立样本。
- 端点统计和分组使用未平滑1–10年及271–280年均值，和图上的单个平滑端点不是同一个估计量。CDR末期十年仍在降低CO₂，不能称为已在C0稳定十年。
- 同一CO₂精确配对：升高阶段第k年与合并第280−k年对应，k=1…140。第140年共用同一峰值。共同浓度域是287.547…1146.515147 ppm；降低末年280处于284.7 ppm，但没有升高第0年AMOC观测，因此保留该真实末点，不向低浓度外推升高分支，也不强制补闭合线。
- 面积：`A = ∫ |AMOC_down(C) − AMOC_up(C)| dC`，单位Sv ppm，只在共同CO₂域积分。将相邻浓度节点间分支差值过零的交点插入，再对绝对差值作梯形积分，避免交叉段面积高估。图中浅灰填色对应11年平滑分支间面积，不是置信区间。
- 归一化：`A* = A / [初始1–10年AMOC均值 × 共同CO₂浓度范围]`。它是平均相对路径分离程度，不是末期未恢复比例，也不照搬参考文献其他变量的归一化定义。
- 平滑敏感性：保存原年值以及5、11、21年平滑的面积。ACCESS的A*为0.0755、0.0600、0.0450、0.0321，对平滑较敏感；其他模式总体面积排序没有因此颠倒。不能用ACCESS的单一平滑面积宣称“无滞回”。
- 不新增显著性检验、闭合等效检验或组间置信区间。每组2或3个不同模式，单成员的有限速率轨迹不能单独证明平衡多稳态或永久不可逆。

## 关键观察

CESM2的面积和末期强度缺口均最大。GFDL与Nor即使达到原有恢复阈值，在同一中间CO₂浓度下，降低阶段AMOC仍普遍弱于升高阶段：端点恢复与路径滞回是两个问题。高初始层结组三个模式也不具有统一的恢复行为，尤其MIROC末期仍明显偏弱；该分组不等于“全部恢复”。

## 文件与运行

依次运行：

```powershell
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s01_amoc/n20_amoc_grouped_hysteresis.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'plot/s01_amoc/n20_amoc_grouped_hysteresis.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s01_amoc/check_amoc_grouped_hysteresis.py'
```

计算产品：`output/s01_amoc/n20_amoc_grouped_hysteresis/`。年值表、同CO₂配对表、填色节点、箭头节点、端点与面积指标、平滑敏感性、元数据和QC分开保存。绘图脚本只读取这些output结果，不重复平滑或积分。

图件：`figure/s01_amoc/n20_amoc_grouped_hysteresis/amoc_co2_grouped_hysteresis.{png,pdf,svg}`。PNG为4320×3210像素、300 dpi、白色不透明背景；PDF和SVG为矢量，SVG文字可编辑。图内没有长方法注释。当前为通用科研/PPT图，不宣称符合任何期刊的最终提交要求。

QC已通过：2240条模式年值、1120对同浓度节点、原始数据一致性、CO₂公式、端点均值、分组数量、计算绘图分离与无自编函数。独立跨零三角形公式复算面积，最大误差9.1×10⁻¹³ Sv ppm。最终PNG已人工查看：箭头沿时间方向，无标签裁切，八幅图使用相同坐标尺度。没有覆盖旧正式图件或结果。

## 参考与工具说明

- Kim, S.-K., et al. (2022). Widespread irreversible changes in surface temperature and precipitation in response to CO₂ forcing. Nature Climate Change, 12, 834–840. https://doi.org/10.1038/s41558-022-01452-z 。用户示意图为其Fig.1；本任务借用“路径面积”和“是否返回初态”分开讨论的框架，不直接套用其具体空间分类和显著性程序。2026-10-20核对。
- Scientific-visualization技能用于数据/变换/绘图分离、统一坐标尺度、颜色加线型冗余编码和实际图像检查，并未用于替代AMOC科学判断。软件方法参考：Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents. https://doi.org/10.48550/arXiv.2609.00065 。2026-10-20核对arXiv当前记录为2026-09-02修订版本。
