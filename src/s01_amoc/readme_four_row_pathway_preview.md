# 四行十四面板AMOC路径试排

当前布局更新：作者要求首行a与b各占两列，已同步两纬度版本。a范围x0.063–0.491，b范围0.539–0.967，宽均0.428；保持首行高度0.230及下三行全部位置。所有曲线和26个箭头数据逐项不变，下方12面板PNG逐像素相同。旧图及布局元数据归档archive/202610292039_equal_two_column_first_row/，18文件源目标与SHA256已核对。布局修订可用n72_equal_first_row_layout.py从该归档元数据复现；默认新生成n68本身也已使用b从第3列起的列位置。

最新修订：按作者要求，26.5°N与40°N两版的k都改成与c相同类型的平滑单调轮廓，保留蓝初态1、橙末态0.860708和共同峰态0.499460；低CO₂端仍开口。只映射示意轮廓端点，不增加实际模式数据平滑。旧两版34个正式产品已移动到archive/202610292033_smooth_open_k_both_latitudes/，源目标、文件数及SHA256核验。新旧PNG像素差仅位于k内部(425,4297)–(1445,4693)，其余13面板数据与显示不变。

作者确认共14个字母面板：a跨首行第1–2列，b跨第3–4列；以下三行各含一示意及三个模式。这里只输出试排，不替换原Fig1、不改分组注册表或下游显著性产品。

计算 `n68_four_row_pathway_preview.py`；绘图 `plot/s01_amoc/n69_four_row_pathway_preview.py`。数据 `output/s01_amoc/n68_four_row_pathway_preview/`，图件 `figure/s01_amoc/n69_four_row_pathway_preview/amoc_pathways_four_rows_preview.{png,pdf,svg}`。完整图片路径在绘图参数区列出并打印。

## 模式和示意

- 第二行c–f：Substantial recovery示意、GFDL、Nor、CAS。
- 第三行g–j：Non-recovery示意、CESM2、CNRM、UK。新图不再使用Incomplete recovery名称。
- 第四行k–n：Larger initial biases示意、ACCESS、CanESM5、MIROC。三个模式末期比例分别100.5%、102.0%、55.7%，因此初始偏差组不是统一的第三种恢复结果。

首行a的26.5°N异常、MMM/SD及CO₂输入逐字节复制原数据；b概念曲线和箭头坐标逐点保留。为容纳四列调整画布布局，首行科学内容不改。下方九个实际模式均取既有40°N/500m及以下AMOC、各自原始初始十年归一化、原11年居中平滑；没有新增平滑、漂移订正或强制闭合。

示意与真实模式严格分开：c保留原闭合基本恢复轮廓，并将概念峰值更新为三成员平均峰值参考；g保留原不恢复轮廓，更新三成员平均概念末态及峰值；k复用c两条平滑轮廓，以各自起点和共同峰值作仿射映射，概念末态仍为平均成员R=0.860708。k保持低CO₂端开口，包含的MIROC未恢复；灰色填充为分支分离，不是置信区间。概念末态取十年均值比，不冒充精确瞬时值。

输入CSV为model×simulation_year；输出path表为panel13×phase2×node（示意601、实际140/141），overview为simulation_year280×model9宽表。没有读取新NetCDF或假设维序。

13.2×12英寸，PNG5940×5400、450dpi，单页PDF950.4×864pt，SVG可编辑。14个主面板、26条实线路径及26个方向箭头全部核对，九模式路径逐点与n22吻合，所有源SHA256不变。已查看完整PNG及独立PDF栅格，无标签、图例或单位裁切。source/plot脚本均无自编函数或类，绘图仅读output。Scientific-visualization技能用于数据与示意区别、共用坐标、实际导出复审；没有声明期刊认证。新输出目录，无覆盖旧图。
