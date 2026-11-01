# 七模式直接mlotst与温盐诊断MLD的差值pattern

## 本轮目的

2026-11-01。把原n081两种方法的并列pattern直接相减，展示方法差异本身。
保留原七模式配对样本、有效格点和显示网格，不改变论文中采用的MLD来源。

计算代码：`src/s02_sia/s02_arctic_region/n095_mld_method_difference_patterns.py`。
绘图代码：`plot/S02_sia/s02_arctic_region/n096_mld_method_difference_patterns.py`。
计算结果：`output/s02_sia/s02_arctic_region/p24_mld_method_comparison/n095_mld_method_difference_patterns/`。
正式图件：`figure/s02_sia/s02_arctic_region/p24_mld_method_comparison/n096_mld_method_difference_patterns/`。

三张七面板图分别为：

- `march_mld_initial_derived_minus_mlotst.png/pdf/svg`：初始十年3月MLD方法差值。
- `march_mld_late_derived_minus_mlotst.png/pdf/svg`：末期十年3月MLD方法差值。
- `march_mld_change_derived_minus_mlotst.png/pdf/svg`：两种末减初响应的方法差值。

## 可用模式与当前原始资料复查

七个可作初末期配对的模式为ACCESS-ESM1-5、CanESM5、CESM2、CNRM-ESM2-1、GFDL-ESM4、NorESM2-LM、UKESM1-0-LL。
CAS仅有升高阶段mlotst，缺少CDR；MIROC两个强迫阶段均缺少本地mlotst。
当前逐文件重新读取NetCDF实际source_id、成员、试验、时间、维序和单位，确认七模式的两段各10个March全部存在。
这只确认本轮初末期比较资料完整，不声称各模式整个280年每月均完整。例如NorESM本地CDR文件始于源年160，但本轮需要的271–280年March完整。
完整覆盖审计见`current_mlotst_availability.csv`，实际文件维序见`current_mlotst_inventory.csv`。

| 模式 | 初始March源年 | 末期March源年 |
|---|---|---|
| ACCESS-ESM1-5 | 101–110 | 371–380 |
| CanESM5 | 1850–1859 | 2121–2130 |
| CESM2 | 1–10 | 131–140，CDR |
| CNRM-ESM2-1 | 1850–1859 | 2120–2129 |
| GFDL-ESM4 | 6–15 | 271–280 |
| NorESM2-LM | 1–10 | 271–280 |
| UKESM1-0-LL | 1850–1859 | 2120–2129 |

每个模式两方法严格使用相同月份；GFDL不把温盐6–15年与直接mlotst1–10年混比。

## 定义、维序与缺测

差值方向统一为温盐计算减直接mlotst，单位m。
初期与末期差值分别为两方法的绝对气候态相减；响应差值为：

`D_response=(derived_late-derived_initial)-(mlotst_late-mlotst_initial)`。

因此响应差值为正只表示温盐产品的末减初变化更正，不表示两方法都预测变深。
温盐方法沿用既有TEOS-10 sigma0、插值至10m的参考密度和0.03kg/m³增密阈值；无交点使用最深有效采样层。
完整配对、原生areacello加权和原始NetCDF维序记录见`readme_mld_method_comparison.md`。

n081输入实际维序：initial/late/change(model7,method2,latitude141,longitude341)，method顺序native、derived。
新差值数组维序：difference(model7,quantity3,latitude141,longitude341)，quantity顺序initial、late、change。
n079输入原生配对表为model×native_cell，34906行，记录真实经纬度和面积。
本轮直接相减已有0.25°共同显示场，不重做插值、平滑或温盐诊断；缺测保留NaN。
0.25°仅为显示间距。源网格的配对统计仍使用原生areacello权重，不从显示网格推断区域误差。
没有进行piControl漂移订正，没有把mlotst当作观测真值。

## 样式与核验

三图均为七模式两行布局，不加顶部总标题，不叠加不相关的初始MLD等值线。
当前显示统一采用Figure2的九锚点红蓝配色、BoundaryNorm及±1000m范围，50m分级，刻度−1000、−500、0、500、1000。
levels/ticks直接读取n093的grouped_mld_sic_fields.npz；颜色锚点与n094绘图脚本完全一致，不额外淡化颜色。
超范围差值原样保存，色标两端箭头表示延伸；原计算产品中的±600m参数是首次导出记录，当前显示以绘图manifest中的Figure2色标为准。
`display_difference_statistics.csv`保留原来的±600m范围诊断，当前±1000m范围的超限计数另外保存，不覆盖科学差值数组。
灰色代表陆地或缺少共同配对支持，不是零差值。此图不添加显著性检验。

计算端核验全部n081源文件SHA256一致，检查相减恒等式以及原生面积加权偏差与n079原统计一致。
`native_difference_statistics.csv`保存7模式×3时期量×3区域共63行的偏差、MAE、RMSE及极值。
绘图仅读取新output；每格式的完整绝对路径在绘图参数区构造并在保存后打印。
此前并排方法图和配对散点图没有被覆盖。

本次色标统一仅修改绘图；初始、末期和响应差值均逐值保留。旧三图及旧绘图代码已按原项目相对路径归档，归档路径与SHA256检查见本次archive目录中的archive_qc.json。
