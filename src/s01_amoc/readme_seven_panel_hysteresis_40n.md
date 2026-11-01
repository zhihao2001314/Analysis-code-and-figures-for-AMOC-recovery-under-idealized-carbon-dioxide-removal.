# 40°N七面板AMOC主图

## 默认版本

作者明确指定40°N为主要分析与稿件图版本，26.5°N保留与RAPID纬度对应的对照。后续未指定纬度的本图修改默认先更新40°N；若作者要求两版则同步更新。

- 计算与坐标准备：`src/s01_amoc/n79_seven_panel_hysteresis_40n.py`。
- 输出数据：`output/s01_amoc/n79_seven_panel_hysteresis_40n/`。
- 绘图：`plot/s01_amoc/n80_seven_panel_hysteresis_40n.py`。
- 主图PNG：`F:/sea/Antic sea ice constraint/figure/s01_amoc/n80_seven_panel_hysteresis_40n/amoc_hysteresis_seven_panels_40n_transparent.png`；同名PDF/SVG。
- 对照PNG：`F:/sea/Antic sea ice constraint/figure/s01_amoc/n78_seven_panel_hysteresis/amoc_hysteresis_seven_panels_26p5.png`。

## 纬度与计算一致性

图a与下方九模式往返轨迹全部为40°N最近原生纬度、500m以下AMOC诊断。原值来自已有40°N年度强迫强度，不使用piControl订正异常列。
坐标取既有n22的模型9×年份280的raw_sv和smooth_11yr_sv；每模型以未平滑初始1–10年均值归一化。图a由同一平滑强度减去同一初始均值得到。
两者均为严格在实验1–280年内的11年居中平均，首尾缩短窗口；图a旧版最后5个平滑年度使用稳定阶段的问题同步排除。前275年与既有40°N图a一致，所有实际滞回坐标保持原40°N值。
不读取原始NetCDF/MAT。输入CSV为model×simulation_year；绘图model_paths为2520行。初始、末期比值与实际滞回面积沿用n22/n73既有诊断。

## 布局与图例

共7面板，布局、组名、模型颜色与26.5°N版一致：a总览，b/c为ACCESS和CAS，d/e为GFDL、NorESM和CanESM，f/g为CESM、CNRM、UKESM和MIROC。只有d标注Hysteresis area(A)及指向阴影的箭头；f保留开放路径与阴影，不显示面积标注。
原独立面积概念图删除。示意b为严格零面积，d/f为已有40°N概念轮廓；只有示意b采用完全重合路径，实际模型值及实际小面积没有置零。
每个实际模型全程同色实线，18条箭头表明时间方向。

图a模式间标准差带alpha为0.18（82%透明）。CO₂强迫采用右轴数值坐标、在主轴最上层绘制，恢复灰色细虚线1.15pt，不加白色描边；图例背景alpha为0.15（85%透明），文字与模式颜色不变；图例层级高于强迫线；强迫数值与九模式数据均不变。

两版图a左下角括号均为各自纬度强迫实验初始1–10年未平滑平均AMOC强度（Sv，一位小数）。读取 `output/amoc_s01_n01/l_40/CMIP6_AMOC_9models_annual_series.csv` 中的 `Forced AMOC strength (Sv)`；26.5°N使用对应目录。按模型筛选Simulation year 1–10并等权平均；九模式均值为各模式初始均值的等权平均。
初始十年值与图a异常及实际轨迹归一化的参考一致，保留在 `initial_legend_values.csv`。本次图例更新不改变强迫曲线。

g图图例放在右上方，避免覆盖CESM2低值路径；图a及c/e/g的图例均自动检查无模式曲线交叠。

## 验证

两版均逐模型核对原始强迫强度、归一化路径、图a异常及初始十年图例值；图例修改前后所有曲线与方向箭头哈希相同。模型颜色匹配图a。
PNG4320×5670、450dpi，PDF单页、无Type3字体。两版PNG及独立PDF栅格均已目视检查。旧40°N十四面板图归档，当前两纬度各保留一套七面板图件。

本版直接复制作者认可的26.5°N绘图样式，仅替换40°N数据路径、纬度标识与输出文件名。文件名含transparent，以区别此前预览；当前图件目录仍仅保留这一套PNG/PDF/SVG。

最新标签修订：b为Recovery with small A，d为Recovery with large A，f仍为Non-recovery。两纬度同步。仅修改标签及图注，不更改曲线、面积或分组；b保留零面积极限示意。旧标签产品和图件归档：`archive/20261031233821_hysteresis_small_large_A_labels`。

2026-11-01示意修订：b更新为略弯曲、两条分支窄幅分离、两端闭合的小面积滞回，不再使用零面积直线；最大分支间距0.044（归一化AMOC），保持本纬度旧示意端点。箭头沿新分支调整；d/f、所有实际模式数据、图a强迫线及半透明样式均不变。修订脚本n82_small_area_hysteresis_schematic_revision.py。
