# 七面板AMOC滞回图

> 作者已明确以40°N为论文主图，26.5°N保留对照。40°N当前主图见 `readme_seven_panel_hysteresis_40n.md`。两版图a图例均采用各自纬度强迫实验初始1–10年未平滑平均AMOC值（Sv，一位小数），与异常及归一化参考一致。

## 当前产物

- 计算与坐标整理：`src/s01_amoc/n77_seven_panel_hysteresis.py`。
- 输出：`output/s01_amoc/n77_seven_panel_hysteresis/`。
- 绘图：`plot/s01_amoc/n78_seven_panel_hysteresis.py`，只读output。
- PNG：`F:/sea/Antic sea ice constraint/figure/s01_amoc/n78_seven_panel_hysteresis/amoc_hysteresis_seven_panels_26p5.png`；同名PDF/SVG。

## 排版与成员

| 面板 | 内容 |
|---|---|
| a | 九模式26.5°N AMOC异常、等权均值±模式间总体标准差、CO₂强迫 |
| b/c | Recovery with small A示意；ACCESS-ESM1-5与CAS-ESM2-0叠加轨迹 |
| d/e | Recovery with large A示意；GFDL-ESM4、NorESM2-LM与CanESM5叠加轨迹 |
| f/g | Non-recovery示意；CESM2、CNRM-ESM2-1、UKESM1-0-LL与MIROC-ES2L叠加轨迹 |

原独立面积概念图删除。只有d把Hysteresis area(A)放在下方空白处，箭头指向阴影面积；f保留开放路径及阴影，移除面积文字与箭头。三幅示意用同一中性线色，靠方向箭头显示时间演变。实际模式各自沿用图a颜色，两阶段串成一条连续同色实线，无阶段配色图例。
图a居中且宽度适中，物理宽高比约2.2，没有横跨整页形成过扁的图框。

## 科学数据保留

所有实际轨迹来自已核对的n75全26.5°N结果。没有更换纬度、平滑、初始均值、CO₂浓度、末期恢复程度或实际滞回面积。
每模式两阶段都含峰值140年，合并时仅去掉同一峰值的重复副本，保留1–280年的全部原坐标。输出长表为model9×simulation_year280=2520行；各组2/3/4模型共九条连续曲线。
18条模式方向箭头的实际起止坐标沿用原值，颜色改为对应模式颜色。图a、CO₂强迫、端点诊断、实际面积文件逐字节复制。
no A是作者指定的概念标签：仅把b示意的两条路径设成完全重合，面积为零。没有把实际ACCESS/CAS的小面积置零，也没有强制实际曲线闭合。
另外两幅示意坐标保留。图a模式间标准差带alpha为0.18（82%透明），CO₂强迫采用右轴数值坐标、在主轴最上层绘制，恢复灰色细虚线1.15pt，不加白色描边；图例背景alpha为0.15（85%透明），文字与模式颜色不变；图例层级高于强迫线；不改任何数值或坐标。

## 验证与归档

九模式280个实际坐标与源序列逐点一致；模式颜色与旧图a逐项一致；绘图曲线哈希已核对。源文件哈希、面板数7、组名、18条方向箭头均已检查。
PNG4320×5670、450dpi；PDF单页、字体非Type3。PNG和独立PDF栅格均已目视检查。
旧26.5°N十四面板图已归档；40°N现已同步生成同样的七面板布局。

图a图例改为左下角单列，标题为Model (initial AMOC, Sv)，各模式与九模式等权平均后附一位小数的初始十年强度；数值保存为 `initial_legend_values.csv`。三组模式的图例中，g统一移至右上空白区，避免40°N下CESM2低值曲线遮挡。所有图例均检查过与实际模式曲线无交叠。

最新标签修订：b为Recovery with small A，d为Recovery with large A，f仍为Non-recovery。两纬度同步。仅修改标签及图注，不更改曲线、面积或分组；b保留零面积极限示意。旧标签产品和图件归档：`archive/20261031233821_hysteresis_small_large_A_labels`。

2026-11-01示意修订：b更新为略弯曲、两条分支窄幅分离、两端闭合的小面积滞回，不再使用零面积直线；最大分支间距0.044（归一化AMOC），保持本纬度旧示意端点。箭头沿新分支调整；d/f、所有实际模式数据、图a强迫线及半透明样式均不变。修订脚本n82_small_area_hysteresis_schematic_revision.py。
