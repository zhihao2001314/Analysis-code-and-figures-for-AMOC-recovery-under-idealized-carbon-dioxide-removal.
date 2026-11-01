> 当前26.5°N图已更新为n77/n78七面板版；十四面板图已归档。见 `readme_seven_panel_hysteresis.md`。下文记录原十四面板的科学计算口径，源数据仍保留。

# 全部采用26.5°N的AMOC滞回图

## 修正内容

此前n75只换图a、下方沿用40°N，是错误的混合纬度版本，已整体归档。
当前n75重新读取26.5°N的九模式年度强迫原值，重算图a、全部九条滞回路径、时间箭头、面积和初末强度。没有复制40°N模式轨迹。
概念图b和三幅示意只复用抽象形状与布局，示意峰值参考改为26.5°N组平均；它们不是模拟数据。

## 输入及维度

- 原始诊断来源：`output/amoc_s01_n01/l_26.5/CMIP6_AMOC_9models_annual_series.csv`。
- 只用 `Forced AMOC strength (Sv)`，不使用piControl订正异常列。
- AMOC定义：各模式最接近26.5°N原生纬度、500m以下大西洋—北极翻转流函数最大值。实际纬度及成员保存在 `latitude_mapping.csv`。
- 输入长表为Model×Simulation year，每模式301个年度值；本次严格选1–280年。不重新读取NetCDF/MAT。
- 年份1–140为CO₂升高，140–280构造降低分支；峰值140为共用转折点。
- 11年居中平均仅在1–280年内计算，首尾缩短窗口。图a与下方路径用同一套平滑强度；末期稳定阶段未参与平滑。
- 初始和末期均值分别来自未平滑的1–10年、271–280年。路径除以各自初始均值；图a减去同一初始均值。
- 因旧图a平滑曾使用280年以后的样本，本次图a最后5个平滑年度同步修正；前275年与既有26.5°N概览逐点一致。

## 分组及结果

保留作者确认的三行2/3/4模式排列，第一组标签为Closed, small area。
三组实际无量纲面积（绝对分支差积分除以CO₂跨度）分别约为0.047、0.127–0.130、0.160–0.231。
这是描述性试排，不据此建立物理阈值，也不强制实际端点闭合。示意c仍是作者指定的小面积闭合直线往返图。

CESM2末期/初始十年均值约0.638，GFDL约1.141，CAS约0.898；26.5°N和40°N对照已保存为 `comparison_with_40n.csv`。
滞回面积在两分支共同浓度范围内计算，使用原CO₂节点和所有交叉点，对绝对差分段线性积分；面积原单位Sv×ppm也已保存。

## 产物与核对

- 计算脚本：`src/s01_amoc/n75_hysteresis_morphology_26p5.py`。
- 数据及QA：`output/s01_amoc/n75_hysteresis_morphology_26p5/`。
- 绘图脚本：`plot/s01_amoc/n76_hysteresis_morphology_26p5.py`，只读output。
- 实际PNG：`F:/sea/Antic sea ice constraint/figure/s01_amoc/n76_hysteresis_morphology_26p5/amoc_hysteresis_morphology_26p5.png`；同名PDF和SVG。
- 9模式×2分支共18条曲线均已与26.5°N年度原值逐点核对，且18条都与40°N对应路径不同。与图a间的初始强度归一化关系逐点一致。
- PNG5940×5400、450dpi；PDF为单页。PNG与PDF渲染均已目视检查。保持十四面板，无总标题。
