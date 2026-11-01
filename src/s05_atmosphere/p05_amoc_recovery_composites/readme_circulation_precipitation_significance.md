# n71：年平均流函数强度核验与降水/流函数显著性（2026-10-27）

最新显示更新：显著性图例移至图注，图内黑点及网纹完全保留，其他数据与标记不变。图注用语保存在绘图manifest的 `significance_caption` 中。继续收紧顶部空白，当前尺寸10.4×9.6英寸、4680×4320像素，PDF748.8×691.2pt；面板物理大小不变。上一版及QC归档到 `archive/202610272055_significance_legend_to_caption/`，源目标状态及哈希已核对。下文尺寸为历史记录，以本段为准。scientific-visualization技能只用于布局/导出检查，来源文献已记录于文末。

本日后续显示更新：依用户要求去掉主图总标题与时期副标题，图注说明背景；保持面板、字号和图例的物理大小，收去顶部0.6英寸标题区域，当前主图10.4×10.0英寸、4680×4500像素。下文4770像素/763.2pt为去标题前记录，当前PDF高度720pt。全部科学数据和显著性计算不变。旧版本及QC核验归档至 `archive/202610272047_remove_figure_title/`；项目AGENTS已记录后续图默认不加总标题的偏好，不批量修改无关旧图。

## 当前范围与文件

本次仅更新 n66 的 annual_circulation_global_precipitation：左下 e 的流函数响应组差和右侧 b/d/f 三幅降水场添加显著性。a/c 未添加检验，DJF 原图未改；global_precipitation_itcz.png 与归档旧版逐字节相同。原始资料、原有填色数值、气候态等值线、雨带格点位置及响应差值方向均未更改。沿用原强迫资料口径，本任务未引入 piControl 订正。

计算：`src/s05_atmosphere/p05_amoc_recovery_composites/n71_circulation_precipitation_significance.py`。
强度核验：同目录 `check_circulation_precipitation_strength.py`。
独立验收：同目录 `check_circulation_precipitation_revision.py`。
绘图：`plot/s05_atmosphere/p05_amoc_recovery_composites/n66_circulation_precipitation_rainband.py`，只读取 output。
正式计算产品位于 `output/s05_atmosphere/p05_amoc_recovery_composites/n71_circulation_precipitation_significance/`，正式图件仍在对应 figure 的 n66 目录。

## 强度核验及边界

本图是大气纬向平均质量流函数，不是海洋 AMOC 体积输送（Sv）：

`Psi = 2*pi*a*cos(lat)/g * integral_[10hPa,p] zonal_mean(v) dp`，存储及色标单位为 `10^9 kg s^-1`。气压积分使用 Pa，初始等值线北半球 Hadley 环流为正、南半球为负。正异常不等于两半球环流都增强。

独立重新读取四模式初期、末期各一个完整年，共96个月原始 va，在南北10°、300 hPa做16次独立复算，与原 n35 最大误差为 `2.84e-14`（10^9 kg s^-1）。NorESM2-LM 原始维度 `(time,lev,lat,lon)`，混合气压为 `a*p0+b*ps`；另三模式为 `(time,plev,lat,lon)`，压力单位 Pa。风速为 m s^-1。独立 Nor 逐柱 log-pressure 插值复现原层选择算法。

初始十年年平均、30°S–30°N、200–700 hPa极值（10^9 kg s^-1）：

| Model | Maximum | Minimum |
| --- | ---: | ---: |
| GFDL-ESM4 | 78.739 | -100.679 |
| NorESM2-LM | 73.637 | -77.900 |
| CESM2 | 77.768 | -112.006 |
| CNRM-ESM2-1 | 95.996 | -116.448 |

绝对极值约 `0.74–1.16e11 kg s^-1`，与下列研究的年平均 Hadley 环流约 `1e11 kg s^-1`量级一致；不能将冬季约 `2e11 kg s^-1`的极值直接当年平均目标。门槛只核对计算、单位、符号及宽泛量级，不是观测偏差检验，也不是模式技巧认证。保留原10 hPa有限顶界、原地形处理及未做整柱质量调整的定义；近地面残余约 `9.1–15.6e9 kg s^-1` 已记录，不能宣称底边界严格闭合。

文献背景：

- Nguyen et al. (2013), *The Hadley Circulation in Reanalyses: Climatology, Variability, and Change*. https://doi.org/10.1175/JCLI-D-12-00224.1
- *Hadley Circulation in the Present and Future Climate Simulations of the K-ACE Model*. https://doi.org/10.1007/s13143-021-00256-z

## 样本、维度和效应

固定模式顺序 GFDL-ESM4、NorESM2-LM、CESM2、CNRM-ESM2-1；前两者为 substantial recovery，后两者为 incomplete recovery，均为两模式等权。每个模式初末各10年。源年起点分别 `[1,271]`、`[1,271]`、`[1,131]`、`[1850,2120]`；这是原环流/降水输入的时间约定，不套用 EEI 的其他窗口。

年度降水原始8个文件均为 `(time,lat,lon)`，kg m^-2 s^-1 乘86400变 mm day^-1，年平均按真实月日数加权。Nor末期使用原先核验过的逐日转月产品。逐年插值采用原 n65 的周期经度、线性纬经度方法，非守恒，无极点外推。

保存 `annual_precipitation_samples.nc`：`pr(model,period,year_in_window,lat,lon)`、`day_weight(model,period,year_in_window)`；4模式、2时期、10年、89纬度（-88:2:88）、180经度（-180:2:178）。所有年度样本再按实际年日数加权，逐点复现旧十年降水气候态，误差小于 `1e-10 mm day^-1`。

流函数读取已有 n35 年度样本：`psi(model,period,year_in_window,plev,lat)`，只检验图示100–1000 hPa的12层和89纬度。`plev`实际存储 Pa；原更完整样本含17层10–1000 hPa。

每模式效应为日数加权的 later−early；组平均为两成员等权平均；第三行为 `Delta incomplete − Delta substantial`，不是两组绝对气候态之差。

## 显著性方法

固定四个模式，不重采样模式身份。在各模式各时期内对年度进行循环分块 bootstrap，同一套年份索引用于全部格点及两个变量；保留空间关联，不把连续年份当作独立模式。

- 9999次，随机种子20261027；块长分别2、3、4年，各重采样10年。
- 每个重采样均重新按所抽年份的实际日数加权，保持原效应的定义。
- 双侧中心化误差检验：`p=(1+count(abs(delta_star−delta_hat)>=abs(delta_hat)))/(B+1)`。
- 每格点取三个块长的最大 p，再实施 Benjamini–Hochberg FDR；阈值 `q<0.05`。
- 四个比较族分别校正：e流函数组差（1068点），b/d/f各自全球降水场（每幅16020点）。不将四个族合并，也未检验雨带峰位或ITCZ位置移动。
- 非有限或两时期均无年度方差的格点不作检验，不应当作“不显著”绘制。

显著格点数量：e为325/1068，b为1648/16020，d为7224/16020，f为4415/16020。这里是格点数，不是面积比例（全球格点面积随纬度变化）。数值和完整 p/q 均保存在 NetCDF 与 test_summary.csv。

方法背景：Wilks (1997) 的时间相关场分块重采样：https://journals.ametsoc.org/view/journals/clim/10/1/1520-0442_1997_010_0065_rhtfaf_2.0.co_2.xml ；Wilks (2016) 的气候场多重检验/FDR讨论：https://doi.org/10.1175/BAMS-D-15-00267.1 。本实现具体为上述中心化 bootstrap，并非声称逐项照搬论文算法。

重要限制：每端点仅10年，不能充分估计低频内部变率；循环连接和窗口内近似平稳是方法假设。每组仅两个固定模式，因此检验只描述这四个模式中的时间采样不确定性，不可外推为模式总体分组差异，更不能据此证明AMOC因果效应。不显著不代表无变化、等效或已恢复。

## 展示、版式与验收

显著处黑点，已检验但 q≥0.05 处淡灰交叉网纹；为可读性，显著点按固定格点索引抽稀（流函数每3个纬度；降水纬经均每3点），不按效应大小筛选，完整检验掩膜未抽稀。

雨带仅使用真实圆 `o`（初期）和叉 `x`（末期），路径已独立核对。原“方框”是相同位置圈叉叠加，并没有方形标记；现改为圆18 pt²、叉8 pt²。重叠处仍允许同心圈叉，不移动或错开真实格点；原10°经度显示子集、各36点及24条同经度连接全部不变。

在图的实际物理尺寸上，地图宽由3.464增至4.222英寸（约22%），两列间距由1.925降至0.451英寸；保持原投影、坐标范围、配色及色标限值。主图4680×4770像素、450dpi；PDF为748.8×763.2pt、字体嵌入，白色不透明。另一全球曲线PNG不变。

独立QC复核144个格点/块长 p值，最大差为0；四族BH、显著/不显著掩膜、点坐标哈希、全部原场/曲线/雨带/连接记录和源哈希均通过。已目视实际PNG和独立Poppler渲染PDF，面板/图例/坐标/色标无裁切，圆叉可辨；网纹反映不显著而不是缺测。尚未指定稿件内最终图宽，缩到单栏将影响细点和标记可见性，投稿前应以最终双栏或整页尺寸复审；不宣称已满足GRL最终制图要求或完整无障碍认证。

旧图件、QC及代码说明备份在 `archive/202610271913_circulation_precipitation_significance/`，清单含SHA256，按原相对路径保留。独立PDF审阅图也存于此处，正式figure仅保留当前活动版本。

## 工作流来源

scientific-visualization 技能用于数据不变校验、显著性与缺测区分、真实尺寸和独立PDF审阅；这不是气候结论的证据。软件/工作流参考记录：Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. https://doi.org/10.48550/arXiv.2609.00065 。2026-10-27核对arXiv当前记录；只记录于项目说明，不自动改动用户论文或Zotero库。
