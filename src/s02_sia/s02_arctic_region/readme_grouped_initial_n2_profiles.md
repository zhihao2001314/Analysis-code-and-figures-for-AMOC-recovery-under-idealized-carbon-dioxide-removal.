# n059：按用户修订三组叠加初始3月N²剖面

## 当前合图入口（2026-10-26）

用户已要求与MLD序列按海区合并。当前活动图为n073两张2×3六联图，见readme_grouped_stratification_mld_composite.md。上排N²仍逐值使用本目录结果和IAPv4，下排MLD观测改读n070 Argo，不用ELG。此前两张独立N²图已归档；以下保留原产品处理历史，科学CSV和QC不变。

2026-10-22显示名称统一为Better-simulated / Incomplete recovery、Better-simulated / Substantial recovery、Larger-bias models；不再将模拟表现名称限定为stratification。原成员CSV和数值不变，历史panel_title字段由绘图命名表覆盖，见`docs/group_naming_convention.md`。

## 本次明确采用的成员

2026-10-21用户明确指出UKESM层结模拟偏差较大，模拟较好的仅CESM2、CNRM、Nor、GFDL四个，再将这四个分两组，其余四个归同一组。因此本图固定为：

| 面板 | 英文标题 | 模式 | 数量 |
| --- | --- | --- | ---: |
| a | Better-simulated stratification / Non-recovery | CESM2、CNRM-ESM2-1 | 2 |
| b | Better-simulated stratification / Recovery | NorESM2-LM、GFDL-ESM4 | 2 |
| c | Larger stratification bias | ACCESS-ESM1-5、CanESM5、MIROC-ES2L、UKESM1-0-LL | 4 |

CAS不在本次八模式图中。成员由用户指定；本次没有重新检验局地深对流恢复，也没有依据新剖面重新筛选成员。第三组不标“全部深对流恢复”。这些模拟好坏标签对应用户的初始层结分组，不是本次新建立的统计技能检验，也不能把同一成员标签解释为Nordic海区单独检验得到的技能排名。未覆盖此前大气复合成员表、AMOC二分类表或旧模式剖面图。

## 初始窗口与观测时期

- 模式：初始1pctCO2第1–10年的3月，先逐年区域面积平均，再平均10个3月。
- **GFDL例外**：已有温盐剖面从实验第6年开始。用户明确批准本图沿用第6–15年10个3月，不伪称为第1–10年，不补造最初5年，也不改成只有5年的第6–10年平均。
- 观测：IAPv4的2005–2024年3月，共20年；与模式初始气候态不是同一时期。这是结构气候态参照，不是同时期偏差归因或因果检验。

具体calendar起止和成员写入 `membership.csv`：ACCESS101–110、Can1850–1859、CESM2 1–10、CNRM1850–1859、Nor1–10、GFDL6–15、MIROC1850–1859、UK1850–1859。

## 区域和输入

直接复用最新n039模式与n043观测产品：

- `output/s02sia/s02_arctic_region/p14_mld_fixed_regions_raw/coast_region_20261012_initial_n2_profiles.csv`
- 同目录 `coast_region_20261012_iapv4_profiles.csv` 及两者metadata。
- 区域为同目录 `mld_coast_regions_20261012.geojson` 中Labrador-Irminger、Nordic Seas两个手绘海岸闭合多边形；沿用最新区域，不退回旧的矩形选区。

Labrador–Irminger为本次主图，原配套Nordic–Iceland同样按固定三组输出对照。两个区域观测均使用对应区域曲线，在各自三个面板中重复叠加，不把两个区域观测混用。

## 维度与原始N²方法记录

本次只读取已经保存的CSV、JSON，不重新读取或推测NetCDF/MAT维度，不重新计算N²。实际保存维度：

- 原模式CSV：`model × region × depth_m`。含MMM的源文件只提取指定八个模式；MMM不作为成员使用。
- 当前 `member_profiles.csv`：8模式×2海区×100深度，共1600行。
- `group_mean_profiles.csv`：3组×2海区×100深度，共600行。
- `observation_profiles.csv`：2海区×100深度，共200行。
- `membership.csv`：8模式行，明确各自分组、绘图线型/颜色、成员和初始窗口。

既有源记录中的原始温盐维序：ACCESS、Can、UK为(time,lev,j,i)；CESM为(time,lev,nlat,nlon)；CNRM、GFDL、MIROC为(time,lev,y,x)；Nor为(time,rho,j,i)。MIROC原sigma-z深度按CF元数据重建，Nor由thkcello重建几何深度，沿用原产品，未新增处理。IAPv4原N2维度为(depth_mid,lat,lon)，原单位10^-4 s^-2，n043已按单位转为s^-2。

n039原方法为逐水柱温盐插值至0–1000m的10m几何深度层，不作垂向外推，使用TEOS-10的SA、CT、pressure和gsw.Nsquared计算N²，再以原生areacello对区域有效网格面积加权。IAP区域使用同一多边形与cos(latitude)面积权重，先逐年区域平均再平均20个3月；到共同深度的插值已在原产品中完成，本次不再插值。

默认使用原强迫数据，不做piControl漂移订正，不做时间或垂向平滑。原有效率和年际SD字段保留，但不额外筛选或作为置信区间绘制。

## 组平均方法

每个深度先读取各模式已经算好的区域初始N²气候态，再作**模式等权算术平均**：每组成员权重为1/2、1/2或1/4。不是先平均温盐再算N²，不是把各模式所有网格或所有年份混合加权，也没有读取旧八模式MMM充当组平均。

只有所有固定成员均有效的深度才输出组平均，避免随深度改变入组数。5m中点对应0–10m层，所有模式均缺失；保留NaN、不填值。因此模式和组平均有效深度为15–995m，共99层；观测为5–995m，共100层。观测在5m的值可以显示，但不拿这一缺少模式支持的层计算模式与观测的剖面诊断。Nor的Nordic负N²值原样保留，不截断为零。

## 图件与样式

用户明确要求仅保留三个面板，随后补充“模拟不好的就别截取掉了一部分”。现面板a–b（模拟较好两组）仍使用共同放大横轴，Labrador–Irminger为−0.02…0.90，Nordic–Iceland为−0.65…1.70；面板c（偏差较大组）恢复原完整横轴−1…7.1（单位均为10^-5 s^-2），完整展示全部有限成员值、组均值与观测值。三幅深度均为0–1000m，不截短深度轴。c标题增加“Full x-range”，图底明确提示横轴尺度不同；不应按曲线视觉横向宽度直接比较三组振幅，须读刻度。

彩色实线表示各模式，放大图增加稀疏空心圆/三角/方形等冗余点形。灰色细短虚线配空心菱形表示组平均，黑色实线表示IAPv4观测。均值线放在成员线后方，避免粗黑线遮挡。点形均选自已有真实深度节点，连续线保留全部原节点，不平滑、不垂向错位、不增加插值。成员CSV中的旧line_style字段保留为历史显示信息，现绘图线型覆盖记录在figure_manifest.json中。

前两组放大视窗之外的个别浅层观测值仍不显示，原始数值保留在CSV，未删除、未截断数值，也未当作缺失；第三组不再截去峰值。图底说明“Panels a–b are zoomed; panel c shows the full x-range. Note the different x-axis scales.”。图例放在各列下方，不覆盖剖面；无长方法脚注、无统计显著性或置信带。

图件目录：`figure/s02_sia/s02_arctic_region/p20_grouped_initial_stratification/n059_grouped_initial_n2_profiles/`。

- `labrador_irminger_initial_march_n2_grouped.{png,pdf,svg}`
- `nordic_iceland_initial_march_n2_grouped.{png,pdf,svg}`

物理尺寸10.8×5.6英寸；PNG为5400×2800像素、500 dpi、不透明白色背景。PDF、SVG为矢量曲线，SVG文字可编辑。为通用科研/PPT图，不宣称符合特定期刊的最终提交要求。

## 运行与QC

```powershell
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s02_sia/s02_arctic_region/n059_grouped_initial_n2_profiles.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'plot/S02_sia/s02_arctic_region/n059_grouped_initial_n2_profiles.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s02_sia/s02_arctic_region/check_grouped_initial_n2_profiles.py'
```

计算产品、分组、诊断、元数据、QC位于 `output/s02_sia/s02_arctic_region/p20_grouped_initial_stratification/n059_grouped_initial_n2_profiles/`。脚本均为顶层顺序分块，无自编函数、main、lambda或类。绘图只读output，不重算组平均。

QC独立采用权重矩阵乘法复算组平均，逐值核对原成员和观测没有改变、UK只在第三组、固定2/2/4成员、GFDL例外与其余初始十年口径、5m缺测和Nordic负值保留；检查来源与本次绘图输入SHA256、前两组共同放大横轴、第三组全部有限成员/均值/观测值在轴内、尺度差异明确提示、图像尺寸/分辨率、矢量文字与脚本计算绘图分离。两幅最终PNG均已人工查看，无标签裁切，外置图例不覆盖实际剖面。

2026-10-21可读性修订前，旧图件目录全部7文件及旧qc_report.json已移动至`archive/20261021100942_n2_profile_readability/`的原项目相对路径，并核对源路径不存在、归档目标存在、数量和SHA256一致；旧绘图、QC脚本和本文档也已复制保存并核对哈希。科学计算CSV与metadata.json未覆盖，本轮只更新绘图和QC。重新运行已有产品前仍须先归档对应正式输出。

随后按用户要求仅保留三个放大面板前，六面板版全部7个图件/manifest及旧QC已移动至`archive/202610211027_n2_three_panel_zoom/`，源、目标、数量和SHA256均已核对；绘图、QC脚本及本文档修订前副本同样已保存并核对哈希。科学数据继续原样保留。

本次第三组恢复完整横轴前，上一版全部7个图件/manifest及旧QC移动至`archive/20261021103222_n2_bias_full_range/`，源、目标、数量与SHA256均已核对；代码和文档修订前副本也已保存并核对哈希。科学CSV与metadata.json不变。

## 软件方法参考

Scientific-visualization技能用于原曲线保真、固定成员均值、共同尺度、颜色/线型冗余编码、缺失值和观测时期说明，以及实际导出图检查。软件参考：Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. https://doi.org/10.48550/arXiv.2609.00065 。本会话已核对arXiv记录，修订日期2026-09-02。
