# n23：两类AMOC滞回路径示意（历史版本）

2026-10-26：用户确认改为现有八模式2/2/4三组，并按实际演变调整蓝线下弯。当前图件由n55计算、n56绘图生成，详见readme_amoc_three_group_schematic.md。此处两类图件已归档，计算输出和脚本保留历史复现，不再作为当前分组图。

用户在九模式恢复判据审查后选择只保留两类。当前概念图为 **Substantial recovery** 与 **Incomplete recovery**。英文、无模式名称、无长图注，采用统一尺度、细轴线、低饱和双色与浅灰路径间填色。没有针对具体期刊声明提交合规，“期刊式”仅指简洁矢量版式。

## 采用的定义

`R = mean(AMOC[271:280]) / mean(AMOC[1:10])`，使用未平滑的强迫年均AMOC；指标为最近40°N纬度、500 m及以下Atlantic–Arctic翻转最大值，仍按默认不做piControl漂移订正。

- **Substantial recovery / 基本恢复**：R≥0.8，表示末期达到至少初始强度的80%，不等于完全恢复、不等于精确闭合，也不要求无滞回。
- **Incomplete recovery / 末期未恢复至门槛**：R<0.8，仅描述第280年结束时的窗口状态，不表示之后不会回升、不表示永久不可逆。

按九模式原始年序列审查，本阶段快照成员为：基本恢复ACCESS、CanESM5、CAS、GFDL、Nor；末期未恢复UKESM、MIROC、CNRM、CESM2。来源及数值保存在 `amoc_endpoint_two_class_membership.csv`。这只是本概念框架采用的纯AMOC端点二分类，不覆盖此前高初始层结筛选表、不改变既有四模式大气复合；后续排除CAS的分析须另行明确样本数。

80%是当前研究的操作性门槛，而非物理临界值或参考文献给定的通用界限。依据n22，80%分类对10/20/30年端点窗口稳定；初始早期仍是变化强迫下的参考窗口，不等于独立的平衡piControl。

## 示意图如何制作

计算侧 `src/s01_amoc/n23_amoc_two_class_schematic.py` 读取n22已经保存的端点诊断，另行生成手工概念坐标。二维坐标维度为 `panel × drawing node`，箭头维度为 `panel × phase`，成员表为model行。未读取新的NetCDF/MAT，不重新提取AMOC，也不对模式轨迹拟合或求组平均。

两图共用相同的初始CO₂增加分支，表示AMOC总体减弱；降低CO₂分支均先滞后减弱、随后回升。左图绘图末点选为初始的0.92，右图选为0.43，分别位于0.8界限之上与之下。这两个数字仅用于清楚表达概念，不是组均值、估计恢复率或具体模式结果。两图均保持非零路径分离；左图也不强行闭合到初态，避免把“基本恢复”错误画成“完全可逆”。数学形状和坐标完整保存在output元数据与CSV。

2026-10-21按用户要求改为圆润的近椭圆环路，随后依据用户指出的蓝线上凸问题修正了上弧。当前参数theta∈[0, pi]，`x=(1-cos(theta))/2`；上弧改为`0.28+0.72*sqrt(1-x**2)`，即以起点为最高点的椭圆四分弧，解析保证最大值为1且全程单调下降。下弧仍为`1-0.72*x-0.22*sin(theta)-gap*(1-x)`，gap为0.08或0.57。没有截平、裁切或隐藏超出1的部分；通过更换概念几何本身消除假初始增强。整体保留椭圆感，但不再宣称上下两弧属于同一个严格椭圆。分类、两个概念末点和实际模式数据均未改变。

绘图侧 `plot/s01_amoc/n23_amoc_two_class_schematic.py` 只读取保存的坐标、箭头、名称和阈值。蓝灰实线为CO₂增加，朱红虚线为CO₂降低，箭头沿时间方向；空心圆为初始状态，方形为概念末期状态，黑小点为共享CO₂峰值。浅灰填色表示路径分离，不是置信区间；0.8灰虚线为采用的端点门槛，1.0灰实线为初始强度参考。横轴C0、Cmax为概念浓度端点，没有添加实际ppm数值。R的正式定义在图外说明，不在图内追加方法段落。

## 输出与核对

计算结果：`output/s01_amoc/n23_amoc_two_class_schematic/`。正式图件：`figure/s01_amoc/n23_amoc_two_class_schematic/amoc_two_class_hysteresis_schematic.{png,pdf,svg}`。

物理尺寸7.2×3.45英寸；PNG为4320×2070像素、600 dpi、白色不透明背景，另有PDF与可编辑文字SVG。期刊与最终版宽尚未指定，不能声称已满足某一期刊最终要求。

执行顺序：

```powershell
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s01_amoc/n23_amoc_two_class_schematic.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'plot/s01_amoc/n23_amoc_two_class_schematic.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s01_amoc/check_amoc_two_class_schematic.py'
```

脚本拒绝覆盖既有正式文件，重做须先归档。QC核对1202个概念节点、4个箭头、九模式5/4成员、共享转折点、末点与80%定义一致性、上下弧各自的椭圆几何恒等式、蓝线最大值严格不大于1且单调下降、无自编函数、英文SVG及PNG尺寸/分辨率/背景。最终PNG已人工查看，蓝线无上凸，标签无裁切，尺度统一，近椭圆环路与两种缺口清楚。

2026-10-21更新前的5个output文件、4个figure文件及4份源码/说明快照保存在 `archive/20261021092025_amoc_elliptical_schematic/`，保留原项目相对路径。移动后核对源路径消失、目标文件数量及全部SHA256；源码副本也逐一核对哈希。旧版本可恢复，当前图件目录仅有更新后的单一活动版本。

2026-10-21消除蓝线上凸之前的对应5个output、4个figure及4份源码/说明快照保存在 `archive/20261021092553_amoc_schematic_no_overshoot/`，同样逐文件核对数量、源路径与目标路径及SHA256。档案保留错误形状的历史记录，不作为当前图件使用。

n21旧三类概念图在新图验证后移入archive，旧计算坐标与脚本保留作可追溯记录，活动图件只保留n23两类版。n20各模式实值曲线、n22诊断图是不同用途的科学产品，本次未改。

## 参考和软件

概念框架参考：Kim, S.-K., et al. (2022). *Widespread irreversible changes in surface temperature and precipitation in response to CO₂ forcing*. Nature Climate Change, 12, 834–840. https://doi.org/10.1038/s41558-022-01452-z 。借用路径分离与端点恢复分开描述的思路，不借用其具体空间分类方法或赋予80%普遍含义；当前会话2026-10-20已核对原文。

Scientific-visualization技能用于概念与数据分离、统一尺度、颜色加线型编码以及图件导出检查。方法软件参考：Kassis, T., Agarwal, V., He, Y., Patel, D., & Brueckner, A. M. (2026). *Scientific Agent Skills: A Library of Procedural Knowledge for Research Agents*. https://doi.org/10.48550/arXiv.2609.00065 。当前会话2026-10-20已核对arXiv记录，修订日期2026-09-02。
