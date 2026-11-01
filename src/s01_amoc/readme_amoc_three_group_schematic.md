# n55 / n56：三组AMOC–CO₂路径示意

当前正式图已按用户要求改为n57/n58三行九联图，见readme_amoc_schematic_model_panels.md；旧n56三联图和其QC/临时诊断归档至F:/sea/Antic sea ice constraint/archive/20261026141408_amoc_nine_panel_layout/。n55六份科学坐标/成员/元数据继续作为九联图输入，未改变；历史脚本保留，不重新运行覆盖结果。

用户要求将旧n23两类示意图改为现有八模式三组，并将原上拱椭圆蓝线调整为更接近实际演变的下弯路径。三组经用户确认是GFDL/Nor、CESM2/CNRM、ACCESS/CanESM5/MIROC/UK，CAS不加入。不修改任何既有分组注册表。

## 分组不是三种恢复结局

a Better-simulated / Substantial recovery：NorESM2-LM、GFDL-ESM4，2个；各自R≥0.8。

b Better-simulated / Incomplete recovery：CESM2、CNRM-ESM2-1，2个；各自R<0.8。

c Larger-bias models / Mixed recovery：ACCESS-ESM1-5、CanESM5、MIROC-ES2L、UKESM1-0-LL，4个；前两者达到恢复门槛，后两者没有。第三组的筛选依据不是AMOC恢复率，不能把它统一画成一种恢复结局。

## 数据与示意几何分开

本次读取n22已经核查的nine_model_normalized_series.csv、nine_model_endpoint_diagnostics.csv及当前n49分组注册表。输入维度model×simulation_year，保留1–280年；未读新NetCDF/MAT、未改原始模拟、未新增漂移订正。

AMOC指标仍为最近40°N、深度≥500m的Atlantic–Arctic翻转最大值。先对每个模式以自身原始初始1–10年均值归一化，再对同组模式等权平均。曲线取用n22已有11年居中平滑；不重新选择平滑窗口。

横轴为线性归一化浓度x=(C−C0)/(Cmax−C0)，C0=284.7ppm，Cmax=C0×1.01^140。不是时间轴，也不是logCO₂。

上升阶段在x=0、0.12、0.28、0.48、0.72、1设置稀疏锚点；x>0锚点由实际平滑组路径按浓度线性取样，x=0取概念初态1。降低阶段同样取样，增加实际平滑组最低点对应浓度作为锚点；x=1与上升分支共用转折点，x=0取真实成员R的平均。PCHIP只用于这些锚点之间的简化绘图，不被当作新模式数据。

R=原始271–280年均值/原始1–10年均值。三组平均成员R分别约0.9304、0.3554、0.8056。第三组成员范围0.5568–1.0203；灰色竖线及短横刻度表示四模式实际末期R范围和各成员，**不是置信区间**。红方块表示成员R的平均，不是组平均AMOC之比，也不是实际第280年瞬时值。

这是一张**Data-informed schematic**，不是完整的模式逐年轨迹，也不是拟合出的平衡分支。第一点x=0的归一化1和两个十年端点R是示意端点口径，不能假称它们都是C0下的同步瞬时模拟值。实际逐年曲线仍由n20/n22科学产品保存，未覆盖。

## “下凹”的口径

不再采用n23的椭圆上弧。采用实际数据锚点后，蓝线中段低于初态与CO₂峰值之间的直线，表达较早明显减弱、随后下降趋缓。未强行给所有模式或整个浓度范围施加同号二阶导数，也不将局地弯曲差异改成直线来追求外观。

红线沿时间从右向左，均保留降碳后继续减弱、随后回升的有限时段响应。箭头沿时间方向。浅灰填色只是分支分离，不是置信带、显著性或新增环面积指标；未声称有限实验能证明永久不可逆。

## 文件与核查

计算：src/s01_amoc/n55_amoc_three_group_schematic.py。

绘图：plot/s01_amoc/n56_amoc_three_group_schematic.py，只读output。

独立QC：src/s01_amoc/check_amoc_three_group_schematic.py。复核实际组路径锚点、端点均值和范围、2/2/4成员、无CAS、蓝线单调下降且初始不超过1、中段下弯、红线最低点/箭头、PNG实际尺寸、PDF嵌入字体和SVG文字。

结果：output/s01_amoc/n55_amoc_three_group_schematic/，含稀疏锚点、1803个示意节点、6个箭头、8成员表、真实端点范围、元数据及QC。

图件：figure/s01_amoc/n56_amoc_three_group_schematic/amoc_three_group_hysteresis_schematic.png/pdf/svg。当前精修版7.5×4.0英寸（190.5×101.6mm），600dpi PNG4500×2400。目标期刊/最终宽度未指定，属于通用科研工作稿。

旧n23两类图在新图验收后归档；旧计算输出和脚本保留历史复现，其他实际模式路径、AMOC/温度图及大气复合均不动。

Scientific-visualization技能影响共同尺度、示意与实值的区别、第三组混合结局和范围而非CI的表达、颜色加线型与实际导出核查。软件引用沿用本阶段readme_amoc_two_class_schematic.md已有记录，不新增未经核实的期刊绘图要求。

最终QC通过：实际组路径锚点独立复核最大差1.11×10⁻¹⁶，2/2/4成员及第三组2恢复/2未恢复均核对；PNG与独立Poppler渲染的PDF均已目视，未见裁切，PDF字体嵌入、SVG文字可编辑。QC使用viper环境的SciPy和配套Python的纯Python pypdf：以viper启动后将配套Python的Lib/site-packages追加到sys.path，再用runpy运行check脚本，无需安装或修改任一环境。

旧n23四份图件/清单及一份已审阅PDF预览移动归档，并复制保存旧README；源目标状态、5个移动文件及备份文件SHA256均已核对。归档清单位置：F:/sea/Antic sea ice constraint/archive/202610261309_amoc_three_group_schematic/archive_manifest.json。

## 2026-10-26：仅精修图形表达

用户要求更精致的版式。本次未重新计算曲线，6份科学CSV/元数据逐字节与精修前一致；三个面板的xlim(-0.07,1.06)、ylim(0,1.12)、所有1803个节点、6个箭头、分组和端点都不变。

改为Arial字体，统一总标题与图例；每幅依次排列分组依据及n、恢复状态、实际模式名称。蓝色#2B5F7A、朱红#BE553D，灰填色#F0F3F5；保留实线/虚线、初始圆点/末期方块的冗余编码。y轴增加0.2间距刻度，不改变尺度。CO₂下标使用mathtext，避免Arial原生缺字。

第三组灰色范围符号仅在图面上向左错开至x=-0.035，以免遮住x=0的初始圆点；纵向范围和每个成员的R不变，并用细线连接至末期均值。这个横向偏移不是另一个CO₂浓度，不用于任何计算。

QC通过：全部模式名存在于实际PDF文字中，字体嵌入且非Type3；PNG、独立PDF渲染和灰度诊断已目视，无裁切和下标缺字。蓝/红对白色背景对比度约6.95/4.61，对灰填色约6.24/4.14；这只是对比度与冗余编码筛查，不是完整色觉或期刊合规认证。

精修前6个图件/QC移动归档、3份代码和说明备份均核对哈希，位于archive/20261026132934_amoc_schematic_visual_refinement/。首轮字体缺字试稿4图件另归档至archive/20261026133201_amoc_schematic_subscript_repair/。活动figure目录仅保留最终PNG/PDF/SVG和清单。

## 2026-10-26：移除顶部总标题

用户随后要求删去各面板具体模式名；当前图内已不画8个模式名称，仅保留分组、n和恢复状态。具体成员在本文分组列表、membership.csv及图件manifest中保留。旧带模式名版本归档至archive/20261026140746_amoc_remove_member_labels/，所有科学输入和路径节点不变。

应用户要求，删去顶部AMOC response to CO₂ reversal及右上Data-informed schematic；保留三幅分组标题、模式名称、图例、坐标与所有曲线。Data-informed schematic意为“依据数据绘制的示意图”：从实际组平均路径抽取稀疏锚点，用PCHIP连接概括演变，而不是完整逐年模拟轨迹或平衡滞回分支。此科学口径保留在本文及图件manifest，不能因为图内文字移除而改变其解释。

本次仅删两处文字，六份科学输入不重算、不改数值。旧图件/QC及代码、说明备份归档至F:/sea/Antic sea ice constraint/archive/20261026140609_amoc_remove_top_title/，两端状态及SHA256已核对。
