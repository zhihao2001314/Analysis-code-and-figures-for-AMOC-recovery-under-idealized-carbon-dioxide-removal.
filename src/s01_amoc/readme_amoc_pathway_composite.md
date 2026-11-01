# n62 / n63：九模式异常、面积概念与两类恢复路径八面板

## 图a改为26.5°N（2026-10-29，当前活动版本）

按作者确认，图a采用与RAPID对应的26.5°N目标纬度，使用各模式最接近该目标的原生网格行。新准备入口`src/s01_amoc/n66_amoc_overview_26p5.py`只读取既有`output/amoc_s01_n01/l_26.5/`的Forced AMOC strength年度原值，不读取原始NetCDF/MAT、不使用piControl趋势或控制态均值。先减各自原始初始1–10年均值，在完整来源年度序列上做原定义11年居中平均、min_periods=1，再显示1–280年；九模式等权均值和总体SD随之更新。

当前绘图入口仍为`plot/s01_amoc/n63_amoc_pathway_composite.py`，读取`output/s01_amoc/n66_amoc_overview_26p5/`；正式图件仍在`figure/s01_amoc/n63_amoc_pathway_composite/amoc_overview_and_pathways_8panels.{png,pdf,svg}`，路径在绘图参数区列出，保存后打印。横纵轴、版式、CO₂强迫及b–h输入逐字节保留。**d/e/g/h实际路径及原恢复分组仍采用40°N定义，本次只换图a**；当前caption.txt/metadata明确这一区别。图a数据范围−12.716至2.831Sv，完全位于保留的−21至3Sv纵轴内。

新CSV实际维度为simulation_year280×model9，另附CO₂、MMM、总体SD和上下界。原年度输入长表按Model/Simulation year转宽表；无新NetCDF维序假设。`panel_a_latitude_mapping.csv`记录各模式实际原生纬度、变量、盆地、单位和原索引来源；目标26.5°N并非将每个模式空间插值到恰好26.5°N。图a端点表单独保存，不据此重分组或改下排40°N路径。现有质量流函数沿用除以10⁹kg/s的等效Sv换算约定。

旧n63图件、清单及编辑前代码/说明已在`archive/*_amoc_panel_a_26p5/`按相对路径归档并核对SHA256。下文有关图a40°N的数据说明属于旧版历史，以本节为准。

## 当前图例内置（2026-10-27）

依用户明确要求，a的九模式全称和9-model mean以紧凑两列放在左下角图内空白区，b的CO₂增减两行图例放在示意图左下角，无白底遮罩。保留加高首行、全部科学坐标、轴范围及其他面板；不重算数据。记录实际图例的画布边界及数据坐标边界，并核对图例矩形不与任何模式曲线、均值或SD边界交叠。旧图/QC和代码说明备份归档`archive/202610272213_amoc_legends_inside/`，源目标和哈希核验。scientific-visualization技能用于遮挡检查和实际PNG/PDF审阅，沿用文末来源记录。

## 当前c曲线舒缓修订（2026-10-27）

仅按用户要求修改c概念示意锚点，使两分支沿CO₂更均匀、舒缓地下降，减弱中段陡降和尾端回弯。保持初始/末期/峰值、共同CO₂范围、所有轴范围和布局；a/b/d/e/f/g/h曲线与箭头完全不变。此变更只属于理想化示意，不修改真实模拟数据。独立QC比较新旧c的最大局部斜率及端点，确认两阶段局部斜率峰值都降低；不是对模式资料做更强平滑。上一版全部产品及代码说明快照核验归档`archive/202610272209_amoc_gentler_c_schematic/`。scientific-visualization技能用于示意与数据区分及PNG/PDF复核，来源引用沿用文末。

## 当前首行加高与下置图例（2026-10-27）

首行a/b绘图区高度由0.205画布高度增至0.265，约加高29%，横向宽度和全部轴范围不变，不压缩纵轴来夸大响应。a模式图例改为下方左侧留白内的两行五列，b阶段颜色图例同步下置，避免覆盖曲线。下两行全部面板的坐标、位置和示意面积形状不变，计算产品未重算。上一布局四份图件/清单和QC/目视记录已核验归档`archive/202610272204_amoc_taller_first_row/`，源目标状态和哈希逐文件检查。scientific-visualization技能用于版式与实际PDF审阅，沿用文末工作流来源记录。

## 最新对称轮廓和端点标注（2026-10-27）

当前b使用中心对称轮廓：601节点上逐点满足`M_up(x)+M_down(1−x)=1.30`，保留上下错开的短平台，避免上版长平台/短平台不均衡。两端明确标为Initial CO₂、Peak CO₂，峰点另用深色圆点。填色内部标为Hysteresis area (A)，中文用“滞回面积”比“滞后区域”准确；这里特指暂态分支分离，不是以年计的时间滞后，也不声称平衡多稳态。下文A_path是历史命名，现统一A。

f高位平台由归一化CO₂约0.50缩短至0.08，与c一致，之后平顺下降，低位返回线保持规整；仍保留原末期恢复率和CO₂峰值的状态值。c及所有真实模式/九模式概览均不变。旧产品和代码说明快照已核验归档`archive/202610272157_amoc_symmetric_hysteresis_guides/`。此次scientific-visualization技能用于中心对称、标签留白和导出复核，工作流引用沿用文末记录。

## 最新手绘形状修订（2026-10-27）

用户明确要求b为偏菱形/斜四边形，c为窄而倾斜的近闭合路径，f为“高位平台—下降—低位返回”的开口路径。当前b/c/f全部改为按手绘概念构造的PCHIP锚点示意，不再沿用n55/n57数据引导曲线形状，也不将平台转折解释为模式临界CO₂。c/f只保留原组平均十年恢复率端点及原共用峰值，形状本身不拟合模型。明确标注Conceptual/Schematic；下文原示意算法段落均为首版历史，以本段及新metadata为准。锚点保存为`schematic_anchors.csv`，不可混同为模拟年份或观测。

a的年度异常表逐字节不变，d/e/g/h的实际曲线及箭头不变；所有路径保持实线，无总标题。旧n62全部结果、旧n63图和代码/说明快照已逐文件核验归档到`archive/202610272144_amoc_diamond_schematic_revision/`。沿用scientific-visualization技能核对示意/真实数据的区分及实际PNG/PDF导出，工作流来源仍保留在本文末段。

最终b只在填色内部标`A_path`，全称与定义放图注，避免说明文字压住轮廓；中间标签版归档`archive/202610272150_amoc_area_label_clearance/`。独立QC及最终PNG/独立PDF目视通过，当前manual_review记录对应本次手绘形状版。

## 设计与范围（2026-10-27）

用户指定首行a占两列放九模式AMOC异常、b放面积示意；第二行c/d/e为一类示意和两个对应模式，第三行f/g/h为另一类，共8个字母面板。没有总标题；必要模式名、坐标、方向图例及Schematic/Conceptual标识保留。蓝/陶红两条CO₂分支均为实线，并以箭头区分时间方向。这里只是新布局预览，不新增模式分组、不根据新图重选模式、不做漂移订正或显著性检验。

原三行九面板中第三个larger-bias业务分组未作为第三种滞回类型保留；其成员仍全部出现在a的九模式概览。下两行只展示原先better-simulated的四模式：Substantial recovery为GFDL-ESM4、NorESM2-LM；Incomplete recovery为CESM2、CNRM-ESM2-1。原恢复率R为末期271–280年/初期1–10年，阈值0.8；不把临时绘图命名改成全新物理分类。

## 数据及维度

计算`src/s01_amoc/n62_amoc_pathway_composite.py`只读取既有output，无原始NetCDF/MAT读取。a取`output/amoc_s01_n01/l_40/amoc_initial_reference_plot_data.csv`的1–280年：约40°N最近模式纬度、Atlantic–Arctic盆地、500m以下最大翻转。异常为原强迫值减每模式自身原始初期十年均值，不使用piControl趋势或控制态均值。这与早期piControl-relative同名anomaly图不同，当前采用初始参照以保持与下方归一化路径的参照一致。

a复用原产品11年居中平滑；该产品在完整上游年度序列上平滑，本次只裁取显示时间，不重算年280附近窗口，因此边缘可能使用280年之后已有的真实年份。d/e/g/h复用n22在原280年序列上的既有11年平滑，末端为缩短窗口；此差异属于原产品的平滑支撑范围，不合并、补造或强行匹配边缘。黑线为九模式等权平均，灰带为模式间总体标准差±1SD，不是CI。原始序列按`Simulation year × Model`转宽表；a保存为`simulation_year280 × model9`的宽CSV（另含MMM/SD/上下界/CO₂）。

c/f完全复用n57前两行601节点PCHIP数据引导示意，横坐标由归一化浓度转为物理ppm；概念初态为1，概念末态取组成员十年恢复率均值，不能把它称为实际瞬时终点。d/e/g/h完全复用n57四模式路径：每模式升碳140点、降碳141点，转折年140共用，不改曲线形状、不额外平滑，保留原有非单调变化。路径长表维度`panel7 × phase2 × node`，列含年、CO₂ ppm和相对AMOC。箭头表`panel7 × phase2`；降碳箭头由高浓度指向低浓度。c/f展示标识避免与真实模式数据混淆。

## 面积命名与解释

b为独立解析概念图，不是模式拟合、原论文数据数字化或平衡解。两分支在起点/CO₂峰值闭合只为解释面积，不能据此说实际模式都完全恢复。初态1、峰态0.34，蓝分支`1−0.66*x^0.75`，陶红分支再减`0.25*sin(pi*x)`，x为0–1归一化浓度；无经验面积统计被重新估计或列成真实模式指标。

图中用 **Path-dependence area, A_path（路径依赖面积）**。可在正文限定为“transient hysteresis area（暂态滞回面积）”，但不直接写成恢复时间或临界点证据。定义为共同线性CO₂区间上两分支绝对差的积分：`A_path = integral |M_up(C)−M_down(C)| dC`。本图纵轴若为相对AMOC、横轴ppm，面积单位为相对AMOC×ppm；若要跨研究量化比较，须先明确标准化方式和共同积分范围。面积衡量路径分离，并非以年计的时间滞后，也不能区分慢调整与真正平衡多稳态。末期是否恢复与路径面积是不同问题；两个真实路径即便交叉也不可人为修成闭环。

图外完整英文图注草稿在`output/s01_amoc/n62_amoc_pathway_composite/caption.txt`及metadata；尚未编辑Word或Zotero。

## 导出与验收

绘图`plot/s01_amoc/n63_amoc_pathway_composite.py`仅读新output，图件在`figure/s01_amoc/n63_amoc_pathway_composite/`。10.2×8.7英寸、600dpi，PNG6120×5220；PDF734.4×626.4pt、TrueType嵌入，SVG文字保留。没有总标题，不声称已按最终GRL版心排版。

独立核验脚本`check_amoc_pathway_composite.py`从原年度Forced强度列重新计算a的基准/平滑/平均/SD，并复核六个原路径面板逐点不变、示意末态均值、24条实际绘线全部实线、14个箭头方向和来源、轴内完整覆盖数据/灰带、源SHA256、字体与尺寸。初次试版因a纵轴低端截断而未采用，图件与PDF预览已核验归档`archive/202610272135_amoc_eight_panel_first_render/`，最终a范围−21至3Sv完整包含最低−19.374Sv。

scientific-visualization技能用于区分示意/数据、同轴尺度、面积语义、透明记录基准与平滑支撑范围、实际PNG和独立PDF渲染审阅；软件工作流引用沿用阶段已有Kassis等（2026）Scientific Agent Skills，https://doi.org/10.48550/arXiv.2609.00065，不当作气候物理证据，不写入用户论文文献库。

最终验收：独立QC通过，实际最终PNG与独立Poppler PDF渲染目视无裁切；面板字母和图例完整。被替代的n58旧九面板四个图件/清单已按原相对路径移至`archive/202610272137_amoc_nine_to_eight_panel/`，源端为空、目标数量及SHA256逐文件核验；旧n57计算结果和n58绘图代码不改，随时可恢复旧图。当前活动交付为n63八面板。
