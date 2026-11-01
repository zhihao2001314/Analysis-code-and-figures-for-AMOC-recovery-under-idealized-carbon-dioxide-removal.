# 两套海区分区的N²与MLD既有结果对照

## 当前选定版本 2026-10-27

用户已选择Rectangular regions矩形版。当前活动图仍在n078目录，stem为`labrador_irminger_9panels_rectangles`。本次按用户要求去掉总标题及顶部模拟表现、样本数分组栏，每个主面板右下角只标恢复类型：第一列Incomplete recovery，第二、三列Substantial recovery；第三列名称沿用用户此前指定，不重新计算或改变成员。图注仍须交代三列模拟表现和2/2/4成员数，以区分两个Substantial recovery组。

分区小图移到a图右侧偏上的空白区，相对父坐标范围为[.46,.26,.52,.64]；分区边界和科学数组不变。顶部留白收去1.2英寸，PNG为5670×5760、450dpi，PDF为907.2×921.6pt；九个数据面板的实际物理尺寸和距底部的位置不变。恢复文字8pt，第二行适当贴近下沿，实际曲线路径与标签字框交叉检查为零。白色细描边用于文字可读性，不用不透明色块遮住曲线。a–i仍与本幅实际图框左上角对齐，水平偏移0pt、向上3pt，12pt粗体。

当前核验入口为`check_rectangular_panel_labels.py`；独立结果为n077下`selected_rectangular_qc.json`及`manual_review.json`。核验实际SVG中的9个序号及9个恢复标签、字框范围、PDF中旧标题缺失、九个面板物理尺寸不变，以及25条科学曲线坐标哈希、三张地图、IAP灰带、源文件和分区几何不变。最终PNG和独立Poppler PDF预览均已目视；本次因新增标签和缩减画布，不再声称整图像素与旧版相同。

本次旧正式图与QC在`archive/202610271330_remove_overall_title/`，仅去总标题的中间版在`archive/202610271340_recovery_labels_layout/`，首次恢复标签排版在`archive/202610271350_recovery_label_spacing/`；均按项目相对路径移动并核对源、目标及SHA256。最终独立PDF预览位于最后一个归档目录。以下早期序号修改和两套分区对照记录保留为历史说明。

修改前矩形版、未选的海岸版及两份旧QC记录共10文件已移动至`archive/202610271300_rectangular_panel_labels/`，另复制2份代码/说明备份；保留项目相对路径，哈希及两端状态逐个核验。n078脚本当前只生成选定矩形版。以下两版对照流程与初次QC作为历史记录保留，不再把未选版作为活动图交付。底层两套科学数据均完整保留；初次`check_regional_definition_comparison.py`对应的是已归档的两版初始快照，不是当前序号布局的核验入口。

用户要求先展示两套分区各自对应的N²和MLD，第三行pattern不改，再比较分区效果。本次新增两张九联对照图，不覆盖n075原图，也不擅自确定最终采用哪套区域。

## 对照内容

1. Rectangular regions：上排读旧矩形N²及同矩形的IAP观测；中排读n061旧矩形MLD。
2. Coastal polygons：上排读现有n039/n043海岸多边形N²及同多边形IAP；中排读n035/n036海岸多边形MLD。

两版均保留各模式实线和灰色虚线组均值，IAP均值±1年际样本SD，横轴Year，无Argo参照，无海冰线。固定三组2/2/4成员，不根据新曲线重新分组。观测为2005–2024年3月，不是试验同期观测。SD不是测量误差、SE或CI。

第一行a、b共用−0.5–3.4，c共用−1–8.7（10⁻⁵s⁻²），深度0–1000m。第二行两版共用0–900m、Year1–280，覆盖两套数据全部曲线，不单独按各版最大值缩放。模式5m层结缺测和11年完整居中平滑的首尾缺口保留。N²初始窗口两版相同：除GFDL6–15年外，其余模式均前10个March。

a图右侧空白嵌入各自区域定位图。LI表示Labrador–Irminger，NS在矩形版为Nordic Seas，在多边形版为Nordic–Iceland sector（南伸至冰岛盆地）。边界直接使用既有顶点，不手工重新绘制，不与海冰线混淆，不计为第10个主面板。嵌图不遮挡任何模式、组均值或IAP SD范围。

第三行原样复用n074地图数组与显示样式。两版填色均为各组March MLD late−early，黑等值线均为初期组均值MLD，不根据选区裁剪或重算。固定共同±1100m色标。两版第三行与此前n075图逐像素核对。

## 重要局限 不属于严格的仅改变区域的敏感性试验

现有两套MLD产品并未完全统一算法和物理年份：旧矩形以最接近10m的有效原生层为参考，海岸版对sigma0插值到10m；两者都用增密0.03kg/m³的首次可解析交点。GFDL沿用既有旧索引/日历对齐方式，海岸源将up日历6–145标为index1–140，而CDR日历141–280标为index141–280；本次不声称修复了实际试验年的错位。矩形N²上界排除，旧MLD部分东西/南北边缘包含规则也不同。边界图只展示几何范围，不把两者冒充严格逐格相同掩膜。

因此两图可以查看各自既有区域诊断，但中排差异不能全部解释为分区本身造成。若要进行严格区域敏感性归因，须在相同原始场、10m处理、有效/底部受限规则和正确物理时间轴上统一重算两套MLD，再比较。不能仅凭某一分区让恢复分组分离得更明显就认定该分区更科学。区域还应有独立物理或地理依据。

## 流程与维度

- n076准备两套既有边界输入，仅复制output，不读取原始气候数据。
- n077整理各版剖面、观测和MLD，并计算固定成员等权均值；观测上下界为原均值±SD。
- n078只读取n077，绘制两张图；不重新计算N²、MLD、组均值、平滑或插值。
- check_regional_definition_comparison独立核验逐模式源数据、实际绘图哈希、权重均值、11年平滑、边界WKB、无曲线遮挡以及地图逐像素相同。

实际维度：version2×model8×depth100，观测version2×depth100；MLD version2×model8×year280；固定组均值version2×group3。共同地图group3×latitude141×longitude341（0.25°既有显示映射，不代表实际原生分辨率）。区域顶点longitude,latitude，矩形region2×vertex_order405，多边形FeatureCollection含2个Polygon。未重新推测NetCDF/MAT维序。

数据在output/s02_sia/s02_arctic_region/p23_stratification_mld_composite/n077_regional_definition_comparison/；图件在figure/s02_sia/s02_arctic_region/p23_stratification_mld_composite/n078_regional_definition_comparison/，每版PNG/PDF/SVG及独立manifest。均为新路径，无覆盖已有正式图。默认无piControl漂移订正。scientific-visualization技能用于共同坐标、保留原始曲线、揭示方法混杂、精确分区边界和实际导出核验；软件参考沿用n059记录，不更改论文参考文献。

## 验收

独立QC已通过：50组实际曲线坐标、各分区源曲线和观测、固定权重组均值、11年平滑、嵌图几何WKB及不遮挡数据检查均通过。两版第三行及原n075相同区域逐像素一致。两张实际PNG与两张独立Poppler PDF预览均已目视，无标题、坐标、图例或嵌图裁切重叠。记录见同阶段manual_review.json；预览位于archive/202610271230_regional_comparison_review/。该归档亦保存了澄清manifest措辞前的两份副本，源/目标哈希核对通过，图片未改。未依据哪版视觉分离更明显替用户选择最终分区。
