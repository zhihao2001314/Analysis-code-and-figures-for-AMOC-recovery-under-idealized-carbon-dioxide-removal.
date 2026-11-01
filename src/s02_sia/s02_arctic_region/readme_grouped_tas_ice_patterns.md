# n066：三组气温变化与海冰初末边缘叠加

## 重要：旧三联图已由经端点核验的六面板替代

后续制作全年北半球降水—急流六面板时，查明本旧版依赖中的UKESM气温CDR窗口实际错位，以及CanESM5初始SIC旧CSV与原始120个月不匹配。两项已在`src/s05_atmosphere/p05_amoc_recovery_composites/n37_annual_nh_precipitation_jet.py`及其n38六面板中修正，详见同目录`readme_annual_nh_precipitation_jet.md`。旧n066科学数组和代码保留供追溯，不代表可继续作为已核验的同窗结果。

旧三联图PNG/PDF/SVG及manifest共4文件已移至`archive/20261021163828_legacy_pattern_cache_superseded/`，本说明旧副本1文件同时保留，5文件源/目标、数量和SHA256均核对。当前图件位于`figure/s05_atmosphere/p05_amoc_recovery_composites/n38_grouped_climate_precipitation_jet/`；不要重新运行本旧绘图脚本恢复成未修正的正式图。以下章节为旧版方法和样式历史记录。

## 图的目的与分组

用户决定暂不继续混合层观测核查，将原海冰、海温、气温pattern合并表达，随后明确只保留气温填色。因此本次只输出一张三个面板的气温—海冰复合图，不再增加海温面板。

| 面板 | 英文分组 | 成员 | n |
| --- | --- | --- | ---: |
| a | Better-simulated stratification / Non-recovery | CESM2、CNRM-ESM2-1 | 2 |
| b | Better-simulated stratification / Recovery | NorESM2-LM、GFDL-ESM4 | 2 |
| c | Larger stratification bias | ACCESS-ESM1-5、CanESM5、MIROC-ES2L、UKESM1-0-LL | 4 |

固定成员来自n059用户已确认的三组；CAS不参与。组名不是本次新做的统计技能检验，也没有将第三组重新认定为全部深对流恢复。n是模式数，不是年度独立样本数。

## 图外解释：填色、线条与时间

- 范围：45–90°N，北极立体投影，经度0°在图下方，180°在上方。
- 填色：组平均近地面气温tas变化，later减initial，单位°C；赭红色增暖，蓝色冷却，米白色接近零。陆地和海洋上的气温都显示，不用灰色陆地覆盖气温场。
- 三组共用−7至+7°C线性色标；本次三组所有有限值均在范围内，无数值因色标被截去。组a范围约−6.20至+3.34°C，组b约−1.62至+1.56°C，组c约+0.51至+6.64°C。
- **当前版本**灰色虚线：initial年平均海冰浓度15%边缘；黑色实线：later年平均海冰浓度15%边缘。线条外的细白色描边只用于提高对比，不表示误差或不确定性。该线型与归档旧版的黑色初始/绿色末期不同，不应混用图例。
- initial：1pctCO2实验第1–10年；later：1pctCO2-cdr降碳实验第131–140年，对应名义拼接实验第271–280年。所有场沿用既有年平均端点，不混用3月气候态，也不读取超长CDR文件的最后10年。
- 这不是3月海冰边缘，也不是某一年或冬季最大海冰范围。浓度达到15%的一侧为达到该阈值的海冰区，但该线不展示区内浓度变化的幅度，也不能替代ΔSIA面积图。
- 图内只保留标题、必要图例和色标，不放计算方法长脚注；本文件提供完整图外说明。没有观测边缘、显著性点、置信带或误差棒。

## 计算步骤与缺测口径

1. 读取各模式保存的`raw_difference`气温场，逐值验证等于`cdr_final-co2_initial`，不用`corrected_difference`。默认保持原强迫数据，不进行piControl漂移扣除。
2. 海冰读取原三变量pattern使用的CSV中的`edge_raw_co2_sic_fraction`和`edge_raw_cdr_sic_fraction`，保留原源值；不改用已订正海冰列或新建的其他海冰插值版本。
3. 温度逐格点按模式等权算术平均；只有所有固定成员均有效才保存组均值。海冰两端也先平均浓度，并进一步使用initial/later共同支持。
4. 从组平均海冰浓度取15%等值线；不是先对每个模式取边缘再平均位置，不是成员边缘的并集或多数表决。组平均边缘不代表每个成员都在该位置有海冰。
5. 原CSV中未提供的海冰格点保留NaN，不沿用旧绘图将其补0的步骤；不把插值缺口解释成无冰。为避免海冰轮廓穿过制图陆地，单独保存陆地掩膜后的轮廓显示场，不修改原始成员或组平均浓度。
6. 没有新增时间平滑或空间平滑，没有重新计算原生网格MLD、SIC或tas。图件保留上游原有1°重网格结果；不是原生模型网格轮廓。循环经度列只用于闭合±180°绘图接缝，不新增观测或模式数据。

注意：上游tas的`raw_difference`虽然没有扣除漂移，其保存时使用了CO2、CDR、piControl共同有效支持。本次沿用该支持，不把这一点误称为重新计算的纯CO2/CDR双阶段掩膜；当前三组气温全部16200格点有效。旧N²剖面中GFDL初始6–15年例外不适用于本次气温和海冰端点。

## 输入、维度与输出

- 气温：`output/s03tem_n01/arctic_tas_<model>_pattern_1deg.npz`，以及同目录summary/QC；所有本次8成员均为standard_tas。
- 海冰：`data/s02_sia/s02_arctic_region/arctic_sia_pattern_values_1deg.csv`。
- 分组：`output/s02_sia/s02_arctic_region/p20_grouped_initial_stratification/n059_grouped_initial_n2_profiles/membership.csv`；只提取模式、组别和排序，不继承其温盐初始年份。
- 底图：n067从本机已缓存的Natural Earth 50m海岸线制作GeoJSON，绘图只读其output；仅绘制细海岸线，不覆盖气温。原n066科学数组中的110m陆地轮廓显示掩膜保持不变，未用新海岸线重新筛选海冰或重算任何均值。
- 实际维序：上游tas数组为`latitude45,longitude360`，lat/lon坐标均为二维45×360；SIC CSV按明确经纬度重索引。纬度45.5…89.5°N，经度−179.5…179.5°E，间隔1°。
- 本次NPZ成员数组为`model8,latitude45,longitude360`，组数组为`group3,latitude45,longitude360`。温度单位°C差值，SIC单位0–1比例。没有本次新增的NetCDF/MAT读取或按位置猜测维度。
- 计算结果：`output/s02_sia/s02_arctic_region/p22_grouped_temperature_ice/n066_grouped_tas_ice_patterns/`，包括成员及组场、分组、支持诊断、来源哈希、QC与底图。
- 图件：`figure/s02_sia/s02_arctic_region/p22_grouped_temperature_ice/n066_grouped_tas_ice_patterns/grouped_air_temperature_sea_ice_edges.{png,pdf,svg}`及manifest。
- 导出：13.2×5.25英寸，PNG7920×3150像素、600dpi、不透明白底；PDF/SVG气温栅格化而文字/线条为矢量，SVG文字可编辑。未指定目标期刊，属于通用科研展示版，不宣称特定期刊合规。

## 当前视觉重设计：完整北极、简洁配色与单行标题

用户反馈前一版整体风格仍不好看，并明确选择保留完整北极三联图、全部经度；因此没有改成北大西洋局地图，也没有裁掉北太平洋。投影仍为中央经度0°的北极立体投影，范围保持45–90°N。

当前配色为显式记录的蓝—米白—赭红连续发散色带，7个颜色锚点依次为`#1D4C68,#4E8EAB,#A9CCD8,#F7F6F2,#E7C9AF,#C28561,#9D4C3D`。颜色在同一归一化下插值，三组共同线性色标仍为±7°C，中心仍为0；没有gamma增强、非线性归一化或单独调整组b。不能将配色变化误认为数据数值变化，也不宣称自定义色带通过完整色觉无障碍认证。

标题简化为a Non-recovery、b Recovery、c Larger bias，n分别为2、2、4；前两组的better-simulated stratification与第三组的larger stratification bias完整定义保留在本文件表格，移出图内重复长标题。c的Larger bias专指原分组的初始层结偏差，不指所有气候变量的综合模拟技能。圆框淡灰、海岸细线，坐标只在首图标注，三图范围/方向完全相同；图例与色标并列在底部。初始为细灰虚线，末期为较深黑实线；细白描边宽度较旧版减小，没有人为平移或分离重合边缘，也没有为美观平滑海冰线。

更细的50m海岸线由新增`src/s02_sia/s02_arctic_region/n067_pattern_cartography.py`单独准备，保存到`output/s02_sia/s02_arctic_region/p22_grouped_temperature_ice/n067_pattern_cartography/`。GeoJSON坐标维序为longitude,latitude、单位度；455个LineString/MultiLineString要素覆盖所有经度，筛选规则仅为原要素北界到达44°N，作图再裁到45–90°N。它是制图资产，不是新的TAS/SIC数据或海冰掩膜；未联网下载，也没有改变1°科学场的分辨率。

重设计前4个图件/manifest和1个QC已移动到`archive/20261021155222_tas_ice_visual_redesign/`，旧代码/说明3文件及科学NPZ/metadata/membership3文件已复制备份。共11文件均核对源目标状态、数量和SHA256；当前科学结果与这些快照逐字节相同，旧版可恢复。

## 2026-10-21可读性修订

按用户要求改为更简洁的论文式版面：去掉原大总标题；主标题缩短为Non-recovery、Recovery、Larger stratification bias；前两组的层结分类放在次级标题中，n仍明确为模式数。字体改为Arial，地图直径由原约3.48英寸增至约4.00英寸（约15%），经度标签移至圆框内侧，缩减色标刻度，减轻经纬网，略增强海岸与海冰线对比。初始仍为黑色实线，末期仍为绿色虚线，白色描边仍仅为可读性处理。没有新增插值、平滑或改变色标映射，±7°C共同色标不变，不能为了突出组b而单独增强其颜色。

旧图件及manifest共4文件、旧QC共1文件已移动到`archive/20261021154208_tas_ice_publication_layout/`，旧绘图/QC/说明和科学结果快照共6文件已复制保留。归档11文件均核对源/目标状态和SHA256；新QC还逐字节核对科学NPZ、metadata、membership与归档前一致。旧版可恢复，活动目录只保留新版图件。

首次版面检查发现色标说明靠近画布下沿，随后上移色标，为文字留足边距。该中间候选图件/manifest4文件和绘图代码副本1文件保存于`archive/20261021154503_tas_ice_colorbar_margin/`；5文件哈希、数量及移动后的源/目标均核对。最终PNG再次人工查看，色标说明完整，科学数组未重新计算。

## 可以怎样解读

共同色标下，组a北大西洋冷却与局地海冰边缘外移空间上对应；组b气温和两端海冰边缘变化较小；组c以增暖和海冰边缘向极侧退缩为主。不同海区可以出现不同方向变化，不概括为所有经度一致。空间对应关系只能用于描述和提出机制假设，不能仅凭叠加图证明气温变化由海冰变化导致，或反之。每组仅2、2、4个模式，需保留个体模式差异和既有分组性质的限制。

## 运行与核查

首次准备时依次运行`src/s02_sia/s02_arctic_region/n066_grouped_tas_ice_patterns.py`、`src/s02_sia/s02_arctic_region/n067_pattern_cartography.py`、`plot/S02_sia/s02_arctic_region/n066_grouped_tas_ice_patterns.py`、`src/s02_sia/s02_arctic_region/check_grouped_tas_ice_patterns.py`。本次仅运行新的制图资产准备、绘图及QC，未重跑已有科学计算。已有正式结果或图件时先按项目规则归档，不直接覆盖。代码全部顶层中文`# %%`顺序分块，无自编函数、类或lambda。

QC以独立矩阵权重复算三组均值，逐值核对原始成员温度和SIC，验证NaN未填0、陆地仅影响轮廓显示、所有来源哈希未变、6条初末轮廓均非空、输出尺寸/DPI/不透明背景、SVG可读文字和计算绘图分离。最终PNG另行人工查看。

scientific-visualization技能用于统一线性色标、初末海冰的灰度和线型冗余区分、保留缺测、避免用陆地遮挡气温，以及实际导出后的视觉检查；没有改变科学数值。
