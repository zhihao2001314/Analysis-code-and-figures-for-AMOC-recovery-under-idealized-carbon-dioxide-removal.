# CAS温盐诊断混合层和层结

2026-10-29按作者要求新增CAS-ESM2-0诊断。仅计算CAS，八模式、IAP观测和既有分组结果不改动。计算入口`n083_cas_mld_stratification.py`；绘图入口`plot/S02_sia/s02_arctic_region/n084_cas_mld_stratification.py`仅读output。

## 原始维度和时间

thetao、so均为标准原生月文件，实际维序`(time,lev,lat,lon)`；30×196×360，深度单位m、向下为正，最浅中心5m、下一层15m。thetao单位degC，so单位0.001，按潜温及Practical Salinity解释。areacello为`(lat,lon)`、m²。两实验、两变量、面积的每个经纬度及温盐深度逐值相同。输入空间截取50–85N、65W–25E，保留整个既有海岸分区；最终地图展示仍为65W–20E、50–85N。

只读r1i1p1f1/gn，严格文件名筛选排除yearmean及r360x180派生文件。升CO₂来自3个文件，200–339年140个March；降CO₂来自4个文件，1–140年140个March。完整月轴逐文件连续、目标March无重复无缺失，温盐真实year-month逐一匹配。保存`input_inventory.csv`及`matched_marches.csv`。两段统一实验年1–140、141–280，不把CDR稳定阶段包括进来，不由文件尾部选窗。

## MLD两种既有处理口径

共同TEOS-10链：几何深度及纬度转压力，SP→SA，潜温→CT，再计算sigma0。当地密度首次达到参考值+0.03 kg/m³的位置为MLD，在相邻层之间线性插值。

1. 正式空间图和海岸分区：sigma0线性插值到恰好10m。`mld_crossing_m`只保留真实交点；`mld_all_m`将有有效10m参考、但直到最深有效层仍不达阈值的柱记为最深采样深度，并以no_cross标记。这是采样/海底受限值，不解释为已解析的深混合层。
2. 矩形旧序列的直接对照：使用最近10m原生层，CAS的5m/15m等距时选5m，且无交点记缺测。单独保存在`mld_nearest_crossing_m`，不混入正式地图或10m参考结果。

每March先逐水柱诊断，再对区域有效交点用areacello权重平均。统计使用未平滑年值；曲线显示完整11年居中平均，首尾各5年保持缺测。地图按四个十年窗逐格10/10有效年取平均，再取四时期共同有效支持：up_initial1–10、up_final131–140、down_initial141–150、down_final271–280。

## N²计算

每March每水柱先将温盐线性插值到0–1000m、10m间隔，无垂向外推，再转SA/CT并用`gsw.Nsquared`。N²中点为5–995m；0m无温盐外推，因此5m层N²保留NaN，15–995m有效。原生温盐有效层检查从表面向下连续，仅允许海底以下缺测，符合本CAS固定z坐标；没有人为填补内部缺层。无额外垂向平滑或piControl漂移订正，不截负N²。

同一March在各深度用原生areacello加权，再取每个十年窗的区域剖面均值/年际样本SD；需要十年均有效，不以nanmean改变时间样本数。逐年区域剖面及有效面积比例也保存，未绘制成未定义的垂向积分指标。

## 分区与比较边界

- coastal：完整沿用作者手绘、经海岸裁剪的两个既有GeoJSON多边形，CAS N²及MLD都在相同多边形计算。
- rectangles：LI=50≤lat<65、65W–20W；NS=65–80N、30W–20E。为对应旧结果，MLD保留旧包含外边界的选点，N²使用原半开边界；两种掩膜都保存在数组中。

既有八模式N²及IAP参照直接读取各自同版分区结果，未重算。八模式MLD也直接读取已存绝对序列；旧GFDL索引/真实年对齐局限不因叠加CAS而消失。本CAS的280个March全部使用已核对真实标签。IAP为2005–2024年March现代气候参照；初始十年模式与IAP不是同期观测。八模式GFDL初始N²仍为原批准的6–15年例外。本次不新增CAS区域技能分类或改变2/2/4组。

## 输出与图注

计算结果：`output/s02_sia/s02_arctic_region/p25_cas_hydrography/n083_cas_mld_stratification/`。

- `cas_native_mld.npz`：MLD实际维序year280×latitude36×longitude91；四时期MLD及三个变化场。区域底部限制标记单独保存。
- `cas_regional_n2.npz`：逐年N²为year280×region4×depth100，原生四时期N²为period4×depth100×latitude36×longitude91。
- `cas_march_mld_timeseries.csv`、`cas_march_n2_profiles.csv`、`cas_endpoint_summary.csv`。
- `metadata.json`、原文件清单、逐月配对、独立原柱样本和QC。

正式图件：`figure/s02_sia/s02_arctic_region/p25_cas_hydrography/n084_cas_mld_stratification/`。

`cas_coastal_n2_mld_comparison`及`cas_rectangular_n2_mld_comparison`：a/b为两海区初始/末期CAS N²（橙实线/虚线），黑线IAP、浅灰原八模式初始N²。c/d为CAS逐年March MLD（浅橙）及完整11年居中均值（深橙），浅灰为原八模式相应平滑MLD，竖点线为CO₂反转。coastal下排用10m参考真实交点，rectangular下排用最近10m层真实交点；两者不能与all/bottom-limited混称。

`cas_march_mld_changes`：a升高期末减初、b降低期末减初、c全实验末减初。蓝为变浅，红为加深，三幅共用±1500m色标，两端箭头表示超出色限。黑线分别为该变化所对应初期的MLD，100、200、500、1000、1500、2000m；地图不显示统计显著性结论。不加总标题/时期副标题，含义放在本说明/后续图注。

图件均输出PNG/PDF/SVG，源码路径区列完整文件路径，运行后逐个打印。底图只从output中已备好的GeoJSON读取，没有绘图时下载底图或读取原始温盐。

## 独立验证

`check_cas_mld_stratification.py`独立从原NetCDF读取四个阶段各12个水柱（48柱），逐层检查原温盐，使用标量阈值/交点方法核对MLD，并用TEOS-10离散压力公式、`gsw.specvol_alpha_beta`复算N²，而非再次调用Nsquared。原始月配对280个、旧产品源SHA256、原文件大小及mtime、区域面积权重、11年平滑缺测边缘、N²5m缺测均通过。

MLD最大误差4.62×10⁻¹¹m，N²最大误差1.09×10⁻¹⁹s⁻²；区域MLD矩阵点积复核最大误差2.85×10⁻¹³m。科学绘图技能用于共同色标、初末/模式/观测可辨识、缺测和无交点限制的如实显示，以及实际导出复核；不宣称特定期刊合规。

## 九模式总览入口

在保留本页CAS单独诊断图的基础上，九个单模式分别与IAP对照、九模式MLD空间变化及九模式海岸分区MLD序列由n085准备、n086绘图。具体口径和三个正式图名见`readme_nine_model_hydrography_panels.md`；图件在`figure/s02_sia/s02_arctic_region/p25_cas_hydrography/n086_nine_model_hydrography_panels/`。
