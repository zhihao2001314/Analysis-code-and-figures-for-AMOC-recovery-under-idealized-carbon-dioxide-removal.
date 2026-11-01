# 直接mlotst与thetao/so诊断MLD的配对比较

2026-10-28。按作者要求，检查当前温盐MLD空间产品与模式直接输出mlotst的差异，不覆盖原图、不替换论文MLD来源。

## 范围与输入

包含7模式：ACCESS-ESM1-5、CanESM5、CESM2、CNRM-ESM2-1、GFDL-ESM4、NorESM2-LM、UKESM1-0-LL。当前本地MIROC没有强迫mlotst文件；CAS仅有升高阶段mlotst，缺少本轮需要的CDR文件及既有同口径温盐端点产品，因此两者不配对。

温盐端点直接读取`output/s02sia/s02_arctic_region/p06_mld/n02_na_deep_mixing_three_phase/<model>_na_march_mld_three_phase_native.npz`。它们由TEOS-10 sigma0、插值到10m的参考密度及增密0.03kg/m³阈值计算，无交点使用最深有效采样层；不是另一个最近10m参考、无交点排除的旧区域时序算法。

mlotst从原始Omon NetCDF读取，文件清单来自既有输入QA，并再次检查成员、试验、维度、单位和实际年月。每个模式各自使用与温盐产品完全相同的初期和末期10个March；GFDL初期明确为第6–15年，不能与mlotst第1–10年混比。其余模式的日期由原始温盐metadata读取，不假设所有模式日历起点相同。

本比较限于初末十年气候态及late−initial变化；不声称完成了独立旧版280年区域时序的算法或时间轴验证。

## 配对与统计

1. 对同成员、同试验、同年月的mlotst取10个March等权平均；要求两时期各10个值均完整，不跳过缺年。
2. 对照温盐产品的原生格点逐点查找坐标对应关系，本轮所有配对距离均为0m。没有插值或重网格。
3. 两种方法、两个时期使用同一有效支持；权重为同网格原生areacello，而非cos(latitude)。area原始格点也逐点核验坐标。
4. 区域为LI：50≤lat<65、−65≤lon<−20；NS：65≤lat<80、−30≤lon<20。另保存50≤lat<85、−65≤lon<20的完整温盐诊断域，作为背景统计，勿与海冰NA扇区混称。
5. 误差符号统一为温盐诊断减mlotst。逐区域保存初期、末期和变化量的面积加权均值偏差、MAE、RMSE、空间相关，以及50/100m误差分箱的面积比例。分箱仅作描述，不是显著性或验收阈值。
6. 空间变化符号一致率另保存仅至少一种方法|ΔMLD|≥50m的描述性子集；主结论13/14指14个区域均值变化的符号一致率，不是逐格显著性检验。

## 已确认的实际维序

|模式|mlotst|areacello|经纬度变量|
|---|---|---|---|
|ACCESS-ESM1-5|time,j,i|j,i|latitude,longitude|
|CanESM5|time,j,i|j,i|latitude,longitude|
|CESM2|time,nlat,nlon|nlat,nlon|lat,lon|
|CNRM-ESM2-1|time,y,x|y,x|lat,lon|
|GFDL-ESM4|time,y,x|yh,xh|lat,lon|
|NorESM2-LM|time,j,i|j,i|latitude,longitude|
|UKESM1-0-LL|time,j,i|j,i|latitude,longitude|

温盐NPZ：latitude、longitude及各MLD场均为同形二维原生裁剪网格。正式维度记录在grid_inventory.csv，原始mlotst的月份和文件记录在matched_months.csv与mlotst_input_inventory.csv。配对表为model×native_cell，共34906行；月份为7模式×2时期×10个March，共140行；统计为7×3区域×3变量，共63行。

## 结果

两种MLD的多数区域变化方向与空间格局接近，但不能概括为所有模式、所有海区均差异可忽略：

- 14个核心模式—海区组合中，13个区域均值变化同号。
- 区域平均变化量之差的绝对值中位数为8.8m，最大64.8m；8/14小于等于10m，10/14小于等于30m，13/14小于等于50m。这些是本次样本描述，不是预设通过标准。
- 变化场的面积加权空间相关为0.786–0.998，中位数0.973。
- CanESM5北欧海初期均值由mlotst的257.8m变为温盐诊断369.7m；末期335.1m与382.2m，对应变化+77.3m与+12.4m，变化量差64.8m。
- CESM2北欧海变化由+4.2m变为−44.4m，符号发生变化；其mlotst区域平均响应本身很弱，不能描述为强烈变深与强烈变浅的对立。
- NorESM2-LM的LI末期均值分别为477.1m和376.9m，绝对差约100.2m；变化量分别−99.4m和−149.0m，均为变浅，但幅度差约49.6m。

|模式|LI变化：mlotst / 温盐诊断(m)|NS变化：mlotst / 温盐诊断(m)|
|---|---:|---:|
|ACCESS-ESM1-5|+62.0 / +64.3|+28.9 / +22.8|
|CanESM5|+88.6 / +57.0|+77.3 / +12.4|
|CESM2|−360.0 / −361.8|+4.2 / −44.4|
|CNRM-ESM2-1|−222.7 / −216.2|−258.8 / −250.1|
|GFDL-ESM4|+304.7 / +281.5|+24.5 / +30.6|
|NorESM2-LM|−99.4 / −149.0|+18.5 / +46.1|
|UKESM1-0-LL|−191.9 / −183.1|−253.5 / −251.8|

本次是两套诊断产品的一致性比较，不以原生mlotst作为绝对真值；也不将差异全部归因于某一个阈值、密度参考或时间平均步骤。论文仍应准确描述实际采用的温盐诊断方法。

## 代码、正式结果及图件

计算：`src/s02_sia/s02_arctic_region/n079_mld_method_comparison.py`

绘图：`plot/S02_sia/s02_arctic_region/n080_mld_method_comparison.py`

QC：`src/s02_sia/s02_arctic_region/check_mld_method_comparison.py`

结果：`output/s02_sia/s02_arctic_region/p24_mld_method_comparison/n079_mld_method_comparison/`

正式图件：

- `F:/sea/Antic sea ice constraint/figure/s02_sia/s02_arctic_region/p24_mld_method_comparison/n080_mld_method_comparison/native_mlotst_vs_temperature_salinity_mld.png`
- `F:/sea/Antic sea ice constraint/figure/s02_sia/s02_arctic_region/p24_mld_method_comparison/n080_mld_method_comparison/native_mlotst_vs_temperature_salinity_mld.pdf`
- `F:/sea/Antic sea ice constraint/figure/s02_sia/s02_arctic_region/p24_mld_method_comparison/n080_mld_method_comparison/native_mlotst_vs_temperature_salinity_mld.svg`

图注：Comparison of March mixed-layer depth from native CMIP6 mlotst and temperature–salinity diagnostics. Panels show regional mean initial MLD (a), late MLD (b), and late-minus-initial change (c) on identical native-grid support. Circles and squares indicate the Labrador–Irminger region and Nordic Seas, respectively; colours distinguish seven models. Dashed lines indicate equality. Each paired comparison uses the same ten March samples, including years6–15 for the GFDL initial window.

## QC及导出检查

独立使用np.average重算63条均值与误差统计，最大差1.14×10⁻¹³m。另从原始NetCDF独立重读7模式×2区域×2时期的28个深混合代表点，十年均值逐一核对通过。所有配对格点坐标距离为0m；在温盐产品保留的有效面积内，mlotst配对覆盖均为100%。

最终PNG为5580×2160、450dpi，PDF为矢量单页，SVG保留文字。已查看最终PNG和独立Poppler PDF预览，标签、两组图例及数据点完整，横纵轴范围相同，未显示总体标题。预览位于`archive/20261028_mld_comparison_label_review/review/`。该归档还保留了仅将坐标轴mld改为MLD之前的3份初次导出，源路径、目标路径、文件数及SHA256已核对。

图形规范沿用scientific-visualization技能；其软件引用沿用`readme_grouped_initial_n2_profiles.md`中的既有记录，不在本轮自动修改论文参考文献。
