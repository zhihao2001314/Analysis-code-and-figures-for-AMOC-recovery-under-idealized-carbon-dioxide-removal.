# 六模式净降水与大气环流

## 口径与文件

- 计算脚本：`src/s05_atmosphere/p05_amoc_recovery_composites/n81_six_model_net_precipitation.py`。
- 计算产品：`output/s05_atmosphere/p05_amoc_recovery_composites/n81_six_model_net_precipitation/`。
- 绘图脚本：`plot/s05_atmosphere/p05_amoc_recovery_composites/n82_six_model_circulation_net_precipitation.py`。
- 最终 PNG/PDF/SVG：`figure/s05_atmosphere/p05_amoc_recovery_composites/n82_six_model_circulation_net_precipitation/annual_circulation_global_net_precipitation_six_models.*`。
- 原图4降水版本保留。本图为净降水这一不同变量的正式结果，没有覆盖原降水图。

恢复组：GFDL-ESM4、NorESM2-LM、CAS-ESM2-0；未恢复组：CESM2、CNRM-ESM2-1、UKESM1-0-LL。
每组固定三模式等权；差值为未恢复组减恢复组。初末十年为实验年1–10和271–280。

净降水为 P−E；E 使用 evspsbl，向上为正，包含地表蒸发、升华和蒸腾的模式总通量。
原始单位 kg m⁻² s⁻¹ 乘86400转换为mm day⁻¹。净降水变化为 Δ(P−E)=ΔP−ΔE，
正异常只表示相对于初态更偏向地表淡水获得，不等于绝对P−E必然为正。
没有采用潜热通量替代蒸发，没有进行piControl漂移订正。

## 原始维序、插值和时间

实际原始evspsbl维序为(time,lat,lon)，所有文件已查验属性并保存清单。
CNRM升高实验standard_name为water_evaporation_flux，降低实验为water_evapotranspiration_flux，
二者说明均为地表向大气的水汽通量；其余五模式为water_evapotranspiration_flux。
CESM2文件有两种数值类型的1e20填充值属性；按Xarray标准解码后，所选120个月均无缺测。

每个时期完整120个月，按真实日历逐月天数加权成年平均，再对原生年平均场作周期经度双线性插值。
目标沿用已有全球降水网格：纬度−88:2:88共89点，经度−180:2:178共180点，不外推。
十年平均采用实际年度天数加权。降水直接读取n72已核验的逐年结果和日权重，不重算降水。

| 模式 | 初始源年 | 末期源年 |
|---|---|---|
| GFDL-ESM4 | 1–10 | 271–280 |
| NorESM2-LM | 1–10 | 271–280 |
| CAS-ESM2-0 | 200–209 | 470–479 |
| CESM2 | 1–10 | 131–140（CDR） |
| CNRM-ESM2-1 | 1850–1859 | 2120–2129 |
| UKESM1-0-LL | 1850–1859 | 2120–2129 |

保存样本维序：(model6,period2,year_in_window10,lat89,lon180)。
模式气候态为(model,period,lat,lon)，模式变化为(model,lat,lon)；
组平均气候态为(product2,period,lat,lon)，组响应及显示掩膜为(product3,lat,lon)。
product顺序为recovered、not_recovered、contrast。

## 检验与图件

组内两个P−E面板采用三模式变化严格同号判据，零值不通过；该判据不属于p值显著性检验。
组差重新对P−E计算中心化双侧循环年块bootstrap，9999次，块长2/3/4年，取最大p后在该全球场上作BH校正，q<0.05。
复用n73保存的随机年索引，但没有复用降水p、q或掩膜；每模式每时期独立抽样，空间共享抽样索引。
推断针对六个固定模式的时间抽样不确定性，不推断总体模式抽样不确定性或因果关系。

左列环流场、初始等值线和组差显著性完全沿用n74，未重新计算。
右列换成P−E。按用户后续要求，b/d叠加组平均绝对净降水初期/末期峰值圈叉，
从20°S–20°N范围内逐经度取最大P−E所在的原始2°网格纬度，再每10°经度显示一次。
蓝色空心圆为Initial，品红叉为Later，同经度有位移时加细灰连线；f不叠加峰值。
峰值取自初末期绝对P−E场，不是Δ(P−E)异常场，也不复用原降水最大值。
该指标供观察净降水雨带位置，与降水最大值或降水中位数ITCZ指标不等同。
计算脚本n83_net_precipitation_peak_markers.py仅读取n81已有气候态；保存全180经度峰值、
峰值大小、是否恰在搜索边界、并列最大值和36个显示点的选择；全部保留，不删掉边界峰。
本次没有重算科学场或bootstrap。旧无峰值图件已按原项目相对路径归档到
`archive/20261031161804_net_precipitation_peak_markers/`，逐文件SHA256与源路径移出已核验。
沿用原图发散配色及±3 mm day⁻¹范围；所有P−E组响应值均在色标范围内。
绿色为净降水增加，棕色为减少。点状标记通过对应判据，网纹未通过；点状标记稀疏采样以保证可读性。

独立检查脚本 `check_six_model_net_precipitation.py` 核验P−E恒等式、分组、六模式原始首年蒸发格点、
27个直接年度索引bootstrap结果、BH校正、图件数据哈希及环流不变。
数据清单、时段QC、来源哈希、统计诊断和最终图注均随输出保存。
