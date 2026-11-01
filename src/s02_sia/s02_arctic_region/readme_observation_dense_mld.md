# n063：按观测覆盖确定ELG密集片区

## 用户要求与当前结果

分别为Labrador–Irminger附近与Nordic Seas附近确定观测密集区，不再按用户旧矩形或手绘分区算观测。模式曲线保持旧矩形不动。因此更新后的参照与模式空间支持不匹配，只能作为背景，不能严格计算模式偏差。

主结果（至少11/21个现代3月有效）为西侧210格点，118.346±18.044m；Nordic64格点，121.825±36.276m。±表示21个经过上游处理的区域年度值的样本SD，非置信区间或测量误差。现有产品仍是clim+经过时间平滑及残差中心化的异常重构，不等同于原始绝对MLD。

## 预先声明的工作规则

1. 时段固定2005–2025的21个3月；只显示ELG。以有值年数衡量处理后格点覆盖，不把它等同于原始独立剖面数量。
2. 气候态必须有限，ocean_mask=1、profile_count>0、first_guess_only=0；这仅说明该格点有气候态剖面支持，不保证每年都有直接剖面。
3. 主门槛为至少11年重构有限；另保存至少8和14年的敏感性。该门槛是本轮透明的工作选择，不是外部标准或最优性验证。
4. 为区分两海域，使用宽定位窗口：西侧45≤lat<65°N、70–10°W；Nordic65≤lat<85°N、40°W–30°E。窗口由分析者指定，不声称完全没有人为地理先验；不是用户旧框，不是最终平均范围。
5. 在各定位窗口取满足条件的最大8邻域连通片（以格点数最大判定，相同时NumPy标签顺序优先）；对角相邻也相连。不填缺测孔洞，不扩张掩膜，不按MLD深浅选择，不按模型结果选择。
6. 对每年的固定主片区内当年有限重构值用cos(latitude)面积加权；NaN不作0，负重构值保留。再对21个年均值等权平均，SD用ddof=1。因此仍会受年间空间支持变化影响。

## 实际范围与覆盖

- 西侧主片区中心坐标外包范围45.5–63.5°N、59.5–26.5°W，210格点；**超出严格Labrador–Irminger海域，不能称局地深对流核心区。** 片区内各年可用面积平均61.88%，最低44.01%、最高73.56%。
- Nordic主片区中心坐标外包范围65.5–74.5°N、0.5–20.5°E，64格点；片区内各年可用面积平均63.16%，最低27.78%、最高83.71%。

以上坐标仅是已选离散格点的外包范围，不是把整个矩形都纳入。完整格点以selected_cells.csv和coverage_and_footprints.npz为准。覆盖率分母为各自固定密集片区面积，不能与原矩形44%/28%直接当成同分母改善幅度。

## 门槛敏感性

| 至少有效年数（共21年） | 西侧格点数 | 西侧均值±SD (m) | Nordic格点数 | Nordic均值±SD (m) |
| --- | ---: | ---: | ---: | ---: |
| 8 | 603 | 166.01 ± 19.20 | 114 | 121.11 ± 26.64 |
| 11（当前展示） | 210 | 118.35 ± 18.04 | 64 | 121.83 ± 36.28 |
| 14 | 31 | 80.27 ± 17.76 | 21 | 91.03 ± 59.90 |

门槛同时改变片区与可用支持。敏感性明显，不能只报道11年而隐去其他结果，更不能选择让观测看起来更深/更接近模式的门槛。改变观测区域并未解决上游算法、方法定义及绝对MLD有效性问题。

## 来源与实际维度

- 输入n060的march_reconstructed_grid.npz：product4,year66,longitude360,latitude180；product0为ELG_100m，year45:66为2005–2025年3月。场单位m。
- ELG/clim/MLD_clim_month_03.mat：mld_obs、ocean_mask、profile_count、first_guess_only均为longitude360,latitude180。坐标longitude为0.5–359.5°E、latitude为−89.5–89.5°N；不转置解释场。
- 新masks为region2,longitude360,latitude180；valid_year_count为longitude360,latitude180，整数0–21。
- 地图底图来自本机已缓存Natural Earth110m land，由计算脚本转为正式output内map_land.geojson；绘图只读output，不联网或读原始观测。
- 不读取用户分区几何以确定掩膜；旧观测CSV仅读取作数值对照，不参与筛选。全部来源SHA256写入metadata。

## 输出、图件与说明

output/s02_sia/s02_arctic_region/p21_grouped_mld/n063_observation_dense_mld/包含主参照2行、敏感性6行、逐年126行、主格点274行、掩膜NPZ、地图底图、元数据及QC。原n061的矩形观测CSV保留为历史产品，当前n061绘图只读取n063主参照。n061所有模式CSV及分组不改。

地图仅显示有效重构年数与橙色主片区轮廓，说明在图外。序列图蓝色空心圆为新的ELG参照，各分组重复同一值；无二阶导数法。黑虚线仍为模式等权均值、彩色为模式、11年居中平滑和名义140.5年阶段线不变。完整逐图解释见docs/grouped_march_mld_figure_explanations.md。

科学可视化技能影响范围图、共同尺度、说明位置和导出核对，不据此宣称观测物理有效或期刊合规。软件技能参考沿用n059/n061已有软件说明。

## 运行顺序

```powershell
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s02_sia/s02_arctic_region/n063_observation_dense_mld.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'plot/S02_sia/s02_arctic_region/n063_observation_dense_mld.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s02_sia/s02_arctic_region/check_observation_dense_mld.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'plot/S02_sia/s02_arctic_region/n061_grouped_mld_timeseries.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s02_sia/s02_arctic_region/check_grouped_mld_timeseries.py'
```

已有正式输出或图件须先按archive/README归档。旧n061图件7项和QC1项已移至archive/20261021124214_mld_observation_dense_footprints/原相对路径，数量、源目标状态和8项SHA256全部核对，另保存4个代码/说明及5个原科学CSV备份。新旧观测取样方案均有明确来源，不覆盖原数据。
