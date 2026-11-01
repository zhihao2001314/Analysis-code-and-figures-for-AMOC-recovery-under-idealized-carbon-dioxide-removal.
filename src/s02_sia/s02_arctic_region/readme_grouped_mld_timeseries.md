# n061：三组混合层序列与观测参照

## 当前合图入口（2026-10-26）

用户要求与初始N²图按海区合并。当前活动图为n073两张2×3六联图，见readme_grouped_stratification_mld_composite.md；本目录模型序列原样使用，MLD观测仍来自n070 Argo hybrid。此前两张独立MLD图已归档；以下为原单图处理与图注历史，科学CSV、原始ELG历史结果和QC不变。

## 当前版：2026-10-24，观测改用Argo

两张活动图仍为n061原路径。只将ELG观测标记替换为n070保存的Argo作者密度混合算法（da_mld）结果；原模式CSV、分组、颜色、线型、11年平滑、纵轴与旧GFDL时序均不改。蓝色空心圆及垂直棒在三个面板左纵轴重复同一海区的均值±样本年际SD，不与实验年份对齐。

Argo沿用n069已核查的2°固定观测覆盖区、2005–2022年3月及9/18年格点支持门槛；年度片区面积覆盖≥50%才计入主参照。Labrador–Irminger16年：207.98±33.87m；Nordic13年：163.97±44.48m。不是原始剖面最大值，也不是完整海区气候态。具体平均、排除年份、阈值/混合算法敏感性、空间与算法不匹配限制均见`docs/grouped_march_mld_figure_explanations.md`当前版章节。

计算：`n070_grouped_argo_mld_reference.py`，输出n070独立目录，依次运行`check_grouped_argo_mld_reference.py`、n061绘图、`check_grouped_mld_timeseries.py`，导出检查使用配套Python的`check_grouped_argo_mld_exports.py`。正式重跑先归档。n061原始计算脚本及旧观测CSV/metadata保留为历史，不能再用它们描述当前Argo参照；绘图manifest指向n070。旧7个图件和1个QC已归档，4份修改前代码/说明已备份，清单在`archive/20261024104704_grouped_mld_replace_elg_argo/archive_manifest.json`。

Scientific-visualization技能影响观测显示和导出核查；软件参考记录在图外说明。此图保留模型/观测区域及MLD算法不一致、旧GFDL时序等局限，不宣称严格技能验证或期刊合规。

---

## 以下为ELG版本历史记录，不再是当前观测输入

**当前观测参照已改由n063独立观测密集片区提供，模式仍保留旧矩形；两者空间支持不同，不是严格技能验证。** 原n061 observation_reference.csv及下文旧矩形观测数值属于处理历史，当前绘图不再读取。当前ELG参照为西侧118.35±18.04m、Nordic121.83±36.28m；具体选择、门槛敏感性、范围图和局限见readme_observation_dense_mld.md及docs/grouped_march_mld_figure_explanations.md。

## 本次用户要求

延续n059八模式固定2/2/4分组，展示分组混合层序列。根据最新用户要求，取消原右侧窄栏，ELG观测均值和SD直接重复标注于每个分组面板的左侧纵轴，二阶导数法不再显示。仅保留三个模式组主面板，共用相同MLD纵轴；三个面板为同一海区的同一参照，不是独立样本，观测不是模式时间曲线，也不与实验年份对齐。

现用户明确改回9月14日旧八模式矩形序列；当前主图Labrador–Irminger，配套图Nordic Seas。沿用以下分组，不依据旧矩形N²重新筛选。三组成员：

- Better-simulated stratification / Non-recovery：CESM2、CNRM-ESM2-1。
- Better-simulated stratification / Recovery：NorESM2-LM、GFDL-ESM4。
- Larger stratification bias：ACCESS-ESM1-5、CanESM5、MIROC-ES2L、UKESM1-0-LL。

分组由用户指定，本次不重新检验深对流恢复，不把第三组统一称为恢复组，不覆盖其他大气或AMOC分组。

## 模式序列来源与处理

仅读取旧`output/s02sia/s02_arctic_region/p06_mld/n05_na_region_absolute_timeseries/na_march_mld_region_absolute_timeseries.csv`。八模式全部复用其原值和原已保存11年曲线，保持旧图数值，不再混用手绘沿岸MLD。Nordic矩形为65≤lat≤80°N、−30≤lon≤20°；Labrador–Irminger为50≤lat<65°N、−65≤lon≤−20°。

沿用原native areacello面积权重；旧实现以最接近10m的有效模式层sigma0为参考，增加0.03kg/m³后寻找首次向下crossing，不是新算法严格插值10m。无交点格点不作精确MLD混入平均。本轮不重算水柱、不取异常、不做漂移订正。**旧GFDL仍有已知索引/真实年份对齐局限，本轮为用户指定的历史曲线复现，未修正也未验证其时序推断。** 不用n055沿岸序列替换，因为掩膜不同；局限记录在图外说明，正式定量时序结论前需单独解决。

单模式使用旧表已保存11年居中平均，要求完整11年，绘图索引6–275，不填首尾。原280年索引表完整保留；GFDL首段是旧索引值，不代表找回了真实物理年1–5。组均值先对固定2/2/4成员逐年等权平均，再作相同平滑；不改变成员数。黑短虚线为组均值，不是观测。

## 历史：旧矩形观测纵轴标记的数值定义

遵循气候态加异常口径，src读取n060已完成的`march_reconstructed_grid.npz`，重新在上述两个旧矩形内面积平均，不能复用原手绘多边形观测均值。为核对完整气候态区域值，src额外读ELG和second_derivative三月clim MAT；plot只读取正式output。

```text
每格点重构值 = 同方法3月clim + 该年3月异常
每年区域值 = 当年有效重构格点的cos(latitude)面积加权均值
标记中心 = 2005–2025共21个三月区域值的等年权平均
误差棒上下端 = 中心 ± 21年样本标准差（ddof=1）
```

计算输出保留ELG与二阶导数法，不把两法平均成一个“真值”；当前按用户要求仅绘制ELG。本期本区域anchor100m与300m区域统计完全相同，因此每个面板显示100m结果一次。图例`Obs. ELG`对应蓝色空心圆；观测衍生重构参照的2005–2025时段写在图外，不称作原始剖面绝对MLD。

| 海区 | ELG均值 ± 年际SD (m) | 二阶导数均值 ± 年际SD (m) |
| --- | ---: | ---: |
| Labrador–Irminger | 148.86 ± 20.10 | 111.72 ± 14.79 |
| Nordic Seas | 122.52 ± 14.50 | 95.36 ± 8.68 |

**SD仅描述已处理区域序列的年际变率，不是测量误差、独立样本标准误或气候态置信区间。** 观测已上游平滑，不能再把21年当作21个独立样本作普通显著性检验。

现代三月旧框覆盖约为Nordic28%、Labrador–Irminger44%，负重构格点不截零；方法与空间支持不同于模式。虽然地理边界已一致，支持并未完全匹配。观测标记不能直接作为严格模式偏差检验。本图不改成完整区域clim均值（统计支持不同）。换回旧框后依然浅，不能将浅值完全归因于先前手绘范围。详见n060说明。

## 实际维度、文件与绘图

本轮不读原始模式NetCDF。观测clim MAT实际维序longitude,latitude，(360,180)，坐标0.5–359.5°E与−89.5–89.5°N，单位m；`loadmat(squeeze_me=True)`后一维坐标不改变场维序。n060 NPZ维序product4,year66,longitude360,latitude180。现有与新增CSV：

- 旧源表：model × simulation_year，两个区域以li/nordic前缀字段保存，单位m。
- `member_timeseries.csv`：8模式×2海区×280年，4480行，含原值与11年平滑。
- `group_mean_timeseries.csv`：3组×2海区×280年，1680行，含固定成员数量和当年有效数量。
- `observation_reference.csv`：2方法×2海区，4行，含均值、SD、覆盖率、负值比例及完整clim对照值。
- `observation_box_march_annual.csv`：2方法×2旧框×66个三月，264行。
- `membership.csv`、`metadata.json`、`qc_report.json`：分组、输入SHA256、物理口径与QC。

计算目录：`output/s02_sia/s02_arctic_region/p21_grouped_mld/n061_grouped_mld_timeseries/`。

绘图仅读取上述结果，不执行平滑、组均值或SD计算。三个主面板使用共同线性y范围：Labrador0–700m、Nordic0–600m，完整包含所有绘制曲线。各面板观测标记用`ax.get_yaxis_transform()`设置x轴坐标0、y数据坐标，故标记位于各自左纵轴对应深度，不属于实验第1年。第140.5年虚线分开名义CO₂增加与CDR阶段；旧GFDL的真实时序仍有局限。

图件目录：`figure/s02_sia/s02_arctic_region/p21_grouped_mld/n061_grouped_mld_timeseries/`。

- `labrador_irminger_grouped_march_mld.{png,pdf,svg}`
- `nordic_iceland_grouped_march_mld.{png,pdf,svg}`
- `figure_manifest.json`

尺寸13.8×5.8英寸；PNG6900×2900、500dpi、不透明白底；PDF/SVG为矢量，SVG文字可编辑。通用研究/PPT版式，没有指定或认证特定期刊要求。Scientific-visualization技能用于共同尺度、颜色/线型冗余、观测独立侧栏和导出后目视检查；软件参考沿用`readme_grouped_initial_n2_profiles.md`所列Scientific Agent Skills文献。

## 运行与QC

```powershell
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s02_sia/s02_arctic_region/n061_grouped_mld_timeseries.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'plot/S02_sia/s02_arctic_region/n061_grouped_mld_timeseries.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s02_sia/s02_arctic_region/check_grouped_mld_timeseries.py'
```

QC逐值核对旧源原值及原平滑曲线，用卷积复核平滑、权重矩阵复核组均值，从网格独立复算旧框观测区域值、均值和样本SD，检查来源哈希、行数、无自编函数及计算绘图分离、PNG尺寸/分辨率/不透明度和矢量文字。QC通过仅说明历史曲线复现、重构算术和显示通过，**不代表旧GFDL时序已验证**。两幅最终PNG均已查看，无标签裁切；SVG核查图内无方法脚注及阶段文字。

第一版7个图件/manifest的图例间距试排已归档至`archive/20261021110615_n061_legend_spacing/`原相对路径；核对源不存在、目标存在、数量7及各文件SHA256一致。现图件目录仅保留最终活动版；此前科学数据、原N²图、n060核查结果均不覆盖。

改回旧矩形前，n061手绘版全部6个计算/QC文件与7个图件/manifest已移动到`archive/20261021114414_grouped_mld_return_legacy_boxes/`原相对路径；源/目标、数量、SHA256全部核对，四个代码和说明文件也已复制保存并核对哈希。旧n05矩形科学源表保持不动，n059/N²及n060数据均不改。

## 无图内说明版

根据用户要求，两张图移除底部两行方法/限制脚注及主面板的阶段文字，只保留标题、坐标、图例、观测标记和读数。灰色竖虚线仍保留，其含义在图外解释。每张图的完整解释单独保存为`docs/grouped_march_mld_figure_explanations.md`，包含分组、区域、线型、时间处理、观测方法、SD和局限。科学CSV及metadata未修改；只重新导出图件和QC。

此前7个图件/manifest及1个QC文件已移动到`archive/20261021120643_mld_external_captions/`，逐一核对源路径不存在、目标存在、文件数量及8个SHA256；绘图代码、QC代码、本文档另复制备份并核对哈希。科学可视化技能在本次仅影响说明位置、底部留白和导出检查，不影响数据值或分组。

## 版式记录：观测先放在a面板左侧纵轴

右侧独立观测栏已取消，两个观测均值和垂直SD棒置于a面板左纵轴，上方保留方法符号图例，不在图内写观测数值、时段或方法说明。曲线、分组、深度范围、观测值、SD与平滑均不改。metadata只更新显示位置说明，具体当前显示以figure_manifest.json为准；计算代码中同步更新该显示字段，不重跑科学计算。

修改前7个图件/manifest及1个QC文件移动至`archive/20261021121918_mld_observation_left_axis/`原相对路径，8个文件的数量、源目标状态与SHA256全部核对。另保存绘图、QC、计算代码、README、图外说明及metadata共6个备份并核对哈希。所有科学CSV应与备份前逐字节一致；QC额外核对仅三个轴面板、SVG包含两种观测标记且无第四面板。科学可视化技能用于混合坐标放置、精简布局和最终图件目视检查。

## 版式记录：每个分组面板均放两种观测

用户进一步要求每个分组小图都有观测，故Labrador与Nordic两图的a、b、c面板左纵轴均重复同一海区的ELG和二阶导数法均值及SD。不是分组分别计算的观测，不改变统计样本量；同一海区三处中心和误差棒完全一致。科学数据、分组、纵轴和平滑均不改，计算代码和metadata不修改。

此前7个图件/manifest和1个QC已移至`archive/20261021122357_mld_observation_each_panel/`，源目标状态、7+1数量和8个SHA256已核对。另复制备份绘图、QC、README、图外说明和5个科学CSV共9个文件并核对哈希。QC核对每张SVG三组面板均有两种观测标记、同一方法的标记深度一致，以及全部既有科学算术检查；导出后逐图目视核查。科学可视化技能本轮仅影响观测重复显示和核验。

## 最新版：每个分组仅显示ELG

根据用户要求，Labrador和Nordic两图共六个面板移除二阶导数法标记及SD棒，上方对应图例一并移除。只保留ELG蓝色空心圆及垂直SD棒，数值与位置不变。二阶导数计算结果仍在正式CSV中，不物理删除；科学CSV、计算代码和metadata均不变。当前图外说明只描述实际绘制的ELG，其余方法仍作为处理历史记录。

旧7个图件/manifest和1个QC归档至`archive/20261021122755_mld_elg_only/`，源目标状态、7+1数量和8个SHA256全部核对，另备份绘图、QC、README、图外说明及5个科学CSV共9个文件并核对哈希。QC继续核对全部科学算术，并检查每个面板都有同深度的ELG标记、无二阶导数标记和图例。科学可视化技能仅影响显示选择与导出检查。

## 当前版：观测独立按覆盖定义范围

用户确认分别确定两个观测密集片区，不按用户分区算观测。本轮n063生成新参照，不覆盖n061旧科学CSV，n061绘图改读n063的observation_reference.csv；导出manifest明确记录新来源和模式/观测掩膜不匹配。n061计算代码仍可复现旧矩形产品，不重跑模式。新的检查先运行n063 QC，再核对新观测年度均值、SD和SVG中圆点实际深度位置；原矩形观测检查仅作为历史产品检查，不能误称当前显示同支持。

旧n061图件7项与QC1项归档至archive/20261021124214_mld_observation_dense_footprints/，8项哈希及源目标状态已核对，4个代码/说明和5个原CSV复制备份并核对。主片区门槛11/21年，敏感性8/14年均保存，不以MLD值挑选区域。完整方法及实际维度记录在同阶段readme_observation_dense_mld.md。
