# 两区九模式与三组地图的紧凑组合图

## 当前图件

- 计算：`src/s02_sia/s02_arctic_region/n093_compact_regional_hydrography.py`。
- 输出：`output/s02_sia/s02_arctic_region/p23_stratification_mld_composite/n093_compact_regional_hydrography/`。
- 绘图：`plot/S02_sia/s02_arctic_region/n094_compact_regional_hydrography.py`，只读取output。
- PNG：`F:/sea/Antic sea ice constraint/figure/s02_sia/s02_arctic_region/p23_stratification_mld_composite/n094_compact_regional_hydrography/li_nordic_nine_models_grouped_mld_ice.png`；同名PDF/SVG。

## 七面板布局

第一行a/b分别为LI和Nordic Seas初始三月N²，九模式合并显示，黑线和灰带为对应IAP均值与年际样本标准差。第二行c/d分别为两区九模式MLD时间序列。上述四幅图不含任何模式组平均线，模型颜色与40°N AMOC主图一致。
IAP仅筛选10–1000m以内已有中间层点，即15–995m；没有把5m值外推到10m或人为补1000m边界。模型剖面深度及数值不变。
N²横轴采用分段线性坐标，说明放在图注与元数据中，图内不再显示Piecewise-linear x-axis文字，物理刻度0、1、2、4、6、8等距显示。原值节点[-1.6,0,1,2,4,6,8,9.2]对应显示坐标[-0.22,0,1,2,3,4,5,5.6]，使用SciPy interp1d与Matplotlib function scale转换；曲线及灰带同步变换，未改原始N²值。负值域压缩到0.22个显示刻度间隔，零值更靠近左边界；负值灰带仍完整保留。两图原值范围仍为−1.6到9.2，以完整保留Nordic Seas在35m处约−1.416的IAP均值减1标准差阴影；该下界是年际离散带，不是把模式N²置负。

c/d纵轴统一为0–700m，每100m一刻度，所有11年平滑曲线完整保留。横轴物理刻度在最终画布上的像素间距已逐项核对相等。

第三行e/f/g为三组残余MLD地图：

| 组 | 成员 |
|---|---|
| Non-recovery | CESM2、CNRM-ESM2-1、UKESM1-0-LL |
| Recovery | NorESM2-LM、GFDL-ESM4、CAS-ESM2-0 |
| Larger bias | ACCESS-ESM1-5、CanESM5、MIROC-ES2L |

这是作者指定的新3/3/3分组，UK移入Non-recovery、CAS保留Recovery。原三模式Recovery的MLD地图与n091中列逐元素相同；另外两组重算。

## 空间与数据口径

区域平均均沿用矩形：LI50–65°N、65–20°W；NS65–80°N、30°W–20°E。N²采用半开边界，MLD沿用原存档矩形的外缘取格约定。
N²读取n085/n087现成九模式剖面，MLD读取n089统一矩形结果。全套1800个剖面记录和5040个MLD记录只整理、不重新计算或平滑，不引入piControl订正。GFDL仍以源6年为分析1年。
地图填色为每模式末十年减初十年MLD，再作固定三成员等权平均；黑线为同三成员初始MLD平均，线级与标签均为200、400、600、800、1000m，间隔统一200m。填色范围改为−1000到1000m，分级间隔50m，标注−1000、−500、0、500、1000m；超范围值用两端延伸色表示，原值完整保留。沿用原0.25°显示网格，所有成员同时有效才平均，不进行skipna或新增插值。

## 海冰边缘线

紫虚线/绿虚线分别是初期/末期组内等权平均年均SIC场的15%等值线。先平均浓度、再取等值线；没有平均边缘线坐标。
九模型SIC均为同一1°显示网格35×360，三成员共同有效才给出均值。旧八模型继承的稀疏表零值恢复保持原样，CAS已有NaN保留，没有新增补零、插值或变动成员数平均。
粉色观测冰缘为三图共用的NOAA/NSIDC1988–2020年平均场15%等值线，读取旧EPSG3411网格448×304及已有显示平滑，不新增平滑。
海冰年均边缘与三月MLD的月份不同，图注中分别说明。

## 维度与验证

不新增NetCDF/MAT读取。剖面长表region2×model9×depth100，IAP筛选后region2×depth99；MLD长表region2×model9×year280；地图各组141×341，组内海冰各时期35×360。
38条实际绘制曲线逐条哈希核对；N²均值及完整±1SD阴影未被坐标范围截断。三个组的MLD和SIC固定成员平均已独立重算核对，六条模式组冰缘均有有效等值线段，观测场逐元素不变。
共用图例与坐标文字自动检查无重叠，PNG和独立PDF栅格均已目视检查。PNG5130×4770、450dpi，PDF单页且无Type3字体。
旧n092九面板组合图已归档，目前活动组合图为n094；单模式诊断图不在本次修改范围。
