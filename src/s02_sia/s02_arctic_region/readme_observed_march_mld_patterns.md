# n065：SEANOE与ELG重构的3月MLD pattern

## 目的和固定口径

用户要求先查看刚下载的SEANOE v2023空间分布，再查看现有重构观测分布。两张独立图均为3月，地图目标范围75°W–30°E、40–85°N，相同Lambert投影、相同线性0–1400m色标，不按旧矩形或手绘区域裁取数据，也不套用n063密集区。仅作独立产品空间诊断，未替换既有混合层序列或观测参照。

此任务使用scientific-visualization技能检查相同尺度、缺测与覆盖标记、来源和实际导出；未指定期刊，不声称期刊合规。技能软件参考沿用n059/n061已有说明。

## 图1：SEANOE v2023

文件：figure/s02_sia/s02_arctic_region/p21_grouped_mld/n065_observed_march_mld_patterns/seanoe_v2023_march_mld_pattern.png。

输入为用户下载并解压到F:/sea/mldobs/SEANOE-obs/mld_dr003_ref10m_v2023.nc。实际MLD变量mld_dr003维序time12,lat180,lon360，单位属性为unit=m（不是units）；n_obs同维序，mask为lat180,lon360。lon为−179.5到179.5°E，lat为−89.5到89.5°N。time不是逐年序列，而是12个月的月中心日序[14,45,73,104,134,165,195,226,257,287,318,348]；索引2对应3月。FillValue=−99999，由Xarray解码为NaN，不作0。

变量元数据定义MLD为相对于10m位势密度增加0.03kg/m³的阈值深度。资料为1970-01-01至2021-12-31实测剖面构建的月气候态，经OI/kriging空间映射，来源包括Argo、WOD内多种观测，不是仅Argo，也不是未映射原始观测。DOI：10.17882/91774。

填色为MLD，单位m；越蓝越深，浅黄越浅。灰色为陆地，淡蓝灰背景为缺测。黑点表示该1°三月格点n_obs=0，但映射气候态仍有限，即无该格点本地剖面支持，不是MLD=0或显著性检验。这样的映射值保留，不用n_obs筛掉；不将n_obs当面积权重。

## 图2：ELG重构观测

文件：figure/s02_sia/s02_arctic_region/p21_grouped_mld/n065_observed_march_mld_patterns/ELG_reconstructed_march_mld_pattern.png。

使用既有n060网格中的ELG_100m方案，逐格点取2005–2025年21个3月中可用值的等年权平均。每年的值为ELG三月clim+已平滑并再次月中心化的异常。不把NaN当0，不截负重构值，不取深度门槛，不按覆盖率删除任何有限格点。当前地图范围内多年平均场没有负值或零值；部分单年重构可以为负，其已保留在计算中。

输入NPZ维序product4,year66,longitude360,latitude180；选product0、year45:66。气候态MAT维序longitude360,latitude180，坐标0.5–359.5°E与−89.5–89.5°N。明确按经度排序并转置为lat180,lon360，与SEANOE坐标逐点一致，不是靠维数顺序猜测。

填色与图1使用相同单位、范围和色标。黑点表示21个3月中只有1–10年有效；这些格点仍绘出，缺测年份不补0。无点只能说明至少11年处理后格点有值，不保证有11个独立直接剖面。黑点含义与图1不同，不可把两图的点数直接当成同一覆盖指标。灰陆地、淡背景和无显著性检验同图1。

## 当前结果与解释边界

| 图示范围内诊断 | SEANOE | ELG重构 |
| --- | ---: | ---: |
| 有限格点数 | 3468 | 2193 |
| 最大值 | 1221.18m | 350.90m |
| 超过500m的格点数 | 89 | 0 |
| 超过1000m的格点数 | 9 | 0 |

SEANOE最大值位置57.5°N、54.5°W在Labrador海；对应ELG格点约79.72m，SEANOE该格点24个三月剖面，ELG有16个可用三月。此点按SEANOE图内最大值选取，仅用于检查相同坐标定位，不是海区平均，也不能推广为全区域误差。

在两者均有限的2193个共同格点上，面积加权均值为SEANOE192.83m、ELG133.44m；该数值仅供内部诊断，不能称偏差，因资料期、MLD定义、映射和抽样均不同。网格坐标一致也不等于观测空间支持和取样历史一致。

ELG重构多年平均与其输入3月气候态在图示范围内最大绝对差约1.53×10^-5m，符合上游按2005–2025逐月残差中心化的构造。因此重构平均图主要反映原气候态，浅值不是本轮区域平均分母造成的。图示差异提示需回查ELG上游温盐输入、MLD算法/参数、质控和映射，但这两幅图本身不证明具体哪一步有错。不能通过改变分母或只选深值来修正产品。

## 结果与复现

计算src/.../n065_observed_march_mld_patterns.py，绘图plot/S02_sia/s02_arctic_region/n065_observed_march_mld_patterns.py，独立检查check_observed_march_mld_patterns.py。output同阶段n065目录保存全域两产品地图数组、各自覆盖信息、原ELG气候态、诊断表、地图底图、完整源属性和SHA256。底图由n063已保存Natural Earth110m land复制，绘图仅读正式output，不联网。

PNG为4050×3780、450dpi白底；PDF/SVG为矢量导出，SVG文字可编辑。图内只保留标题、时期、坐标、MLD色标、覆盖符号图例和Greenland地名；所有详细解释在本文档。相同0–1400m范围完整包含两图当前数值，未为突出ELG局部变化而单独缩色标。

```powershell
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s02_sia/s02_arctic_region/n065_observed_march_mld_patterns.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'plot/S02_sia/s02_arctic_region/n065_observed_march_mld_patterns.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s02_sia/s02_arctic_region/check_observed_march_mld_patterns.py'
```

本轮全部为新正式路径，无覆盖已有output/figure。重新运行已有路径前须按archive/README归档。QC核查原文件哈希、月份和维序、经度排序、独立masked-array时间均值、NaN、支持标记、同点定位、量级统计、共同色标及实际导出尺寸；通过不代表两种方法等价或ELG绝对MLD物理有效性已验证。
