> 当前组合图已改为n093/n094七面板紧凑版：两区九模式N2/MLD分别叠加，UK移入Non-recovery，地图3/3/3，含组平均海冰边缘。见 `readme_compact_regional_hydrography.md`。下文为已归档旧九面板版本。

# 加入CAS的九模式层结与混合层组合图

## 当前图件与代码

- 计算：`src/s02_sia/s02_arctic_region/n091_grouped_hydrography_with_cas.py`。
- 输出：`output/s02_sia/s02_arctic_region/p23_stratification_mld_composite/n091_grouped_hydrography_with_cas/`。
- 绘图：`plot/S02_sia/s02_arctic_region/n092_grouped_hydrography_with_cas.py`，只读取output。
- 图件：`F:/sea/Antic sea ice constraint/figure/s02_sia/s02_arctic_region/p23_stratification_mld_composite/n092_grouped_hydrography_with_cas/labrador_irminger_9panels_with_cas_rectangles.png`，以及同名PDF/SVG。

## CAS归组依据

本图沿用初始LI层结表现与末期AMOC恢复程度分组，不是Fig1的纯滞回形态分组。
CAS初始LI剖面相对IAP的RMSD为0.2671×10⁻⁵s⁻²，深度15–995m；其末期/初始十年AMOC比值为26.5°N约0.8984、40°N约1.0186，均满足既有0.8恢复判据。
所以CAS加入中间基本恢复组。这里IAP是现代气候参照，与模式初始气候并非同时期观测，不将RMSD解释为严格统计优劣检验。

| 列 | 成员 | 图内标签 |
|---|---|---|
| 左 | CESM2、CNRM-ESM2-1 | Non-recovery |
| 中 | NorESM2-LM、GFDL-ESM4、CAS-ESM2-0 | Substantial recovery |
| 右 | ACCESS-ESM1-5、CanESM5、MIROC-ES2L、UKESM1-0-LL | Larger stratification bias |

原八模式保留既有分组。右列包含不同恢复结果，因此删除原图将其笼统标成Substantial recovery的错误标签。未迁移其他图件、正文统计或全项目分组注册表。

## 数据与区域口径

- 剖面与时间序列均为LI矩形区域，N²采用50≤lat<65、−65≤lon<−20；MLD遵循已保存旧矩形的外边界取格口径。
- 原八条N²和MLD曲线保持数据不变，CAS取既有n083矩形初始三月剖面，MLD取n089统一矩形序列。
- n089沿用最近有效原生10m层参考密度、sigma0增加0.03kg/m³的已解析首次阈值穿越，排除无穿越柱；没有改用海岸多边形或原始mlotst。
- 原样保留三月MLD完整11年居中平均和两端缺测；组均值为固定成员等权平均。IAP三月2005–2024均值和年际样本标准差文件逐字节不变。
- GFDL源第6年为分析第1年；初始源6–15，不再使用“缺失初始五年”的说明。
- 地图为完整北大西洋逐格点场，不做矩形区域平均。复用原Nor/GFDL空间场并加入n085保存的CAS地图，按固定三成员共同有效掩膜等权平均。初始、末期和差值满足平均差值恒等式。
- 地图和区域时间序列空间支持与既有格点MLD参考深度处理不同；没有宣称地图可直接重建中排的区域均值，也没有新增插值或漂移订正。

## 维度及验证

不读取新的NetCDF/MAT。现成CSV：模型9×深度100，模型9×年份280；地图数组每组141×341，初始/末期/差值分别保存。
中列三行的组均值均已更新。两侧地图数组逐元素不变，27条实际绘制曲线和源数据哈希已核对，标签与曲线未重叠。
PNG5670×5760、450dpi，PDF单页且无Type3字体，PNG及独立PDF栅格均已目视检查。旧八模式n078图已按原相对路径归档，当前正式图为n092。
