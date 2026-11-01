# Nordic Seas九模式与IAP的初始March N²

计算准备`n087_nordic_initial_n2_nine_models.py`只提取已保存的八模式初始剖面、n083 CAS矩形初始剖面及IAP，未重新读取原始温盐。绘图`plot/S02_sia/s02_arctic_region/n088_nordic_initial_n2_nine_models.py`仅读取output中的这套提取结果。

## 分区与统计

矩形Nordic Seas：65≤latitude<80°N，30°W≤longitude<20°E，只统计有效海洋格点。与上张Labrador–Irminger矩形图（50≤latitude<65°N，65°W≤longitude<20°W）对应，不采用海岸多边形。此前九模式MLD时间序列用的是海岸多边形，两种分区不混称。

九个模式均为原始初始十年March N²剖面，继承既有GFDL第6–15年的例外；CAS为实验1–10年（日历200–209）。每格彩线为该模式，黑线为完全相同的IAP 2005–2024年March平均剖面。IAP是现代气候参照，与理想试验初始期并非同期。无新piControl订正、垂向平滑或外推，缺测保留。

输入长表维度为model9×depth100，IAP为depth100；深度中点5、15、…、995m。沿用LI图横轴−0.6至9.2（单位10⁻⁵s⁻²）及0–1000m深度轴，便于两海区直接比较，已检查有效值不会被横轴裁切。CAS原始NetCDF实际维序time/lev/lat/lon及GSW计算见`readme_cas_mld_stratification.md`，此脚本不读取NetCDF。

## 正式输出

数据：`output/s02_sia/s02_arctic_region/p25_cas_hydrography/n087_nordic_initial_n2_nine_models/`，含提取表、元数据、源SHA256及图件核对。

图件：`figure/s02_sia/s02_arctic_region/p25_cas_hydrography/n088_nordic_initial_n2_nine_models/nordic_seas_initial_n2_9models_iap.{png,pdf,svg}`。

3×3、CAS在第九格，无顶部总标题或时期副标题；完整图片路径在绘图参数区明示，保存后打印。旧八模式加单独IAP的Nordic图在新版核验后归档，保留源相对路径及SHA256核对记录。
