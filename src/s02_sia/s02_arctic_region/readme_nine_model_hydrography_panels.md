# 九模式单模式N²、MLD空间变化与海岸分区序列

当前MLD单模式图已统一到作者选定的矩形口径，由n089计算准备、n090绘图负责，见readme_rectangular_mld_nine_models.md。原coastal_march_mld_timeseries_9models图已归档；n086当前只生成本页前两张N²与MLD空间图。沿岸科学数据保留用于比较，下文关于沿岸单模式图属于历史版本，不作为当前论文MLD序列。GFDL按作者确认的源第6年为分析第1年，不再据源1–5缺失否定端点诊断。

## 运行与路径

1. 计算准备：`src/s02_sia/s02_arctic_region/n085_nine_model_hydrography_panels.py`。
2. 绘图：`plot/S02_sia/s02_arctic_region/n086_nine_model_hydrography_panels.py`。
3. 核对：`src/s02_sia/s02_arctic_region/check_nine_model_hydrography_panels.py`。

保存数据在`output/s02_sia/s02_arctic_region/p25_cas_hydrography/n085_nine_model_hydrography_panels/`；正式PNG/PDF/SVG在`figure/s02_sia/s02_arctic_region/p25_cas_hydrography/n086_nine_model_hydrography_panels/`。绘图只读取前者，不读取原始温盐，不计算新平滑或分组。

## 三张图口径

- `labrador_irminger_initial_n2_9models_iap`：3×3，每个模式的初始March N²彩线叠加同一黑色IAP参照剖面。沿用截图的矩形LI分区50≤lat<65、65W≤lon<20W，深度0–1000m；IAP为2005–2024年March，不能解释为同期验证。八模式原始结果不改，原GFDL初始6–15年例外保留；CAS初始1–10年，对应日历200–209年。不做垂向外推或额外平滑，CAS最上层5m N²为缺测。
- `north_atlantic_march_mld_change_9models_ice`：3×3，填色为末期十年减初始十年的March MLD；蓝变浅、红加深，固定原色标±1100m。黑线为初始March MLD，200、500、1000、1500、2000m。紫虚线为初始、绿虚线为末期年平均15%海冰边缘，洋红实线为观测1988–2020年年平均15%冰缘。八模式MLD/冰缘/观测显示数组与原图完全保留；CAS原生MLD使用与原图一致的线性griddata显示到0.25°网格。CAS海冰读取最新完整1°保存产品，保留NaN；与旧八模式历史稀疏表恢复零值的显示流程差异已记元数据。观测保持原显示场，含原sigma0.8网格展示平滑，本次没有新增平滑。地图不作统计显著性推断。
- `coastal_march_mld_timeseries_9models`：3×3，蓝为Nordic–Iceland、橙为Labrador–Irminger海岸多边形区域绝对March MLD；灰虚线为CO₂。温盐TEOS-10 sigma0相对恰好10m增加0.03kg/m³的首次交点，跨层线性插值；无交点不纳入区域交点均值。原生面积加权后显示完整11年居中平均，首尾各5年缺测，不补齐。八模式原序列不改，既有GFDL对齐局限保留；CAS已逐月核对真实年标签。

三图均无顶部总标题或时期副标题；面板编号、模式名、坐标单位及图例保留。未做piControl漂移订正，未改变既有恢复分组。CAS原始计算及独立物理公式复核见`readme_cas_mld_stratification.md`与n083 QC。

## 维度与检查

本次读取保存产品，不重新读取NetCDF。N²表是model/depth的长表，9×100行；IAP100行。MLD是model/region/year长表，9×2×280行。地图MLD数组latitude×longitude=141×341；CAS海冰latitude×longitude=35×360；观测投影网格448×304。CAS原始NetCDF实际维序time/lev/lat/lon及坐标检查见CAS计算说明。

核对旧八模式剖面/序列、同一IAP、源文件SHA256、原地图所有NPZ数组逐项相同、CAS原生网格节点插值守恒及冰缘缺测、九模式面板、实际九个导出和PNG分辨率。PDF另核对单页并渲染预览，三张PNG均目视检查。旧的同用途八模式剖面及MLD序列图在新版确认后归档，未改动其他地区/方法对照图。
