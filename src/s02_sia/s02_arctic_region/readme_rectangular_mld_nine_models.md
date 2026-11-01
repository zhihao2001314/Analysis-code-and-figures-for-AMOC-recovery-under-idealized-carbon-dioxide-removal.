# 九模式混合层序列统一为选定矩形口径

作者指出n078分组图与今天n086九模式单图的MLD数值不同。本次查明：n078中排读取n061旧矩形；n086单图读取海岸多边形（Nordic包含冰岛盆地），并且前者最近有效10m原生层参考密度、后者恰好10m插值参考密度。两者均为TEOS10 sigma0增密0.03kg/m³首次交点、无交点不纳入区域均值、原生面积加权及完整11年居中平滑。不能把差异全部归因于分区，也不是新增piControl订正或CAS数据补充造成。

沿用作者此前明确选择的矩形，本次把九模式单图统一为：LI50≤lat<65、65W≤lon≤20W；NS65≤lat≤80、30W≤lon≤20E。八模式原始序列与选定n078图LI曲线逐点一致；CAS取n083已经计算的矩形 `mld_nearest_crossing_m` 及其已存11年曲线，保持参考层处理与原八模式一致。GFDL严格按作者确认的源6年为分析1年，6–15为初始、源145为分析峰值140；没有再次移5年或把源1–5误作缺失。

计算准备 `n089_rectangular_mld_nine_models.py`；绘图 `plot/S02_sia/s02_arctic_region/n090_rectangular_mld_nine_models.py`。正式数据 `output/s02_sia/s02_arctic_region/p25_cas_hydrography/n089_rectangular_mld_nine_models/`，正式图 `figure/s02_sia/s02_arctic_region/p25_cas_hydrography/n090_rectangular_mld_nine_models/rectangular_march_mld_timeseries_9models.{png,pdf,svg}`。保存路径在绘图参数区明确列出并打印。

例如GFDL LI初始/末期十年均值：矩形293.799/594.344m，沿岸366.922/765.476m；NS为矩形215.833/324.575m，沿岸289.243/448.629m。完整九模式、两海区对照保存 `rectangular_vs_coastal.csv`。这些是未平滑十年均值，图显示完整11年平滑，不能直接把曲线可见端点当十年平均。

4480条旧矩形数值保留；其中2240条LI与n078选定分组图源表逐点一致，NS来自同一n061矩形产品。18组初期、末期和变化均与当前Table1核对通过。没有重读原始NetCDF，没有改MLD空间pattern、初始N²图、原分组图或恢复成员。CSV实际维度为model9×region2×simulation_year280，具名字段读取；首尾各5个平滑年份缺测保留。

旧沿岸单图3格式及旧总清单已按archive_record.json记录归档，逐文件核验源消失、目标存在和SHA256；其他两张N²/MLD空间图共6个导出文件哈希不变。n086绘图当前只生成前两张图，单模式MLD统一由n090负责。旧科学沿岸数据保留为区域敏感性比较来源，不作物理删除。

绘图保持原九面板、0–900m共享坐标、蓝NS/橙LI/灰CO₂布局；图例改为真实矩形海区名称，无总标题。PNG5760×4500、450dpi，PDF与SVG导出，实际PNG与独立PDF栅格均已目视。Scientific-visualization技能用于相同口径、缺测、图例及实际导出检查，不作期刊合规声明。
