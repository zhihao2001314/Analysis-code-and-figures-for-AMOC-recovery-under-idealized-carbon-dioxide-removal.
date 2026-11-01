# n062：旧版矩形MLD区域示意图

用户要求调出当前已切回旧序列的分区图。检索到的已有独立地图均是后来的手绘版；因此本图为依据旧n04正式元数据补绘的示意图，不声称是历史旧图原件。

来源：`output/s02sia/s02_arctic_region/p06_mld/n04_na_region_timeseries/na_march_mld_region_metadata.json`；原代码`src/s01_amoc/n06_mldtime.py`第358–359行。

- 1，橙色虚线：Labrador–Irminger，50≤lat<65°N，−65≤lon≤−20°E。
- 2，蓝色实线：Nordic Seas，65≤lat≤80°N，−30≤lon≤20°E。

旧代码选择框内原生海洋格点中心，要求areacello有效且为正；区域均值还需该年有有效MLD。区域不随时间变化，也不按MLD>300/500m筛选。地图填色不表示MLD大小；灰色陆地不参与海洋统计。65°边界是否包含按上述严格不等式，不因地图边线宽度而改变。

src仅读取旧元数据，将矩形的四条等经纬边加密成每组405个顶点，以忠实表现Lambert投影下的弯曲边界，不改选区。`vertices.csv`实际维度region2×vertex_order405，字段longitude/latitude单位为度E/度N；`regions.csv`两行。结果保存`output/s02_sia/s02_arctic_region/p21_grouped_mld/n062_legacy_mld_regions/`。plot只读取这些output和已缓存Natural Earth110m底图，不读模式原始数据，不重算MLD，不调用下载。

图件保存`figure/s02_sia/s02_arctic_region/p21_grouped_mld/n062_legacy_mld_regions/legacy_rectangular_mld_sampling_regions.{png,pdf,svg}`；8.8×7.8英寸，PNG4400×3900、500dpi、不透明白底，PDF/SVG矢量。现有MLD序列、N²图和手绘区域产品均未改动。

科学可视化技能用于准确边界、颜色/线型/编号冗余编码与导出检查，软件参考沿用`readme_grouped_initial_n2_profiles.md`所列Scientific Agent Skills文献。不宣称特定期刊合规。

```powershell
& 'D:/Download/Anaconda/envs/viper/python.exe' 'src/s02_sia/s02_arctic_region/n062_legacy_mld_regions.py'
& 'D:/Download/Anaconda/envs/viper/python.exe' 'plot/S02_sia/s02_arctic_region/n062_legacy_mld_regions.py'
```

输出已存在时须先归档再重跑。本次为新路径，无覆盖。
