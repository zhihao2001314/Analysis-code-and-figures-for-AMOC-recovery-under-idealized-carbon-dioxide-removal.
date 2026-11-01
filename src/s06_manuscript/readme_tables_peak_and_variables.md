# 九模式峰值实际状态、末期变化及变量时段表

## 当前作者确认口径 2026-10-29

**GFDL温盐源第6年就是分析第1年。** 初始十年源6–15，升碳源6–145映射分析1–140，峰值分析140对应源145，末期源271–280。不是缺少分析初始五年；此前物理第140年提取及缺测解释均已被作者纠正，不再用于交付。

当前入口为 `n07_gfdl_start_year6.py` → `n05_tables_word.py`（默认读取n07）。正式数据 `output/s06_manuscript/n07_gfdl_start_year6/`；最终Word为 `docs/table1_cmip6_peak_and_residual_gfdl_start6.docx` 与 `docs/table2_variables_and_periods_gfdl_start6.docx`。GFDL峰值MLD为LI48.2858m、NS118.1647m；残余变化LI+300.5452m、NS+108.7419m，表内取整。原模式系列未改，其他八模式表值不变。作者此前打开的旧Word保留，未强行关闭或覆盖。

以下为历史实现与QA记录；涉及GFDL缺测或物理140的内容均已由本节取代。

计算入口 `n04_tables_peak_and_variables.py`，Word排版入口 `n05_tables_word.py`。只读导出核查使用 `check_tables_word_render.ps1`。最终文件为 `docs/table1_cmip6_peak_and_residual.docx` 和 `docs/table2_variables_and_periods.docx`；计算及QA存于 `output/s06_manuscript/n04_tables_peak_and_variables/`。

## 口径

- 用户明确要求峰值为CO₂最高的实验第140年，展示变量实际值，既不是131–140年平均，也不是相对初始的变化。T和SIA为当年年平均；MLD为当年3月。
- 末期栏保留残余变化：271–280年平均减1–10年平均。原八模式T和MLD结果保持；未新增piControl漂移订正。
- CAS气温直接读取标准2m月tas，原生球面格点面积加权、实际月日数年平均；不采用旧最低模式层代理。原生维序 `(time,lat,lon)`，128纬×256经。上升初始200–209年、峰值339年；下降末期470–479年。
- SIA采用九模式同一未订正年序列的 `raw_sia_million_km2`，明确为45–90°N海冰面积，不是海冰范围。旧端点空间产品与新序列的差异保存在CSV，不将区域分区之和误作全部北极面积。原八模式末期数值与旧表一致，CAS末期更新为约−2.03百万km²。
- MLD延续原表矩形、原生面积平均、最近有效10m层作为位密度参考，密度阈值0.03kg/m³、首次交点、无交点排除。CAS使用n083保存的同口径 `mld_nearest_crossing_m`。这与新沿岸图恰好10m插值参考/多边形平均不同，不能静默混合。
- LI：50≤lat<65、65W≤lon≤20W；NS：65≤lat≤80、30W≤lon≤20E。CAS原生结果维序及温盐匹配详见s02的n083说明。
- 旧GFDL矩形MLD真实年份仍未核实，且最早温盐从物理第6年起；表中峰值和残余MLD均不填旧索引值。未借用mlotst、沿岸或别的时期补值。
- Table2第一页列核心分析的变量、模式范围、实验相对时段及两个观测产品。第二页明确列tas源文件日历标签，并另列CAS温盐与海冰时轴例外；没有把tas日历推广到全部海洋变量。实际论文若增删其他探索性图，应对应更新数据清单。

## Word与QA

沿用旧DOCX的Times New Roman黑白三线表，原生Word单元格、固定列宽673pt、取消自动调整、禁止拆分行，表头重复。扩展列后改为Letter横向；复制到稿件时应保留源格式并使用横向页面。不同目标页宽或合并格式可能改变布局，不作跨文档无条件不变的承诺。

内置render_docx因未安装LibreOffice失败；未下载软件。改用已安装Microsoft Word后台只读导出PDF，再以捆绑Poppler渲染。Table1单页、Table2两页，所有最终页已目视核查，无蓝色标题边线、裁切、重叠或孤立表注。旧试排12文件已归档并逐文件核对源消失、目标存在及SHA256，记录于 `archive/202610291704_word_tables_render_refinement/`。


## GFDL第140年补充核查

2026-10-29：原n04序列QA明确GFDL升碳源年6–145、降碳源年141–280。原代码把升碳前140个March逐个标成1–140，因此物理140年对应旧索引135，可以报告峰值LI72.2036m、NS76.5530m。初始物理1–5年仍缺，严格271–280减1–10的ΔMLD继续缺测，不偷偷用6–15代替。旧全序列及图件未修改。

新核查计算n06_gfdl_mld_peak_time_audit.py，正式结果在output/s06_manuscript/n06_gfdl_mld_peak_time_audit/。当前可交付Word以_gfdl140.docx结尾，原Word被用户打开锁定所以没有覆盖。第一次归档只成功复制旧表1，源与目标SHA相同，源仍在；partial_archive_status.json记录此状态，没有反复移动或关闭用户Word。新Word已只读导出PDF并检查全部3页。
