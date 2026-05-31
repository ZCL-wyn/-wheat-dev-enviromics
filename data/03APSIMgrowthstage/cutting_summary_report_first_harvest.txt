=== 基于物候阶段的数据切割摘要报告（以第一个HarvestRipe为准） ===

生成时间： 2026-02-08 11:36:19.596418 

1. 数据基本信息：
   对齐后原始数据环境数： 8 
   对齐后原始数据总行数： 2249 
   切割后数据环境数： 2 
   切割后数据总行数： 498 

2. HarvestRipe情况：
   有多个HarvestRipe的环境数： 2 
   有1个HarvestRipe的环境数： 0 
   有2个HarvestRipe的环境数： 2 
   有3+个HarvestRipe的环境数： 0 

3. 物候边界确定：
   成功确定Sowing日期的环境数： 8 
   成功确定HarvestRipe日期的环境数： 2 
   有效环境数（同时有两者）： 2 

4. 切割效果：
   平均切割前行数： 282 
   平均切割后行数： 249 
   平均删除行数： 33  ( 11.7 %)
   平均切割前天数： 282 
   平均切割后天数： 249 
   平均删除天数： 33  ( 11.7 %)

5. 切割规则验证：
   第一个阶段是Sowing的环境数： 2 / 2 
   最后一个阶段是HarvestRipe的环境数： 2 / 2 
   有多个HarvestRipe的环境中正确选择第一个的比例： 100 %

6. 物候完整性：
   有Sowing阶段的环境数： 2 / 2 
   有HarvestRipe阶段的环境数： 2 / 2 
   有Anthesis阶段的环境数： 2 / 2 
   平均关键阶段完整率： 100 %

7. 生育期长度：
   平均生育期长度： 249 天
   最短生育期： 242 天
   最长生育期： 256 天
   平均播种到开花： 192.5 天
   平均开花到成熟： 55.5 天

8. 输出文件清单：
   - 切割后的气象数据：weather_data_cut_by_phenology_first_harvest.csv
   - 物候边界信息（含计数）：phenology_boundaries_with_counts.csv
   - HarvestRipe计数信息：harvestripe_counts.csv
   - 切割效果统计：cutting_effect_statistics.csv
   - 物候完整性统计：phenology_completeness.csv
   - 生育期长度统计：growth_period_statistics.csv
   - 多个HarvestRipe环境验证详情：multiple_harvestripe_validation.csv

✅ 基于物候阶段的数据切割已完成！
⚠️ 特别注意：对于有多个HarvestRipe的环境，已以第一个HarvestRipe为准进行切割。
