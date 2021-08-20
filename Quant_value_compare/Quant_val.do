drop _all 

*Working directory (cd) must be manually set to where files are save
*cd ""

import delimited "Value_quant_data.csv" 

save "Data_all", replace

cor nurs_value nurs_count if year<2016 
cor seed_value seed_count if year<2016 
cor flower_value flower_count if year>1982 & year<2016  
