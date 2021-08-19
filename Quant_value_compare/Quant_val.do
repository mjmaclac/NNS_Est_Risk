drop _all 

cd "U:\NonNative\GitHub_docs\Quant_value_compare" 

import delimited "U:\NonNative\GitHub_docs\Quant_value_compare\Value_quant_data.csv" 

save "Data_all", replace

cor nurs_value nurs_count if year<2016 
cor seed_value seed_count if year<2016 
cor flower_value flower_count if year>1982 & year<2016  
