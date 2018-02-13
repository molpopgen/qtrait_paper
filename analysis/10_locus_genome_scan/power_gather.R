library(readr)
library(dplyr)
library(tidyr)

mutate_alpha <- function(x)
{
    return(as.numeric(sub('[^0-9]+(.*)', '0.\\1', x)))
}

x = read_delim("nsig_per_window.txt.gz", delim = " ")

x = x %>% gather(stat, value, tajd05:max_abs_nSL001) %>% extract(stat, c('stat', 'alpha'), '([^0-9]+)([0-9]+)') %>% mutate(alpha=as.numeric(paste0('.', alpha)))

f = gzfile("nsig_per_window_tidied.txt.gz","wb")
write_delim(x,f,delim=" ")
close(f)
