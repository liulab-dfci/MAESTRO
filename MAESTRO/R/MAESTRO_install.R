# install required R packages
library(devtools)

install_github("hms-dbmi/pagoda2", upgrade = "never")
install_github("SUwonglab/scABC@v0.1", upgrade = "never")
install.packages("stringi", repos="https://mirrors.tongji.edu.cn/CRAN/")
install_github("chenfeiwang/MAESTRO", upgrade = "never")
