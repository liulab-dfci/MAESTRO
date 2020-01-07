# install required R packages
library(devtools)

# install_github("hms-dbmi/pagoda2", upgrade = "never")
install_github("SUwonglab/scABC@v0.1", upgrade = "never")
install_version(package = "Seurat", version = "3.0.2", repos="https://mirrors.tongji.edu.cn/CRAN/")
install_github("chenfeiwang/MAESTRO", upgrade = "never")
