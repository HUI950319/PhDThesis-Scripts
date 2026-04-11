#!/usr/bin/env Rscript
#' SCENIC Pipeline R Wrapper
#' 
#' 该脚本提供了一个 R 函数来运行 pySCENIC 分析流程
#' 通过调用 conda 环境中的 Python 脚本执行完整的 GRN 分析
#'

#' 检测 conda 安装路径
#' 
#' @return conda 安装路径，如果未找到则返回 NULL
detect_conda_path <- function() {
  # 常见的 conda 安装路径
  possible_paths <- c(
    file.path(Sys.getenv("HOME"), "miniconda3"),
    file.path(Sys.getenv("HOME"), "anaconda3"),
    file.path(Sys.getenv("HOME"), "mambaforge"),
    file.path(Sys.getenv("HOME"), "miniforge3"),
    "/opt/conda",
    "/opt/miniconda3",
    "/opt/anaconda3"
  )
  
  for (path in possible_paths) {
    conda_sh <- file.path(path, "etc", "profile.d", "conda.sh")
    if (file.exists(conda_sh)) {
      return(path)
    }
  }
  
  # 尝试通过 which conda 查找
  conda_bin <- suppressWarnings(system("which conda", intern = TRUE, ignore.stderr = TRUE))
  if (length(conda_bin) > 0 && nchar(conda_bin) > 0) {
    # 从 conda 可执行文件路径推断 conda 根目录
    # 例如: /home/user/miniconda3/bin/conda -> /home/user/miniconda3
    conda_root <- dirname(dirname(conda_bin))
    if (file.exists(file.path(conda_root, "etc", "profile.d", "conda.sh"))) {
      return(conda_root)
    }
  }
  
  return(NULL)
}

#' 验证配置文件
#' 
#' @param config_file 配置文件路径
#' @return 如果验证通过返回 TRUE，否则抛出错误
validate_config <- function(config_file) {
  if (!file.exists(config_file)) {
    stop("配置文件不存在: ", config_file)
  }
  
  # 读取配置文件内容
  lines <- readLines(config_file, warn = FALSE)
  
  # 检查必需的节
  required_sections <- c("[global]", "[preprocess]", "[pyscenic]", "[aucell]")
  for (section in required_sections) {
    if (!any(grepl(section, lines, fixed = TRUE))) {
      stop("配置文件缺少必需的节: ", section)
    }
  }
  
  return(TRUE)
}


#' 运行 SCENIC Pipeline
#' 
#' 该函数通过调用 conda 环境中的 Python 脚本执行完整的 pySCENIC 分析流程
#' 
#' @param config_file 配置文件路径，默认为 "configs.txt"
#' @param conda_path conda 安装路径，NULL 表示自动检测
#' @param env_name conda 环境名称，默认为 "pyscenic-env"
#' @param skip_preprocess 是否跳过预处理步骤，默认 FALSE
#' @param skip_pyscenic 是否跳过 pySCENIC 步骤，默认 FALSE
#' @param skip_aucell 是否跳过 AUCell 步骤，默认 FALSE
#' @param verbose 是否显示详细输出，默认 TRUE
#' @param working_dir 工作目录，默认为 scenic_pipe 目录
#' 
#' @return 返回一个列表，包含执行状态和输出路径
#' 
#' @examples
#' \dontrun{
#' # 使用默认配置运行
#' result <- run_scenic()
#' 
#' # 指定配置文件
#' result <- run_scenic(config_file = "my_config.txt")
#' 
#' # 跳过预处理步骤
#' result <- run_scenic(skip_preprocess = TRUE)
#' }
#' 
#' @export
run_scenic <- function(
  config_file = "configs.txt",
  conda_path = NULL,
  env_name = "pyscenic-env",
  skip_preprocess = FALSE,
  skip_pyscenic = FALSE,
  skip_aucell = FALSE,
  verbose = TRUE,
  working_dir = NULL
) {


  # 记录开始时间
  start_time <- Sys.time()
  
  if (verbose) {
    cat("\n")
    cat("=" %s% rep("=", 59) %s% "\n")
    cat("  SCENIC Pipeline - R Wrapper\n")
    cat("=" %s% rep("=", 59) %s% "\n")
    cat("  开始时间:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
    cat("=" %s% rep("=", 59) %s% "\n\n")
  }
  
  # 确定工作目录
  if (is.null(working_dir)) {
    # 尝试找到 scenic_pipe 目录
    script_path <- tryCatch({
      # 如果是通过 source 运行的
      dirname(sys.frame(1)$ofile)
    }, error = function(e) {
      # 如果是直接运行的，使用当前目录
      getwd()
    })
    working_dir <- script_path
  }
  
  # 检查 main.py 是否存在
  main_py <- file.path(working_dir, "main.py")
  if (!file.exists(main_py)) {
    # 尝试在当前目录查找
    if (file.exists("main.py")) {
      working_dir <- getwd()
      main_py <- "main.py"
    } else {
      stop("找不到 main.py 文件，请确保在 scenic_pipe 目录下运行或指定正确的 working_dir")
    }
  }
  
  if (verbose) cat("工作目录:", working_dir, "\n")
  
  # 检测 conda 路径
  if (is.null(conda_path)) {
    conda_path <- detect_conda_path()
    if (is.null(conda_path)) {
      stop("无法自动检测 conda 路径，请手动指定 conda_path 参数")
    }
  }
  
  if (verbose) cat("Conda 路径:", conda_path, "\n")
  
  # 验证 conda.sh 存在

  conda_sh <- file.path(conda_path, "etc", "profile.d", "conda.sh")
  if (!file.exists(conda_sh)) {
    stop("conda.sh 不存在: ", conda_sh)
  }
  
  # 处理配置文件路径
  if (!file.path.is.absolute(config_file)) {
    config_file <- file.path(working_dir, config_file)
  }
  
  # 验证配置文件
  if (verbose) cat("配置文件:", config_file, "\n")
  validate_config(config_file)
  if (verbose) cat("配置文件验证通过!\n\n")
  
  # 构建命令参数
  cmd_args <- paste0("--config ", basename(config_file))
  
  if (skip_preprocess) {
    cmd_args <- paste(cmd_args, "--skip-preprocess")
  }
  if (skip_pyscenic) {
    cmd_args <- paste(cmd_args, "--skip-pyscenic")
  }
  if (skip_aucell) {
    cmd_args <- paste(cmd_args, "--skip-aucell")
  }
  
  # 构建完整的 shell 命令
  # 使用 bash -c 来执行，确保 source 命令可用
  shell_cmd <- sprintf(
    'cd "%s" && source "%s" && conda activate %s && python main.py %s',
    working_dir,
    conda_sh,
    env_name,
    cmd_args
  )
  
  if (verbose) {
    cat("执行命令:\n")
    cat("  ", shell_cmd, "\n\n")
  }
  
  # 执行命令
  if (verbose) cat("开始执行 SCENIC Pipeline...\n\n")
  
  exit_code <- system(paste("bash -c", shQuote(shell_cmd)), intern = FALSE)
  
  # 记录结束时间
  end_time <- Sys.time()
  elapsed_time <- difftime(end_time, start_time, units = "mins")
  
  if (verbose) {
    cat("\n")
    cat("=" %s% rep("=", 59) %s% "\n")
    cat("  执行完成\n")
    cat("  结束时间:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
    cat("  耗时:", round(as.numeric(elapsed_time), 2), "分钟\n")
    cat("  退出状态:", ifelse(exit_code == 0, "成功", "失败"), "\n")
    cat("=" %s% rep("=", 59) %s% "\n")
  }
  
  # 读取配置文件获取输出路径
  config_lines <- readLines(config_file, warn = FALSE)
  output_path_line <- grep("^output_path", config_lines, value = TRUE)
  if (length(output_path_line) > 0) {
    output_path <- gsub(".*=\\s*", "", output_path_line[1])
    output_path <- trimws(output_path)
    if (!file.path.is.absolute(output_path)) {
      output_path <- file.path(working_dir, output_path)
    }
  } else {
    output_path <- NULL
  }
  
  # 返回结果
  result <- list(
    success = (exit_code == 0),
    exit_code = exit_code,
    output_path = output_path,
    elapsed_time = elapsed_time,
    config_file = config_file,
    working_dir = working_dir
  )
  
  class(result) <- c("scenic_result", "list")
  
  return(invisible(result))
}


#' 辅助函数：字符串拼接
#' @keywords internal
`%s%` <- function(a, b) paste0(a, b, collapse = "")


#' 检查文件路径是否为绝对路径
#' @keywords internal
file.path.is.absolute <- function(path) {
  # Linux/Mac: 以 / 开头
  # Windows: 以盘符开头 (如 C:/)
  grepl("^(/|[A-Za-z]:)", path)
}


#' 打印 scenic_result 对象
#' @export
print.scenic_result <- function(x, ...) {
  cat("SCENIC Pipeline 执行结果\n")
  cat("--------------------------\n")
  cat("状态:", ifelse(x$success, "成功", "失败"), "\n")
  cat("退出码:", x$exit_code, "\n")
  cat("耗时:", round(as.numeric(x$elapsed_time), 2), "分钟\n")
  cat("配置文件:", x$config_file, "\n")
  cat("工作目录:", x$working_dir, "\n")
  if (!is.null(x$output_path)) {
    cat("输出路径:", x$output_path, "\n")
  }
  invisible(x)
}


# 如果直接运行此脚本，则执行示例
if (!interactive() && sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) == 0) {
    cat("用法: Rscript run_scenic.R <config_file>\n")
    cat("示例: Rscript run_scenic.R configs.txt\n")
  } else {
    result <- run_scenic(config_file = args[1])
    if (!result$success) {
      quit(status = 1)
    }
  }
}

