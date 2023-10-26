#ifndef UTILS_LOGGER
#define UTILS_LOGGER
#include <spdlog/spdlog.h>

#include <string>

class Logger {
public:
  Logger(std::string logger_name) {
    logger_ = spdlog::get(logger_name);
    logger_->set_pattern("[%H:%M:%S.%e] [%^%l%$] %v");
  };
  ~Logger();
  template <typename... Args> void info(Args &&...args) { logger_->info(std::forward<Args>(args)...); }
  template <typename... Args> void warn(Args &&...args) { logger_->warn(std::forward<Args>(args)...); }
  template <typename... Args> void error(Args &&...args) { logger_->error(std::forward<Args>(args)...); }
  template <typename... Args> void debug(Args &&...args) { logger_->debug(std::forward<Args>(args)...); }

  void set_level(spdlog::level::level_enum level) { logger_->set_level(level); }

private:
  std::shared_ptr<spdlog::logger> logger_;
};
#endif